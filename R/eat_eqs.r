
eat.from.cluster = function(df, cluster=df$cluster[1], eating.funs = list(eat.single.front,eat.single.back,eat.double.front), mat=NULL, repeated=TRUE, verbose=TRUE) {
  restore.point("eat.from.cluster")

  num.funs = length(eating.funs)
  fun.ind = last.fun = 1

  while(TRUE) {
    fun = eating.funs[[fun.ind]]
    while (TRUE) {
      res = fun(df=df,cluster=cluster, mat=mat, verbose=verbose)
      if (res$changed==0) {
        if (!repeated) last.fun = fun.ind
        break
      }
      df = res$df; mat = res$remaining.mat; cluster = res$remaining.cluster
      last.fun = fun.ind
      if (!repeated) break
    }
    fun.ind = ((fun.ind +1 -1) %% num.funs)+1
    if (fun.ind ==last.fun) break
  }
  df
}

old.eat.from.cluster.df = function(df, cluster=df$cluster[1], var.li = df$vars, syms=df$var, skip.eat = FALSE, only.calls=NULL, repeated = TRUE) {

  nr = length(var.li); nc = length(syms)
  mat = matrix(0,nr,nc)
  colnames(mat) = syms
  for (i in seq_along(var.li)) {
    mat[i,var.li[[i]]] = 1
  }
  mat

  if (!skip.eat) {
    calls = create.eating.calls()
    if (!is.null(only.calls)) {
      calls = calls[only.calls]
    }
    sum.changed = 1
    while(sum.changed>0) {
      sum.changed = 0
      changed = 0
      for (call in calls) {
        eval(call)
        sum.changed = sum.changed + changed
      }
      if (!repeated) break
    }
  }

  # put remaining rows into a cluster
  rows = which(df$cluster==0)
  if (length(rows)>0) {
    remain.syms = setdiff(syms, df$var)

    df$cluster[rows] = max(df$cluster)+1
    df$cluster.size[rows] = length(rows)
    df$var[rows] = remain.syms
    df$level[rows] = max(df$level)+1
  }

  # reindex negative clusters and levels
  # put -1 at end, before it -2 and so on...
  num.cluster = length(unique(df$cluster))
  num.level = length(unique(df$level))

  rows = which(df$cluster < 0)
  new.cluster = num.cluster + df$cluster[rows] +1
  new.level = num.level + df$level[rows] +1
  df$cluster[rows] = new.cluster
  df$level[rows] = new.level
  df = df[order(df$cluster),]
  df
}


make.eating.matrix = function(df, var.li=df$vars, syms=df$var) {
  nr = NROW(df); nc=NROW(df)
  mat = matrix(0,nr,nc)
  colnames(mat) = df$var
  for (i in seq_along(var.li)) {
    mat[i,intersect(var.li[[i]], syms)] = 1
  }
  mat
}

# search for equations that have a single variable
# remove them from the cluster and put them in front
eat.single.front = function(df, cluster, mat=NULL, repeated=FALSE, verbose=TRUE) {
  restore.point("eat.single.front")

  crows = which(df$cluster==cluster)
  if (is.null(mat)) {
    mat = make.eating.matrix(df)
    mat = mat[crows,crows]
  }
  syms = colnames(mat)
  # equations that have a single variable
  rows = which(rowSums(mat) == 1)
  if (length(rows)==0)
    return(list(changed=0, df=df, remaining.cluster=cluster, extracted.clusters=NULL, remaining.mat=mat))

  # variables that correspond to those equations
  cols = max.col(mat[rows,,drop=FALSE])

  esyms = syms[cols]

  # we have variables that are determined by two or more equations
  # decide which equations to take and which to drop
  if (any(duplicated(esyms))) {
    drop.inds = NULL

    dupl.syms = unique(esyms[duplicated(esyms)])

    txt = paste0("The variables ", paste0(dupl.syms, collapse=", ")," are determined by more than one equation:")

    sym = dupl.syms[1]
    for (sym in dupl.syms) {
      inds = which(esyms==sym)
      arows = crows[rows[inds]]
      aeqs = df$eq_[arows]

      # select the equation we want to solve sym for
      ret = suggest.eq.for.var(eqs=aeqs, var=sym)
      # Index of rows that are not taken
      drop.inds = c(drop.inds, inds[-ret$eq.ind])

      txt = c(txt, paste0("  ",sym,":"), paste0("    - ",sapply(aeqs, deparse1),collapse="\n"))
    }
    txt = paste0(txt, collapse="\n")
    if (verbose)
      cat(txt)

    # only keep selected equations
    rows = rows[-drop.inds]
    cols = max.col(mat[rows,,drop=FALSE])
    esyms = syms[cols]


  }

  erows = crows[rows]
  rrows = setdiff(crows, erows)
  rsyms = setdiff(syms, esyms)
  ne = length(rows)
  remaining.mat=mat[-rows,rsyms, drop=FALSE]

  df$cluster[erows] = cluster - (1:ne) / (ne+1)
  df$level[erows] = df$level[erows[1]] - 0.5
  df$cluster = rank(df$cluster,ties.method = "min")
  df$level = rank(df$level,ties.method = "min")
  df$cluster.size[erows] = 1
  df$cluster.size[rrows] = df$cluster.size[rrows[1]]-ne
  df$var[erows] = esyms
  df$var[rrows] = rsyms

  remaining.cluster = df$cluster[rrows[1]]
  extracted.clusters = unique(df$cluster[erows])

  df = df[order(df$cluster),]
  changed = length(rows)

  if (verbose) {
    cat("\neat.single.front:")
    cat("\nvars: ", paste0(esyms, collapse=", "))
    cat("\nrows: ", paste0(erows, collapse=", "),"\n")
  }

  if (repeated) {
    res = eat.single.front(df=df, cluster = remaining.cluster, mat=remaining.mat)
    remaining.cluster = res$remaining.cluster
    remaining.mat = res$remaining.mat
  }

  return(list(changed=changed, df=df, remaining.cluster=remaining.cluster, extracted.clusters=extracted.clusters, remaining.mat=remaining.mat))


}


eat.double.front = function(df, cluster, mat=NULL, verbose=TRUE) {
  restore.point("eat.double.front")

  crows = which(df$cluster==cluster)
  if (is.null(mat)) {
    mat = make.eating.matrix(df[crows,])
  }
  syms = colnames(mat)

  # find equations that contain the same two variables
  vars.id = paste0(as.data.frame(t(mat)), sep="|")

  rows = which((duplicated(vars.id) | duplicated(vars.id, fromLast = TRUE)) & rowSums(mat) == 2)

  if (length(rows)==0)
    return(list(changed=0, df=df, remaining.cluster=cluster, extracted.clusters=NULL, remaining.mat=mat))

  changed = length(rows)
  first = which(duplicated(vars.id, fromLast = TRUE))

  erows = crows[rows]
  df$cluster.size[erows] = 2
  df$level[erows] = df$level[crows[1]]-0.5

  counter = 0
  esyms = NULL
  for (row in first) {
    counter = counter+1
    srows = crows[rows[vars.id[rows]==vars.id[row]]]
    cur.syms = syms[which(mat[row,]==1)]
    esyms = c(esyms, cur.syms)
    df$cluster[srows] = cluster+ counter / (length(rows)+1)
    df$var[srows] = cur.syms
  }

  if (verbose) {
    cat("\neat.double.front: ")
    cat("\nvars: ", paste0(esyms, collapse=", "))
    cat("\nrows: ", paste0(erows, collapse=", "),"\n")
  }

  rrows = setdiff(crows, erows)
  rsyms = setdiff(syms, esyms)
  ne = length(rows)
  remaining.mat=mat[-rows,rsyms]

  df$level[erows] = df$level[erows[1]] - 0.5
  df$cluster = rank(df$cluster,ties.method = "min")
  df$level = rank(df$level,ties.method = "min")
  df$cluster.size[rrows] = df$cluster.size[rrows[1]]-ne
  df$var[rrows] = rsyms

  remaining.cluster = df$cluster[rrows[1]]
  extracted.clusters = unique(df$cluster[erows])

  df = df[order(df$cluster),]
  changed = length(rows)

  return(list(changed=changed, df=df, remaining.cluster=remaining.cluster, extracted.clusters=extracted.clusters, remaining.mat=remaining.mat))


}


# Try to remove equations at the back of a cluster
# search for variables that appear in only one equation
# this equation must then define this variable
# we can put the equation behind the larger cluster
eat.single.back = function(df, cluster, mat=NULL, repeated=FALSE, verbose=TRUE) {
  restore.point("eat.single.back")

  crows = which(df$cluster==cluster)
  if (is.null(mat)) {
    mat = make.eating.matrix(df[crows,])
  }
  syms = colnames(mat)

  # Compute for each symbol the number of equations
  # it is still contained in
  num.eqs = colSums(mat)

  # Find the symbol columns that are only in a single
  # equation
  cols = which(num.eqs==1)
  if (length(cols)==0) {
    return(list(changed=0, df=df, remaining.cluster=cluster, extracted.clusters=NULL, remaining.mat=mat))
  }
  # Identify the equations that have the single columns
  rows = which(rowSums(mat[,cols, drop=FALSE]) >= 1)


  # Error two or more variables are uniquely defined by the same equation:
  # we cannot solve that equation
  if (length(rows)<length(cols)) {
    cols.rows = sapply(cols, function(col) which(mat[,col]==1))
    dup.rows = unique(cols.rows[duplicated(cols.rows) | duplicated(cols.rows,fromLast = TRUE)])

    txt = "In eat.single back:"
    for (dup.row in dup.rows) {
      dup.vars = names(cols.rows)[cols.rows==dup.row]
      dup.eq = crows[dup.row]
      txt = c(txt, paste0("  - The variables ", paste0(dup.vars,collapse=", "), " are only defined in the single equation ", dup.eq, ":\n    ",deparse1(df$eq_[[dup.eq]])))
      if (!identical(df$eq_[[dup.eq]],df$org_[[dup.eq]])) {
        txt = c(txt, paste0("  Original form was:\n    ",deparse1(df$org_[[dup.eq]]) ))
      }
    }
    if (verbose) {
      txt = c(txt,"Skip eat.single.back.")
      cat(paste0(txt,collapse="\n"))
    }
    return(list(changed=0, df=df, remaining.cluster=cluster, extracted.clusters=NULL, remaining.mat=mat))

  }

  ne = length(rows)
  erows = crows[rows]
  later.level = df$level[crows[1]]+0.5
  if (verbose) {
    cat("\neat.single.back: ")
    cat("\nvars: ", paste0(names(cols), collapse=", "))
    cat("\nrows: ", paste0(erows, collapse=", "),"\n")
  }
  df$level[erows] = later.level
  df$cluster[erows] = cluster + (1:ne) / (ne+1)
  df$cluster.size[erows] = 1
  for (col in cols) {
    crow = crows[which(mat[,col] >= 1)]
    df$var[crow] = syms[col]
  }

  remaining.mat = mat[-rows, -cols]

  esyms = syms[cols]
  rrows = setdiff(crows, erows)
  rsyms = setdiff(syms, esyms)
  ne = length(rows)

  df$cluster = rank(df$cluster,ties.method = "min")
  df$level = rank(df$level,ties.method = "min")
  df$cluster.size[rrows] = df$cluster.size[rrows[1]]-ne
  df$var[rrows] = rsyms

  remaining.cluster = df$cluster[rrows[1]]
  extracted.clusters = unique(df$cluster[erows])

  df = df[order(df$cluster),]
  changed = ne

  if (repeated) {
    res = slow.eat.single.back(df=df, cluster = remaining.cluster, mat=remaining.mat)
    remaining.cluster = res$remaining.cluster
    remaining.mat = res$remaining.mat
    extracted.clusters = NA
  }

  return(list(changed=changed, df=df, remaining.cluster=remaining.cluster, extracted.clusters=extracted.clusters, remaining.mat=remaining.mat))


}


# Just a helper function to make the code size
# in equation.system.cluster look smaller. This makes debugging of it easier
old.create.eating.calls = function(env = parent.frame()) {
  calls = list(

    # search for equations that have a single variable
    # remove them from the cluster and put them in front
    eat.single.front = quote({
      changed = 0
      while(TRUE) {
        rows = which(df$active & rowSums(mat) == 1)
        if (length(rows)==0) break
        cols = max.col(mat[rows,,drop=FALSE])
        cur.syms = syms[cols]
        if (any(duplicated(cur.syms))) {
          dupl.syms = unique(cur.syms[duplicated(cur.syms)])
          warning(paste0("The variables "), paste0(dupl.syms, collapse=", ")," are determined by more than one equation!")
        }

        df$active[rows] = FALSE
        df$cluster[rows] = max(df$cluster)+seq_along(rows)
        df$cluster.size[rows] = 1
        df$level[rows] = max(df$level+1)
        df$var[rows] = cur.syms
        mat[,cur.syms] = 0
        changed = changed + length(rows)

        if (verbose) {
          cat("\neat.single.front:")
          cat("\nvars: ", paste0(cur.syms, collapse=", "))
          cat("\nrows: ", paste0(rows, collapse=", "),"\n")
        }
      }
      changed
    }),

    # search for systems of equations that only have two variables
    # (the variables must be the same variables)
    # put those equation systems in the front
    eat.double.front = quote ({
      changed = 0
      while(TRUE) {
        rows = which((duplicated(df$vars.id) | duplicated(df$vars.id, fromLast = TRUE)) & rowSums(mat) == 2)
        if (length(rows)==0) break
        changed = changed + length(rows)
        first = which(duplicated(df$vars.id, fromLast = TRUE))

        df$active[rows] = FALSE
        df$cluster.size[rows] = 2
        df$level[rows] = max(df$level)+1

        for (row in first) {
          crows = which(df$vars.id==df$vars.id[[row]])
          cur.syms = syms[which(mat[row,]==1)]
          df$cluster[crows] = max(df$cluster)+1
          df$var[crows] = cur.syms
        }
        vars = unique(unlist(var.li[rows]))
        mat[,vars] = 0

        level = level+1

        if (verbose) {
          cat("\neat.double.front: ")
          cat("\nvars: ", paste0(vars, collapse=", "))
          cat("\nrows: ", paste0(rows, collapse=", "),"\n")
        }

      }
      changed
    }),

    # Try to remove equations at the back of a cluster
    # search for variables that appear in only one equation
    # this equation must then define this variable
    # we can put the equation behind the larger cluster
    eat.single.back = quote({
      changed = 0
      while(TRUE) {
        # Compute for each symbol the number of equations
        # it is still contained in
        num.eqs = colSums(mat)
        # Find the symbol columns that are only in a single
        # equation
        cols = which(num.eqs==1)
        if (length(cols)==0) break
        # Identify the equations that have a single column
        rows = which(rowSums(mat[,cols, drop=FALSE]) >= 1)
        later.level = min(df$level)-1
        if (verbose) {
          cat("\neat.single.back: ", later.level, "")
          cat("\nvars: ", paste0(names(cols), collapse=", "))
          cat("\nrows: ", paste0(rows, collapse=", "),"\n")
        }
        df$level[rows] = later.level
        df$active[rows] = FALSE
        df$cluster[rows] = min(df$cluster)-seq_along(rows)
        df$cluster.size[rows] = 1
        for (col in cols) {
          row = which(mat[,col] >= 1)
          df$var[row] = syms[col]
        }

        mat[rows,] = 0
        changed = changed + length(rows)
      }
      changed
    }),
    # Try to remove equations from back
    eat.double.back = quote({
      changed = 0
      while(TRUE) {
        # Compute for each symbol the number of equations
        # it is still contained in
        num.eqs = colSums(mat)
        # Find the symbol columns that are in exactly two equations
        cols = which(num.eqs==2)
        if (length(cols)<2) break

        col.rows = lapply(cols, function(col) {
          rows = which(mat[,col] >= 1)
          sort(rows)
        })
        col.rows = do.call(rbind,col.rows)
        dupl = which(duplicated(col.rows, fromLast = TRUE))
        if (length(dupl)==0) break

        level = min(df$level)-1
        for (ind in dupl) {
          inds = which(col.rows[,1] == col.rows[ind,1] & col.rows[,2] == col.rows[ind,2])
          ccols = cols[inds]
          crows = col.rows[ind,]

          df$active[crows] = FALSE
          df$cluster.size[crows] = 2
          df$level[crows] = level

          df$cluster[crows] = min(df$cluster)-1
          df$var[crows] = syms[ccols]

          if (verbose) {
            cat("\neat.double.back:")
            cat("\nvars: ", paste0(syms[ccols], collapse=", "))
            cat("\nrows: ", paste0(crows, collapse=", "),"\n")
          }
          mat[crows,] = 0
          changed = changed + 2
        }
      }
      changed
    })
  )
}

