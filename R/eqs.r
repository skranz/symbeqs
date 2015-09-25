examples.with.econ.models = function() {

  set.restore.point.options(display.restore.point = TRUE)

  library(EconModels)
  # Model builder
  setwd("D:/libraries/EconModels/EconModels")
  initEconModels()
  em = load.model("ThreeEq")

  init.model(em)

  options(warn = 2)
  res = testwise.init.model(em)

}

examples.cluster.equations = function() {


  eqs = alist(y = y == y_eq, pi = pi == piT, pi_mr = pi_mr ==      piT, y_mr = y_mr == y_eq, Epi = Epi == lag_pi, y = y == A -      a * lag_r, pi = pi == Epi + alpha * (y - y_eq), next_Epi = next_Epi ==      pi, pi_mr = y_mr == y_eq - alpha * beta * (pi_mr - piT),      y_mr = pi_mr == next_Epi + alpha * (y_mr - y_eq), rcut = y_mr ==          A - a * rcut, r = r == pmax(rcut, -Epi), i = i == r +          Epi)

  endo = c("y", "pi", "pi_mr", "y_mr", "Epi", "lag_pi", "lag_r", "next_Epi",  "rcut", "r", "i")

  df = cluster.equations(eqs, endo=endo)

  eval.cluster.df(df, make.random.params(df=df))

  df



  eqs = list(
    quote(x+y == z),
    quote(x ==5),
    quote(z + a +y == 3)
  )
  df = cluster.equations(eqs, endo=c("x","y","z"))


  options(warn=2)
  # More equations than variables
  eqs = list(
    quote(x+y == 3),
    quote(x ==5),
    quote(6 == 2*(x+y))
  )

  suggest.eq.for.var(eqs=eqs, var="y")
  df = cluster.equations(eqs, endo=c("x","y"))

  find.check.eqs(df=df)

}


#' Try to solve or simplify system of equations
cluster.equations = function(eqs, var.li=NULL, endo=NULL, exo=NULL, verbose=TRUE, solve.symbolic=TRUE, skip.eat=FALSE, solve.level = 100, skip.big=FALSE, funs=NULL, cluster.df=NULL) {
  restore.point("cluster.equations")


  eqs = compute.equation.funs(eqs,funs)
  # Find for each formula the contained endogenous variables
  if (is.null(var.li)) {
    var.li = lapply(eqs, function(form) {
      vars = find.variables(form)
      vars = setdiff(vars, exo)
      if (!is.null(endo)) vars = intersect(vars, endo)
      vars
    })
  }
  syms = unique(unlist(var.li))

  # deal with free symbols (more variables than equations)
  # add dummy equations
  num.free = length(syms)-length(eqs)
  if (num.free>0) {
    if (verbose)
      cat("\nWe have", num.free, "more variables than equations. Include dummy equations 0==0 to match number of equations with number of variables.")
    free.rows = (NROW(df)+1):(NROW(df)+num.free)
    eqs = c(eqs, replicate(n=num.free, quote(0==0)))
    var.li = c(var.li, vector("list",num.free))
  } else if (num.free < 0) {
    if (verbose)
      cat("\nWe have", -num.free, "more equations than variables. Include dummy variables DUMMY_1_,...")
    syms = c(syms, paste0("DUMMY___", 1:(-num.free)))

  }
  nr = length(eqs)
  vars.id = sapply(var.li, function(vars) {
    if (length(vars)==0) return("")
    paste0(sort(vars), collapse="|")
  })

  df = data_frame(var=syms,expl_ = vector("list",nr),eq_ = eqs,org_=eqs, solved=FALSE, level=1,cluster=1,cluster.size=length(eqs), vars=var.li, vars.id = vars.id, num.vars=sapply(var.li, length),org.ind=seq_along(eqs), val=NA_real_)


  if (!skip.eat) {
    df = eat.from.cluster(df, cluster=1)
  }
  if (solve.symbolic)  {
    df = solve.symbolic.cluster.df(df, skip.big=skip.big)
  }

  df = cluster.df.update.var.eq.info(df)

  df
}

#' Update df$vars, df$vars.id, df$is.free.var, and df$is.check.eq
cluster.df.update.var.eq.info = function(df) {
  restore.point("cluster.df.update.var.eq.info")

  # Update variable list, given the new equations
  endo = df$var
  df$vars = lapply(df$eq_, function(form) {
    intersect(find.variables(form),endo)
  })

  df$vars.id = sapply(df$vars, function(vars) {
    if (length(vars)==0) return("")
    paste0(sort(vars), collapse="|")
  })

  # Find check.eqs and possible free variables
  df$is.check.eq = sapply(seq_along(eqs), function(i) {
    ! (df$var[i] %in% df$vars[[i]])
  })

  # We only have a check.eq if all equations in the cluster are check.eqs
  df = mutate(group_by(df, cluster), is.check.eq = all(is.check.eq))

  # non-singleton clusters of check eqs will be separated
  rows = which(df$cluster.size >= 2 & df$is.check.eq)
  if (length(rows)>0) {
    df$cluster[rows] = df$cluster[rows]+ seq_along(rows) / (length(rows)+1)
    df$cluster = rank(df$cluster,ties.method = "min")
    df$cluster.size[rows] = 1
  }

  # If a check.eq has not a dummy variable, the assigned variable is free
  df$is.free.var = df$is.check.eq & !str.starts.with(df$var,"DUMMY___")

  # Set free variables to zero and mark them as solved
  df$expl_[df$is.free.var] = as.list(rep(0, sum(df$is.free.var)))
  df$solved[df$is.free.var] = TRUE
  df
}

solve.symbolic.cluster.df = function(df, skip.big=FALSE, skip.small=FALSE, eq.solve.fun=sym.solve.eq,  simplify.fun = simplify.eq) {
  restore.point("solve.symbolic.cluster.df")

  if (!skip.small) {
    rows = which(df$cluster.size == 1 & !df$solved)
    row = rows[1]
    for (row in rows) {
      res = eq.solve.fun(df$eq_[[row]], df$var[row])
      if (!res$solved) next
      expl = simplify.fun(res$eq[[3]])
      df$eq_[[row]] = substitute(lhs==rhs,list(lhs=as.name(df$var[row]), rhs=expl))
      df$expl_[[row]] = expl
      df$solved[row] = TRUE
    }
  }

  if (!skip.big) {
    clusters = unique(df$cluster[df$cluster.size>1])

    for (cluster in clusters) {
      res = symbolic.cluster.subst.eqs(df = df, cluster=cluster,rename.clusters = FALSE)
      df = res$df
    }
    if (length(clusters)>0) {
      # now reorder clusters and levels
      df$cluster = rank(df$cluster,ties.method = "min")
      df$level = rank(df$level,ties.method = "min")
      df = df[order(df$cluster),,drop=FALSE]
    }
  }
  df
}

symbolic.cluster.equations = function(eqs, exo=NULL, endo=NULL, eq.solve.fun=sym.solve.eq, level=-1, subst=NULL, skip.big=FALSE, simplify.fun = simplify.eq) {
  restore.point("symbolic.cluster.equations")

  if (length(eqs)==2 & "K" %in% endo) restore.point("nfjdmfgdkgrngdgfgfhzbgvf")
  vars = unique(unlist(lapply(eqs, find.variables)))
  if (is.null(endo)) {
    endo = setdiff(vars, exo)
  } else {
    free.var = setdiff(endo,vars)
    if (length(free.var)>0) {
      warning("The variable(s) ", paste0(free.var)," cannot be determined. Please specify an initial value for them.")
    }
  }

  if (length(endo) != length(eqs)) {
    stop("number of endo variables and number of equations is not identical")
  }

  # check if one equation is already in form of an explicit solution
  solved = FALSE
  if (! (skip.big & length(eqs)>1)) {
    for (eq.ind in seq_along(eqs)) {
      var = is.expl.eq(eqs[[eq.ind]],endo)
      if (!is.null(var)) {
        new.eq = eqs[[eq.ind]]
        solved = TRUE
        break
      }
    }
    # try to solve an equation for an endogenous variable
    # using the symbolic equation solver
    if (!solved) {
      sym.li = lapply(eqs, function(eq) find.variables(eq))
      for (eq.ind in seq_along(eqs)) {
        eq = eqs[[eq.ind]]
        vars = intersect(sym.li[[eq.ind]],endo)
        for (var in vars) {
          res = eq.solve.fun(eq,var)
          if (res$solved) {
            solved = TRUE
            new.eq = res$eq
            break
          }
        }
        if (solved) break
      }
    }
  }
  # could solve an equation
  if (solved) {
    expl_ = new.eq[[3]]

    # Try to simplify the explicit solution
    if (!is.null(simplify.fun))
      expl_ = simplify.fun(expl_)

    # try to solve remaining equations,
    # plugging in the solved value
    if (length(eqs)>1) {
      #restore.point("symbolic.cluster.equations.outer")
      reqs = eqs[-eq.ind]
      child.inds = seq_along(eqs)[-eq.ind]


      # substitute solution into remaining equations
      subst.li = list(substitute((expl_),list(expl_=expl_)))
      names(subst.li) = var
      reqs = lapply(reqs, substitute.call,  env=subst.li)

      rdf = symbolic.cluster.equations(reqs, endo=setdiff(endo,var),eq.solve.fun = eq.solve.fun,level=level-1, subst=c(subst, var))
      rdf$org.ind = child.inds[rdf$org.ind]
    } else {
      rdf = NULL
    }
    restore.point("nhfnfngigtirmgjnng")

    eq.df = data_frame(var=var,expl_=list(expl_), eq_ = list(new.eq), cluster.size = 1, solved=TRUE, level=level, org.ind=eq.ind, subst=list(subst))
    df = rbind(rdf, eq.df)
  } else {
    if (!is.null(simplify.fun))
      eqs = lapply(eqs, simplify.fun)

    df = data_frame(var=endo,expl_=vector("list",length(endo)), eq_ = eqs, cluster.size = length(endo), solved=FALSE, level=level, org.ind = seq_along(eqs), subst=replicate(length(endo),subst,simplify = FALSE) )
  }
  df
}

# is an equation already in explicit form with the variable on the lhs?
is.expl.eq = function(eq, vars=NULL) {
  restore.point("is.expl.eq")
  lhs = eq[[2]]
  if (!is.name(lhs)) return(NULL)

  lhs.var = as.character(lhs)
  # a different lhs variable than the endogenous variables
  if (!is.null(vars)) {
    if (!lhs.var %in% vars) return(NULL)
  }
  # lhs variable is also on rhs
  if (lhs.var %in% find.variables(eq[[3]])) return(NULL)

  # indeed in explicit form
  lhs.var
}


compute.symbolic.ops = function(call, recursive = TRUE) {
  restore.point("compute.symbolic.ops")
  if (is.call(call)) {
    name = as.character(call[[1]])
    if (name=="Deriv") {
      restore.point("compute.symbolic.ops.Deriv")
      call = eval(call)
    }
  }
  if (length(call)>1 & recursive) {
    for (i in seq_along(call)[-1])
      call[[i]] = compute.symbolic.ops(call=call[[i]], recursive=TRUE)
  }
  call
}

replace.variable.funs = function(call, funs, recursive=TRUE) {
  restore.point("replace.variable.funs")
  if (is.call(call)) {
    name = as.character(call[[1]])
    # A call to be replaced
    if (name %in% names(funs)) {
      restore.point("replace.variable.funs.inner")
      args = lapply(seq_along(call)[-1], function(i) call[[i]])
      vars = names(call)
      if (is.null(vars) & length(args)>0) {
        stop(paste0("Error when calling the user defined function ", deparse1(call),":\n You must name all arguments!"))
      } else {
        if (any(vars[-1]=="") & length(args)>0) {
          stop(paste0("Error when calling the user defined function ", deparse1(call),":\n You must name all arguments!"))
        }
      }
      call = funs[[name]]
      if (length(vars)>0) {
        names(args) = vars[-1]
        call = substitute.call(call, args)
      }
    }
  }
  if (length(call)>1 & recursive) {
    for (i in seq_along(call)[-1])
      call[[i]] = replace.variable.funs(call=call[[i]],funs=funs, recursive=TRUE)
  }
  call
}


# replace functional references to other functions
# and calls to Deriv
compute.equation.funs = function(eqs, funs=NULL) {
  restore.point("compute.equation.funs")

  # replace functional references to other functions
  # do not yet check for cycles
  counter = 0
  changed = rep(TRUE,length(eqs))
  while (sum(changed)>0) {
    counter = counter +1
    for (eq.ind in which(changed)) {
      eq = replace.variable.funs(eqs[[eq.ind]], funs)
      changed = !identical(eq, eqs[[eq.ind]])
      if (changed) eqs[[eq.ind]] = eq
    }
    if (counter > 1000) {
      stop("To many function substitutions. You probably have cycles in your function references!")
    }
  }

  # replace Deriv functions
  eqs = lapply(eqs,compute.symbolic.ops)
  eqs
}


extract.explicit.eqs = function(eqs) {
  # extract explicitly defined variables
  # they can be referenced as functions
  is.explicit= sapply(seq_along(eqs), function(i) {
    eq = eqs[[i]]
    var = names(eqs[i])
    if (as.character(eq[[2]]) != var) return(FALSE)
    rhs.vars = find.variables(eq[[3]])
    if (var %in% rhs.vars) return(FALSE)
    return(TRUE)
  })
  funs = lapply(eqs[is.explicit], function(eq) eq[[3]])
  funs
}


extract.cluster.df.dummies = function(df) {
  rows = str.starts.with(df$var,"DUMMY___")
  return(list(df=df[!rows,,drop=FALSE], test.eqs = df$org_[rows]))
}
