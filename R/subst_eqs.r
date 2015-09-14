
symbolic.cluster.subst.eqs = function(df, cluster, rename.clusters=TRUE) {
  restore.point("cluster.subst.var")

  rows = which(df$cluster==cluster)
  res = symbolic.subst.eqs(eqs=df$eq_[rows], vars=df$var[rows])

  erows = rows[res$eq.inds]
  srows = rows[res$subst.inds]
  ns = length(srows)

  # No equation could be substituted
  if (ns==0)
    return(list(df=df,remaining.cluster=cluster, new.clusters=NULL))


  if (length(res$eqs)>0) {
    df$var[erows] = res$vars
    df$cluster.size[erows] = length(erows)
    df$eq_[erows] = res$eqs
  }

  df$var[srows] = res$subst.vars
  df$eq_[srows] = res$subst.eq
  df$expl_[srows] = lapply(res$subst.eq, function(eq) eq[[3]])
  df$solved[srows] = TRUE
  df$cluster.size[srows] = 1
  df$cluster[srows] = cluster + (1:ns) / (ns+1)
  df$level[srows] = df$level[rows[1]] + (1:ns) / (ns+1)

  if (rename.clusters) {
    df$cluster = rank(df$cluster,ties.method = "min")
    df$level = rank(df$level,ties.method = "min")
  }

  remaining.cluster = df$cluster[erows[1]]
  new.clusters = df$cluster[srows]

  df = df[order(df$cluster),]

  return(list(df=df,remaining.cluster=remaining.cluster, new.clusters=new.clusters))

}

symbolic.cluster.subst.var = function(df, cluster, var, eq.ind=1) {
  restore.point("cluster.subst.var")

  rows = which(df$cluster==cluster)
  res = symbolic.subst.var(df$eq_[rows], var=var, eq.ind=eq.ind)
  if (!res$ok) return(list(ok=FALSE,df=df,remaining.cluster=cluster, new.cluster=NA))

  erow = rows[eq.ind]
  df$eq_[[erow]] = res$subst.eq
  if (length(rows)==1)
    return(list(ok=TRUE,df=df,remaining.cluster=NA, new.cluster=cluster))

  syms = df$var[rows]
  rsyms = setdiff(syms,var)
  rrows = setdiff(rows,erow)
  df$cluster[erow] = df$cluster[erow]+0.5
  df$level[erow] = df$level[erow]+0.5

  df$cluster = rank(df$cluster,ties.method = "min")
  df$level = rank(df$level,ties.method = "min")
  df$cluster.size[rrows] = df$cluster.size[rrows[1]]-1
  df$cluster.size[erow] = 1
  df$var[rrows] = rsyms
  df$var[erow] = var

  remaining.cluster = df$cluster[rrows[1]]
  new.cluster = df$cluster[erow]

  df = df[order(df$cluster),]

  return(list(ok=TRUE,df=df,remaining.cluster=remaining.cluster, new.cluster=new.cluster))

}

symbolic.subst.var = function(eqs, var, eq.ind=1, eq.solve.fun=sym.solve.eq,  simplify.fun = Deriv::Simplify, expl=NULL) {
  restore.point("symbolic.subst.var")

  eq = eqs[[eq.ind]]
  if (is.null(expl)) {
    res = eq.solve.fun(eq,var)
    if (!res$solved) {
      return(list(ok=FALSE, eqs=eqs, subst.eq=NULL))
    }
    subst.eq = res$eq
    expl_ = subst.eq[[3]]
  } else {
    subst.eq = substitute(lhs == rhs, list(lhs=as.name(var), rhs=expl))
    expl_ = expl
  }

  # Try to simplify the explicit solution
  if (!is.null(simplify.fun)) {
    expl_ = simplify.fun(expl_)
    subst.eq = substitute(lhs==rhs, list(lhs=as.name(var),rhs=expl_))
  }
  reqs = eqs[-eq.ind]
  # substitute solution into remaining equations
  subst.li = list(substitute((expl_),list(expl_=expl_)))
  names(subst.li) = var
  reqs = lapply(reqs, substitute.call,  env=subst.li)
  if (!is.null(simplify.fun)) {
    reqs = lapply(reqs, simplify.fun)
  }
  return(list(ok=TRUE, eqs=reqs, subst.eq=subst.eq))
}

symbolic.subst.eqs = function(eqs,vars, subst.eqs=NULL, subst.vars=NULL, eq.inds=seq_along(eqs), subst.inds=NULL) {

  restore.point("symbolic.subst.eqs")
  if (length(eqs)==0) {
    return(list(eqs=eqs, vars=vars, subst.eqs=subst.eqs, subst.vars=subst.vars, eq.inds=eq.inds, subst.inds=subst.inds))
  }
  sug = suggest.subst.eq.and.var(eqs, vars)
  if (!sug$ok) {
    return(list(eqs=eqs, vars=vars, subst.eqs=subst.eqs, subst.vars=subst.vars, eq.inds=eq.inds, subst.inds=subst.inds))
  }
  var = sug$var; eq.ind = sug$eq.ind
  res = symbolic.subst.var(eqs = eqs,var = var,eq.ind = eq.ind,expl=sug$expl)

  # Call again with remaining values
  symbolic.subst.eqs(eqs = res$eqs, vars = setdiff(vars,var), subst.eqs=c(res$subst.eq,subst.eqs), subst.vars = c(var, subst.vars), eq.inds=setdiff(eq.inds,eq.inds[eq.ind]), subst.inds=c(eq.inds[eq.ind], subst.inds))
}

# Heuristics that suggest the next variable and equation to be substituted
# when solving a system of equations by substitution method
suggest.subst.eq.and.var = function(eqs,vars,eq.solve.fun=sym.solve.eq) {
  restore.point("suggest.subst.eq.and.var")

  lhs = lapply(eqs, function(eq) eq[[2]])
  rhs = lapply(eqs, function(eq) eq[[3]])
  lhs.single = sapply(lhs, is.name)
  rhs.single = sapply(rhs, is.name)

  lhs.var = sapply(seq_along(lhs), function(ind) {
    if (!lhs.single[ind]) return("")
    as.character(lhs[[ind]])
  })
  rhs.var = sapply(seq_along(rhs), function(ind) {
    if (!rhs.single[ind]) return("")
    as.character(rhs[[ind]])
  })

  lhs.single = lhs.single & lhs.var %in% vars
  rhs.single = rhs.single & rhs.var %in% vars



  identity = which(lhs.single & rhs.single)

  # If we have identities, propose them
  if (length(identity)>0) {
    row = identity[1]
    # substitute the longer variable
    lvar = as.character(lhs[[row]])
    rvar = as.character(rhs[[row]])
    if (nchar(lvar)>nchar(rvar) | (! rvar %in% vars) ) {
      return(list(ok=TRUE,eq.ind=row, var=lvar, expl=rhs[[row]], type="id"))
    } else {
      return(list(ok=TRUE,eq.ind=row, var=rvar, expl=lhs[[row]], type="id"))
    }
  }

  # check for explicit solutions
  var = sapply(c(lhs[lhs.single], rhs[rhs.single]), as.character)
  expl = c(rhs[lhs.single], lhs[rhs.single])
  is.expl = sapply(seq_along(expl), function(i) ! (var[i] %in% intersect(vars,find.variables(expl[[i]]) )))
  if (any(is.expl)) {
    var = var[is.expl]
    expl = expl[is.expl]
    eq.ind = c(which(lhs.single),which(rhs.single))[is.expl]
    # Take shortest expl but prefer longer variables names
    size = sapply(expl, call.size) - (nchar(var) / 100)
    row = which.min(size)
    return(list(ok=TRUE,eq.ind=eq.ind[row], var=var[row], expl=expl[[row]], type="expl"))
  }

  # go through all equations by size
  size = sapply(eqs, call.size)
  eq.inds = order(size)
  eq.ind = eq.inds[5]
  for (eq.ind in eq.inds) {
    eq = eqs[[eq.ind]]
    vt = count.variables(eq)
    evars = intersect(names(vt), vars)
    vt = vt[evars]
    ord = order(vt - nchar(evars) / 100)
    for (var in evars[ord]) {
      res = eq.solve.fun(eq,var)
      if (res$solved) {
        return(list(ok=TRUE,eq.ind=eq.ind, var=var, expl=res$eq[[3]], type="solved"))
      }
    }
  }

  # Could not solve any variable
  return(list(ok=FALSE,eq.ind=NA, var="", expl=NULL, type="unsolved"))
}
