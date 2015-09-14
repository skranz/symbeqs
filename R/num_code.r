

eval.cluster.df = function(clu.df, exo = list()) {
  restore.point("eval.cluster.df")

  exo = lapply(exo, as.numeric)
  vals = rep(NA_real_,NROW(clu.df))
  names(vals)  = clu.df$var
  clusters = sort(unique(clu.df$cluster))

  cluster = 1
  for (cluster in clusters) {
    fun = make.cluster.solve.fun(cluster=cluster, clu.df=clu.df)
    sol = try(fun(vals=c(exo,vals)))
    if (is(sol,"try-error")) break
    if (!is.list(sol)) {
      vals[[names(sol)]] = sol
    } else {
      if (!sol$ok) break
      cat("\nsolved cluster ", cluster)
      vals[names(sol$par)] = sol$par
    }
  }
  vals
}

numerical.solve.eqs = function(eqs,vars, exo) {
  solve.fun = numerical.solve.eqs.fun(eqs,vars)
  solve.fun(exo)
}

numerical.solve.eqs.fun = function(eqs, vars) {
  restore.point("eqs.numerical.solve.fun")

  syms = unique(unlist(lapply(eqs,find.variables)))

  xsubst = lapply(seq_along(vars), function(i) {
    substitute(x[i],list(i=i))
  })
  names(xsubst) = vars

  exo = setdiff(syms, vars)

  vsubst = lapply(exo, function(par) {
    substitute(vals[[par]],list(par=par))
  })
  names(vsubst)=exo

  # Create function
  impl_ = lapply(eqs, function(eq_) {
    call = substitute(lhs-(rhs), list(lhs=eq_[[2]],rhs=eq_[[3]]))
    call = substitute.call(call, c(xsubst,vsubst))
    call
  })
  names(impl_) = NULL
  inner = as.call(c(list(as.symbol("c")),impl_))

  fn = function(x,vals) {}
  body(fn) = substitute({inner}, list(inner=inner))
  fn

  #code = deparse1(body(fn))
  #code = sep.lines(code,",")
  #code

  var.assign.li = lapply(seq_along(vars), function(i) {
    substitute(var <- x[i], list(var=as.name(vars[i]),i=i))
  })
  var.assign = as.call(c(list(as.symbol("{")),var.assign.li))

  code = substitute({
    x <- runif(len.vars)
    fn <- fun
    #sol = mynleqslv(x=x, fn=fn)
    sol = optims(par=x, eq.fn=fn, vals=vals)
    if (!sol$ok) {
      warning(paste0("Could not solve ", paste0(endo, collapse=", "),":\n", sol$message ))
      sol$par[] = NA_real_
    }
    names(sol$par) = endo
    sol
  }, list(fun=fn,len.vars=length(vars),endo=vars, var.assign=var.assign))

  solve.fun = function(vals) {}
  body(solve.fun) = code
  solve.fun
}

make.cluster.solve.fun = function(cluster=clu.df$cluster[[1]], clu.df) {
  restore.point("make.cluster.solve.code")

  df = clu.df[clu.df$cluster==cluster,]

  vars = df$var
  solved = df$solved[1]
  if (solved & NROW(df)>1)
    stop("We have a solved cluster with more than 1 row")

  if (solved) {
    expl = df$expl_[[1]]
    syms = find.variables(expl)
    vsubst = lapply(syms, function(par) {
      substitute(vals[[par]],list(par=par))
    })
    names(vsubst)=syms
    expl = substitute.call(expl, vsubst)
    code = substitute({val <- rhs; names(val)=var; val}, list(var=df$var, rhs = expl))
    solve.fun = function(vals) {}
    body(solve.fun) <- code
  } else {
    solve.fun = numerical.solve.eqs.fun(eqs=df$eq_,vars)
  }

  solve.fun

}



eqs.solution.report = function(eqs, var.li=NULL, exo=NULL, verbose=TRUE, funs=NULL) {
  restore.point("eqs.solution.report")

  txt = NULL
  w = function(...) txt <<- c(txt,...)
  write.cluster.df = function(df) {
    levels = unique(df$level)
    for (lev in levels) {
      clusters = unique(df$cluster[df$level==lev])
      w(paste0("Level ", lev))
      for (cluster in clusters) {
        rows = which(df$cluster==cluster)
        w(paste0("  Cluster ", cluster, " (",length(rows),")"))
        str = paste0("    ",rows,". ",df$var[rows],"=",signif(df$val[rows],3), " : ", sapply(df$eq_[rows],deparse1))
        if (any(!is.finite(df$val[rows])) & min(rows)>1) {
          brows = 1:(min(rows)-1)
          vals = c(exo,df$val[brows] )
          subst = lapply(vals, signif, digits=4)
          names(subst) = c(names(exo),df$var[brows])
          num.eqs = lapply(df$eq_[rows], substitute.call, env=subst)
          str = paste0(str, "  <=>  ", num.eqs)
        }
        w(str)
      }
    }
  }



  w("# Part 1: Original equations\n",
    paste0("  ",seq_along(eqs),". ", sapply(eqs,deparse1))
  )

  # Solve functions and derivatives
  eqs = compute.equation.funs(eqs,funs)
  # Transform into clusters
  df = cluster.equations(eqs = eqs, var.li=var.li, exo=names(exo), solve.symbolic=FALSE, skip.eat=TRUE, skip.big=TRUE)

  df$val = eval.cluster.df(clu.df = df, exo=exo)

  w("\n# Part 2: Substitute functions and derivatives, put in one big cluster and try to solve\n",
    paste0("  ",seq_along(eqs),". ", sapply(eqs,deparse1))
  )
  write.cluster.df(df)

  w("
\n# Part 3: Eat away equations from big cluster

            from front: equations with single cluster variables
            from back:  cluster variables contained in a single equation\n")

  test.df = df;
  cluster = test.df$cluster[which.max(test.df$cluster.size)]
  test.df = eat.from.cluster(df=test.df,cluster=cluster, repeated=FALSE,eating.funs = list(eat.single.front))
  any(duplicated(test.df$var))

  cluster = test.df$cluster[which.max(test.df$cluster.size)]
  test.df = eat.from.cluster(df=test.df,cluster=cluster, repeated=FALSE,eating.funs = list(eat.single.back))
  any(duplicated(test.df$var))



  df = eat.from.cluster(df=df, cluster=1)
  df$val = eval.cluster.df(clu.df = df, exo=exo)

  write.cluster.df(df)

  #cluster.subst.var(df = df, cluster=cluster, var = "Wdem",eq.ind=7)
  w("
\n# Part 4: Try to solve single equations and systems (by substitution method) symbollically\n
    Note: The algorithm does not automatically detect multicolinarity, i.e. underdeterminiation
    of the equation system.
    If the results seem unreasonable, try adding conditions or fix values.\n\n")

  df = solve.symbolic.cluster.df(df)
  df$val = eval.cluster.df(clu.df = df, exo=exo)
  write.cluster.df(df)

  cat(paste0(txt,collapse="\n"))

  res = symbolic.cluster.subst.eqs(df=df, cluster=7)
  ndf = res$df

  rows = which(df$cluster==7)
  eqs = df$eq_[rows]
  vars = df$var[rows]
  symbolic.subst.eqs(eqs=eqs, vars=vars)


  df = df
  vals = df$val
  names(vals) = df$var
  exo.vals = c(exo, vals)
}


cluster.df.nice.code = function(df) {
  restore.point("cluster.df.nice.code")

  clusters = unique(df$cluster)
  clu = 1
  li = lapply(clusters, function(clu) {
    rows = which(df$cluster == clu)
    vars = df$var[rows]
    solved = df$solved[rows[1]]
    if (solved & length(rows)>1)
        stop("We have a solved cluster with more than 1 row")

    if (solved) {
      row = rows
      txt = paste0(
        "\n# Compute ", vars, " from ", deparse1(df$org_[[row]]),
        "\n",vars," = ", deparse1(df$expl_[[row]]),
        "\n",vars
      )
    } else {
      txt = paste0(
        "\n# Solve ", paste0(vars,collapse=", "), " from ",
        paste0("\n#  ", sapply(df$org_[rows],deparse1),collapse="")
      )


      xsubst = lapply(seq_along(vars), function(i) {
        substitute(x[i],list(i=i))
      })
      names(xsubst) = vars
      # Create function
      impl_ = lapply(df$eq_[rows], function(eq_) {
        call = substitute(lhs-(rhs), list(lhs=eq_[[2]],rhs=eq_[[3]]))
        call = substitute.call(call, xsubst)
        call

      })
      inner = as.call(c(list(as.symbol("c")),impl_))
      fn = function(x) {}
      body(fn) = substitute({inner}, list(inner=inner))
      fn

      code = substitute({
        x <- runif(length(endo))
        fn <- fun
        #sol = mynleqslv(x=x, fn=fn)
        sol = optims(par=x, eq.fn=fn)
        if (!sol$ok) {
          warning(paste0("Could not solve ", paste0(endo, collapse=", ")," in steady state.\n", sol$message ))
        }
      }, list(fun=fn,endo=vars))
      code.txt = paste0(sapply(code[-1], deparse1, collapse = "\n"), collapse="\n")
      txt = c(txt,"", code.txt)
      if (length(vars)>1) {
        txt = c(txt,
          "x = sol$par",
          paste0(vars, "= x[", seq_along(vars),"]", collapse="; "),
          paste0("c(",paste0(vars,collapse=","),")")
        )
      } else {
        txt = c(txt,
          paste0(vars, "= sol$par"),
          paste0(vars)
        )
      }
    }
    paste0(txt, collapse="\n")
  })




  code = paste0(li, collapse="\n\n")

  check.eq =  lapply(df$org_, function(eq) {
      deparse1(eq)
  })
  check.res = lapply(df$org_, function(eq) {
      paste0(deparse1(eq[[2]]),"-(",deparse1((eq[[3]])),")")
  })
  code = paste0(code,
"\n\ncheck.df <- data_frame(
  eq_ = qlist(", paste0(check.eq, collapse=",\n    "), "\n  ),

  lhs_rhs = c(",paste0(check.res, collapse=",\n    "), "\n  )
)
check.df$ok = abs(check.df$lhs_rhs) < 1e-8
check.df
")

  #cat(code)
  code
}

# Create fairly fast code that computes a cluster.df
adapt.make.inner.compute.code = function(df,  exo=NULL, subst.li=NULL) {
  restore.point("make.inner.compute.code")

  clusters = sort(unique(df$cluster))
  clu = 2
  li =lapply(clusters, function(clu) {
    rows = which(df$cluster == clu)
    if (all(df$solved[rows]) & (!all.implicit)) {
      code.li = lapply(rows, function(row) {
        code = substitute(lhs<-rhs,
             list(lhs=as.symbol(df$var[[row]]),rhs=df$expl_[[row]]))
        if (!is.null(subst.li))
          code = substitute.call(code, subst.li)
        code
      })
      return(code.li)
    } else {
      endo = df$var[rows]

      xsubst = lapply(seq_along(endo), function(i) {
        substitute(x[i],list(i=i))
      })
      names(xsubst) = endo

      # substitute original impl_ formula with x, par.mat and var.mat
      impl_ = lapply(df$impl_[rows], function(call) {
        call = substitute.call(call, xsubst)
        call = substitute.call(call, subst.li)
        call
      })
      inner = as.call(c(list(as.symbol("c")),impl_))

      assign = lapply(seq_along(endo), function(i) {
        substitute(var.mat[ti,var] <- res[[i]],list(i=i,var=endo[[i]]))
      })
      assign = as.call(c(list(as.symbol("{")),assign))

      if (lag.as.start) {
        #restore.point("ndjnfjdngjg")
        rhs = lapply(seq_along(endo), function(i) {
          substitute(var.mat[ti-1,var],list(var=endo[[i]]))
        })
        rhs = as.call(c(list(as.symbol("c")),rhs))

        start = substitute(start.x <- rhs, list(rhs=rhs))
      } else {
        start = substitute(start.x <- runif(num.endo), list(num.endo=length(endo)))
      }
      code = substitute({
          start
          sol = mynleqslv(x=start.x,ti=ti,par.mat=par.mat, var.mat=var.mat,
            fn = function(x,ti,par.mat, var.mat) {
            inner
          })
          if (max(abs(sol$fvec))> 1e-07) {
            warning(paste0("Could not solve ", paste0(endo, collapse=", "), " for ti = ", paste0(ti, collapse=", ")," max deviation = ",max(abs(sol$fvec)) ))
          }
          res = sol$x
          assign
        },
        list(inner=inner, start=start,assign=assign,endo=endo)
      )
      list(code)
    }

  })
  li = do.call(c, li)

  inner = as.call(c(list(as.symbol("{")),li))
  inner
}
