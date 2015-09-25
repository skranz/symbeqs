simplify.eq = function(expr, env=parent.frame(), scache=new.env()) {
  res = Deriv::Simplify(expr, env, scache)
  if (identical(res,TRUE)) {
    res = quote(0==0)
  }
  res
}
