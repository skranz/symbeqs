# Tools for symbolically simplifying or solving systems of equations

The main function is `cluster.equations` that tries to separate the equation systems into subsets of equations, called cluster, that can be solved sequentially. A singleton cluster consists of a single equation and variable and may also be solved symbollically. Clusters that cannot be reduced to singletons must always be solved numerically in the end.

### Basic example:



```r
  library(symbeqs)
  eqs = list(
    quote(x+y == z),
    quote(x ==5),
    quote(z + a +y == 3)
  )
  df = cluster.equations(eqs, endo=c("x","y","z"))
```

```
## 
## eat.single.front:
## vars:  x
## rows:  2 
## 
## eat.double.front: 
## vars:  y, z
## rows:  2, 3
```

```r
  as.data.frame(df)
```

```
##   var            expl_                   eq_           org_ solved level
## 1   x                5                x == 5         x == 5   TRUE     1
## 2   y -((a + x - 3)/2) y == -((a + x - 3)/2) z + a + y == 3   TRUE     2
## 3   z            x + y            z == x + y     x + y == z   TRUE     3
##   cluster cluster.size    vars vars.id num.vars org.ind val is.check.eq
## 1       1            1       x       x        1       2  NA       FALSE
## 2       2            1 y, a, x   a|x|y        2       3  NA       FALSE
## 3       3            1 z, x, y   x|y|z        3       1  NA       FALSE
##   is.free.var
## 1       FALSE
## 2       FALSE
## 3       FALSE
```


Here we get exlicit solutions for all variables (column "expl_"). Note that the explizit solution for "y" contains the value of "x". That is because "x" can be computed before "y". The solution for "y" also contains the symbol "a", which is considered an exogenous parameter.

We can solve the equations numerically (exploiting all explicit solutions and simplifications into clusters)


```r
eval.cluster.df(df, exo=list(a=5))
```

```
##    x    y    z 
##  5.0 -3.5  1.5
```

Alternatively, we can build a function that allows to quickly compute this solution for given exogenous parameters

... NOT YEt NICELY IMPLEMENTED ...

###  An equation system with one more equation than variables:


```r
  # More equations than variables
  eqs = list(
    quote(x+y == 3),
    quote(6 == 2*(x+y)),
    quote(x ==5)
  )
  as.data.frame(df)
```

```
##   var            expl_                   eq_           org_ solved level
## 1   x                5                x == 5         x == 5   TRUE     1
## 2   y -((a + x - 3)/2) y == -((a + x - 3)/2) z + a + y == 3   TRUE     2
## 3   z            x + y            z == x + y     x + y == z   TRUE     3
##   cluster cluster.size    vars vars.id num.vars org.ind val is.check.eq
## 1       1            1       x       x        1       2  NA       FALSE
## 2       2            1 y, a, x   a|x|y        2       3  NA       FALSE
## 3       3            1 z, x, y   x|y|z        3       1  NA       FALSE
##   is.free.var
## 1       FALSE
## 2       FALSE
## 3       FALSE
```

The last equation has the flag `is.check.eq == TRUE`. This means this equation must hold true given the values of "x" and "y" computed by the earlier equations.

## An example with perfect colliniarity

We now have an example with 3 equations and 3 variables. Unfortunately, equations 1 and 2 are collinear, so that we have effectively only 2 equations. When solving it, we get one check.equation and the assigned variable is a free variable.

```r
  # Collinearity
  eqs = list(
    quote(x+y+z == 3),
    quote(6 == 2*(x+y+z)),
    quote(x ==5)
  )
  df = cluster.equations(eqs, endo=c("x","y","z"))
```

```
## 
## eat.single.front:
## vars:  x
## rows:  3 
## 
## eat.double.front: 
## vars:  y, z
## rows:  2, 3
```

```r
  as.data.frame(df)
```

```
##   var        expl_               eq_                 org_ solved level
## 1   x            5            x == 5               x == 5   TRUE     1
## 2   z            0            0 == 0 6 == 2 * (x + y + z)   TRUE     1
## 3   y -(x + z - 3) y == -(x + z - 3)       x + y + z == 3   TRUE     3
##   cluster cluster.size    vars vars.id num.vars org.ind val is.check.eq
## 1       1            1       x       x        1       3  NA       FALSE
## 2       2            1    NULL                3       2  NA        TRUE
## 3       3            1 y, x, z   x|y|z        3       1  NA       FALSE
##   is.free.var
## 1       FALSE
## 2        TRUE
## 3       FALSE
```

```r
  # Free variables
  df$var[df$is.free.var]
```

```
## [1] "z"
```

```r
  # Value of z is undetermined and by default set to 0
  eval.cluster.df(df)
```

```
##  x  z  y 
##  5  0 -2
```
By default, we set the value of a free variable to 0. But you may pick any other value by overwriting df$expl_.

## An example with more variables than equations


```r
  # Too many variables
  eqs = list(
    quote(x+y+z == 3),
    quote(x ==5)
  )
  df = cluster.equations(eqs, endo=c("x","y","z"))
```

```
## 
## We have 1 more variables than equations. Include dummy equations 0==0 to match number of equations with number of variables.
## eat.single.front:
## vars:  x
## rows:  2 
## In eat.single back:
##   - The variables y, z are only defined in the single equation 2:
##     x + y + z == 3
## Skip eat.single.back.
```

```r
  as.data.frame(df)
```

```
##   var        expl_               eq_           org_ solved level cluster
## 1   x            5            x == 5         x == 5   TRUE     1       1
## 2   z            0            0 == 0         0 == 0   TRUE     2       2
## 3   y -(x + z - 3) y == -(x + z - 3) x + y + z == 3   TRUE     3       3
##   cluster.size    vars vars.id num.vars org.ind val is.check.eq
## 1            1       x       x        1       2  NA       FALSE
## 2            1    NULL                0       3  NA        TRUE
## 3            1 y, x, z   x|y|z        3       1  NA       FALSE
##   is.free.var
## 1       FALSE
## 2        TRUE
## 3       FALSE
```

```r
  # Free variables
  df$var[df$is.free.var]
```

```
## [1] "z"
```

```r
  # Value of z is undetermined and by default set to 0
  eval.cluster.df(df)
```

```
##  x  z  y 
##  5  0 -2
```

## Internals

We always want to have as many variables as equations.

If orginally, we have more variables than equations. We add dummy equations 0 == 0 and continue simplifying. 

If we have more equations than variables, we add DUMMY variables. An equation that remains marked with a DUMMY variable will be a check equation.


