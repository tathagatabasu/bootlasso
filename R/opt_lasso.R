###########################################################################
#
# Optimization Methods
#
# Matthias C. M. Troffaes
# 10 Oct 2018
#
# Additions : Tathagata Basu
# 11 Oct 2018
#
###########################################################################

#' Subgradient optimization.
#' 
#' @param x Starting value.
#' @param f Function to optimize.
#' @param df Any subgradient of f.
#' @param ts Sequence of step sizes.
#' @export

sg_optim = function(x, f, df, ts) {
  fx = f(x)
  x.best = x
  fx.best = fx
  for (t in ts) {
    x = x - t * df(x)
    fx = f(x)
    if(fx < fx.best) {
      x.best = x
      fx.best = fx
    }
  }
  x.best
}

#' Proximal gradient optimization.
#' 
#' We're trying to minimize f + g; prox is the proximal operator for g.
#' 
#' @param x Starting value.
#' @param f Function to optimize.
#' @param df Gradient of f.
#' @param p Proximal operator for g.
#' @param ts Sequence of step sizes.
#' @export

pg_optim = function(x, f, df, g, p, ts) {
  fx = f(x)
  x.best = x
  fx.best = fx
  for (t in ts) x = p(t, x - t * df(x))
  x
}

#' Co-ordinate descent optimization.
#' 
#' @param x Starting value.
#' @param f Function to optimize.
#' @param df Any subgradient of f.
#' @param n_it Number of iterations. Default value 100
#' @export

cd_optim = function(x, f, v, s, n_it) {
  fx = f(x)
  x.best = x
  fx.best = fx
  for (j in 1:n_it) {
    fx.last = fx.best
    for(i in 1:length(x)) {
      x[i] = s(v(i,x))
    }
    fx = f(x)
    if(fx < fx.best) {
      x.best = x
      fx.best = fx
    }
    if(abs(fx.last-fx)<0.00001)
      break
  }
  x.best
}


#' Generate a sequence of step sizes for optimization.
#' 
#' @param t Starting value.
#' @param m Number of constant steps.
#' @param n Number of diminishing steps.
#' @export

opt_ts = function(t, m, n)
  if (n != 0) {
    c(rep(t, m), t / (1:n))
  } else {
    rep(t, m)
  }
