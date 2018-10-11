###########################################################################
#
# Optimization Methods
#
# Matthias C. M. Troffaes
# 10 Oct 2018
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
  for (t in ts) x = p(t, x - t * df(x))
  x
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
