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

#' Sub-gradient optimization.
#' 
#' Function for sub-gradient optimization for non-differentiable functions
#' @param x Starting value.
#' @param f Function to optimize.
#' @param df Any subgradient of f.
#' @param ts Sequence of step sizes. use opt_ts() for stepsize generation
#' @return The function returns the optimal solution
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
#' Function for proximal-gradient optimization for non-differentiable functions
#' We're trying to minimize f + g; where f is the differentiable part and g is the non-differentiable part.
#' @param x Starting value.
#' @param f Function to optimize.
#' @param df Gradient of f.
#' @param pg Proximal operator for g.
#' @param ts Sequence of step sizes. use opt_ts() for stepsize generation
#' @return The function returns the optimal solution
#' @export

pg_optim = function(x, f, df, pg, ts) {
  for (t in ts) x = pg(t, x - t * df(x))
  x
}

#' Co-ordinate descent optimization.
#' 
#' Function for co-ordinate descent optimization for non-differentiable functions
#' @param x Starting value.
#' @param f Function to optimize.
#' @param x_it Iterative x obtained from the subgradient function
#' @param n_it Number of iterations. Default value is 100
#' @return The function returns the optimal solution
#' @export

cd_optim = function(x, f, x_it, n_it = 100) {
  for (j in 1:n_it) {
    x.last = x
    for(i in 1:length(x)) {
      x[i] = x_it(x, i)
    }
    if(sum(abs(x.last - x)) < 0.00001)
      break
  }
  x[abs(x) < 0.00001] = 0
  x
}


#' Sequence of step sizes for optimization.
#' 
#' Function to generate stepsize for sub-gradient optimization and proximal-gradient optimization
#' @param t Starting value.
#' @param m Number of constant steps.
#' @param n Number of diminishing steps.
#' @return The sequence of stepsize
#' @export

opt_ts = function(t, m, n)
  if (n != 0) {
    c(rep(t, m), t / (1:n))
  } else {
    rep(t, m)
  }
