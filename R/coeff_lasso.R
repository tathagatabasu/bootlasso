###########################################################################
#
# Coefficient Paths
#
# Matthias C. M. Troffaes
# 10 Oct 2018
#
# Additions : Tathagata Basu
# 11 Oct 2018
###########################################################################
#source("R/opt_lasso.R")
require(glmnet)

#' OLS term
#' @export
square_f = function(x, y, beta)
  sum((y - x %*% beta)^2) / (2 * nrow(x))

#' OLS differential
#' @export
square_df = function(x, y, beta)
  -t(x) %*% (y - x %*% beta) / nrow(x)

#' LASSO penalty term
#' @export
lasso_f = function(lambda, beta) lambda * sum(abs(beta))


#' LASSO penalty subgradient (Note: not actual gradient.)
#' @export
lasso_df = function(lambda, beta) lambda * sign(beta)

#' Proximal operator for lasso_f.
#' @export
lasso_p = function(lambda) function(t, x) sign(x) * pmax(0, abs(x) - lambda * t)

#' Soft threshold operator
#' @export
soft  = function(lambda)function(x) sign(x) * max(0, abs(x) - lambda)

#' LASSO objective
#' @export
square_lasso_f = function(lambda, x, y, beta)
  square_f(x, y, beta) + lasso_f(lambda, beta)

#' LASSO sub-gradient
#' @export
square_lasso_df = function(lambda, x, y, beta)
  square_df(x, y, beta) + lasso_df(lambda, beta)

#' LASSO optimization (sub-gradient method)
#' @param lambda Penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @export
lasso_optim_sg = function(lambda, x, y, beta0, ts) {
  f = function(beta) square_lasso_f(lambda, x, y, beta)
  df = function(beta) square_lasso_df(lambda, x, y, beta)
  sg_optim(x=beta0, f=f, df=df, ts=ts)
}

#' LASSO optimization (proximal gradient)
#' @param lambda Penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @export
lasso_optim_pg = function(lambda, x, y, beta0, ts) {
  f = function(beta) square_f(x, y, beta)
  df = function(beta) square_df(x, y, beta)
  g = function(beta) lasso_f(lambda, beta)
  p = lasso_p(lambda)
  pg_optim(x=beta0, f=f, df=df, g=g, p=p, ts=ts)
}

#' LASSO optimization using coordinate descent
#' @param lambda Penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @export
lasso_optim_cd = function(lambda, x, y, beta0, ...){
  s = soft(lambda)
  betas = beta0
  for (i in 1:ncol(x)) {
    betas[i]=s(t(x[,i])%*%(y-x[,-i]%*%betas[-i])/(t(x[,i])%*%x[,i]))
  }
  return(betas)
}

#' Coefficient path (sg)
#' @param lambda Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @export
lasso_sg = function(lambdas, x, y, beta0, ts) {
  betas = lapply(lambdas, function(lambda) lasso_optim_sg(lambda, x, y, beta0, ts))
  return(betas)
}

#' Coefficient path (pg)
#' @param lambda Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @export
lasso_pg = function(lambdas, x, y, beta0, ts) {
  betas = lapply(lambdas, function(lambda) lasso_optim_pg(lambda, x, y, beta0, ts))
  return(betas)
}

#' Coefficient path (cd)
#' @param lambda Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @export
lasso_cd = function(lambdas, x, y, beta0, ...){
  betas = lapply(lambdas, function(lambda)lasso_optim_cd(lambda, x, y, beta0))
  return(betas)
}

#' sg Plot
#' @param lambda Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @export
lasso_sg_plot = function(lambdas, x, y, beta0, ts) {
  betas = lasso_sg(lambdas, x, y, beta0, ts)
  matplot(
    log(lambdas), t(do.call(cbind, betas)),
    type="l", lty=1,
    xlab=expression(paste(log(lambda))), ylab=expression(paste(beta)),
    main="Sub-gradient Method")
}

#' pg plot
#' @param lambda Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @export
lasso_pg_plot = function(lambdas, x, y, beta0, ts) {
  betas = lasso_pg(lambdas, x, y, beta0, ts)
  matplot(
    log(lambdas), t(do.call(cbind, betas)),
    type="l", lty=1,
    xlab=expression(paste(log(lambda))), ylab=expression(paste(beta)),
    main="Proximal Gradient Method")
}

#' cd plot
#' @param lambda Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @export
lasso_cd_plot = function(lambdas, x, y, beta0, ...) {
  betas = lasso_cd(lambdas, x, y, beta0, ts)
  matplot(
    log(lambdas), t(do.call(cbind, betas)),
    type="l", lty=1,
    xlab=expression(paste(log(lambda))), ylab=expression(paste(beta)),
    main="Co-ordinate Descent Method")
}

#' glmnet plot
#' @param lambda Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @export
glmnet_optim_plot = function(lambdas, x, y) {
  glmfit = glmnet(x, y, alpha=1, family="gaussian")
  plot(
    glmfit, xvar="lambda", xlim=c(min(log(lambdas)), max(log(lambdas))),
    xlab=expression(paste(log(lambda))), ylab=expression(paste(beta)))
  abline(h=0, col="black", lty=2)
}

#' Examples 1
#' @export
lasso_test_1 = function() {
  x = matrix(data = rnorm(9000), ncol = 3)
  b = as.matrix(c(-3,0,3))
  er = as.matrix(rnorm(3000))
  y = x %*% b + er
  lambdas = exp(seq(-5,2,0.2))
  beta0 = as.matrix(c(-3,0,3))
  ts = opt_ts(0.1, 1000, 1000)
  lasso_sg_plot(lambdas, x, y, beta0, ts)
  lasso_pg_plot(lambdas, x, y, beta0, ts)
  lasso_cd_plot(lambdas, x, y, beta0, ts)
  glmnet_optim_plot(lambdas, x, y)
}
#' Examples 2
#' @export
lasso_test_2 = function()
{
  x = matrix(data = rnorm(1200), nrow = 50, ncol = 24)
  b = as.matrix(rep(c(-3,-2,-1,1,2,3), 4))
  er = as.matrix(rnorm(50))
  y = x%*%b + er
  lambdas = exp(seq(-5,2,0.2))
  beta0 = b
  ts = opt_ts(0.1, 1000, 1000)
  lasso_sg_plot(lambdas, x, y, beta0, ts)
  lasso_pg_plot(lambdas, x, y, beta0, ts)
  lasso_cd_plot(lambdas, x, y, beta0, ts)
  glmnet_optim_plot(lambdas, x, y)
}
