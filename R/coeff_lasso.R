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

#' OLS term
#'
#' For internal use
#' @export

square_f = function(x, y, beta)
  sum((y - x %*% beta)^2) / (2 * nrow(x))

#' OLS differential
#'
#' For internal use
#' @export

square_df = function(x, y, beta)
  -t(x) %*% (y - x %*% beta) / nrow(x)

#' LASSO penalty term
#'
#' For internal use
#' @export

lasso_f = function(lambda, beta, wt) lambda * sum(abs(beta * wt))


#' LASSO penalty subgradient (Note: not actual gradient.)
#'
#' For internal use
#' @export

lasso_df = function(lambda, beta, wt) lambda * sign(beta) * wt

#' Proximal operator for lasso_f.
#'
#' For internal use
#' @export

lasso_p = function(lambda, wt) function(t, x) sign(x) * pmax(0, abs(x) - lambda * t * wt)

#' Soft threshold operator
#'
#' For internal use
#' @export

soft  = function(lambda, wt)function(x, i) sign(x) * max(0, abs(x) - lambda * wt[i])

#' CD LASSO Soft Threshold term
#'
#' For internal use
#' @export

st_f = function(i, x, y, beta)
  t(x[,i]) %*% (y - x[,-i] %*% beta[-i]) / (t(x[,i]) %*% x[,i])

#' LASSO objective
#'
#' For internal use
#' @export

square_lasso_f = function(lambda, x, y, beta, wt)
  square_f(x, y, beta) + lasso_f(lambda, beta, wt)

#' LASSO sub-gradient
#'
#' For internal use
#' @export

square_lasso_df = function(lambda, x, y, beta, wt)
  square_df(x, y, beta) + lasso_df(lambda, beta, wt)

#' LASSO optimization (sub-gradient method)
#'
#' This function is to be used for solving LASSO optimization using sub-gradient method
#' @param lambda Penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @return The function returns the coefficient of the predictors.
#' @export

lasso_optim_sg = function(lambda, x, y, beta0, ts, wt) {
  if ((is.null(wt) == T) | (length(wt) != ncol(x))) 
     wt = rep(1, ncol(x))
  else wt = ncol(x) * wt/sum(wt)
  
  f = function(beta) square_lasso_f(lambda, x, y, beta, wt)
  df = function(beta) square_lasso_df(lambda, x, y, beta, wt)
  sg_optim(x = beta0, f = f, df = df, ts = ts)
}

#' LASSO optimization (proximal gradient)
#'
#' This function is to be used for solving LASSO optimization using proximal-gradient method
#' @param lambda Penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @return The function returns the coefficient of the predictors.
#' @export

lasso_optim_pg = function(lambda, x, y, beta0, ts, wt) {
  if ((is.null(wt) == T) | (length(wt) != ncol(x))) 
     wt = rep(1, ncol(x))
  else wt = ncol(x) * wt/sum(wt)
  
  f = function(beta) square_f(x, y, beta)
  df = function(beta) square_df(x, y, beta)
  pg = lasso_p(lambda, wt)
  pg_optim(x = beta0, f = f, df = df, pg = pg, ts = ts)
}

#' LASSO optimization using co-ordinate descent
#'
#' This function is to be used for solving LASSO optimization using co-ordinate descent method
#' @param lambda Penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param n_it Number of iterations. Default value 100
#' @param wt weights for the coefficients of weighted LASSO.
#' @return The function returns the coefficient of the predictors.
#' @export

lasso_optim_cd = function(lambda, x, y, beta0, n_it = 100, wt){
  if ((is.null(wt) == T) | (length(wt) != ncol(x))) 
     wt = rep(1, ncol(x))
  else wt = ncol(x) * wt/sum(wt)
  
  s = soft(lambda, wt)
  f = function(beta) square_lasso_f(lambda, x, y, beta, wt)
  v = function(i, beta) st_f(i, x, y, beta)
  x_it = function(beta, i)s(v(i, beta), i)
  cd_optim(x = beta0, f = f, x_it = x_it, n_it = n_it)
}

#' Coefficient path (sg)
#'
#' This function is to be used for getting co-efficient path using sub-gradient optimization 
#' @param lambdas Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @return The function returns the coefficient of the predictors for each lambda.
#' @export

lasso_sg = function(lambdas, x, y, beta0, ts, wt, ...) {
  betas = lapply(lambdas, function(lambda) lasso_optim_sg(lambda, x, y, beta0, ts, wt))
  return(betas)
}

#' Coefficient path (pg)
#'
#' This function is to be used for getting co-efficient path using proximal-gradient optimization 
#' @param lambdas Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @return The function returns the coefficient of the predictors for each lambda.
#' @export

lasso_pg = function(lambdas, x, y, beta0, ts, wt, ...) {
  betas = lapply(lambdas, function(lambda) lasso_optim_pg(lambda, x, y, beta0, ts, wt))
  return(betas)
}

#' Coefficient path (cd)
#'
#' This function is to be used for getting co-efficient path using co-ordinate descent optimization 
#' @param lambdas Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param n_it Number of iterations. Default value 100.
#' @param wt weights for the coefficients of weighted LASSO.
#' @return The function returns the coefficient of the predictors for each lambda.
#' @export

lasso_cd = function(lambdas, x, y, n_it=100, wt, ...){
  beta0 = as.matrix(rep(0, ncol(x)))
  betas = lapply(lambdas, function(lambda) lasso_optim_cd(lambda, x, y, beta0, n_it, wt))
  return(betas)
}

#' sg Plot
#'
#' Plots the co-efficient path using sub-gradient optimization 
#' @param lambdas Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @export

lasso_sg_plot = function(lambdas, x, y, beta0, ts, wt) {
  betas = lasso_sg(lambdas, x, y, beta0, ts, wt)
  matplot(
    log(lambdas), t(do.call(cbind, betas)),
    type = "l", lty = 1,
    xlab = expression(paste(log(lambda))), ylab = expression(paste(beta)),
    main = "Sub-gradient Method")
}

#' pg plot
#'
#' Plots the co-efficient path using proximal-gradient optimization 
#' @param lambdas Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @export

lasso_pg_plot = function(lambdas, x, y, beta0, ts, wt) {
  betas = lasso_pg(lambdas, x, y, beta0, ts, wt)
  matplot(
    log(lambdas), t(do.call(cbind, betas)),
    type = "l", lty = 1,
    xlab = expression(paste(log(lambda))), ylab = expression(paste(beta)),
    main = "Proximal Gradient Method")
}

#' cd plot
#'
#' Plots the co-efficient path using co-ordinate descent optimization 
#' @param lambdas Values of the penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param n_it Number of iterations. Default value 100
#' @param wt weights for the coefficients of weighted LASSO.
#' @export

lasso_cd_plot = function(lambdas, x, y, n_it = 100, wt) {
  betas = lasso_cd(lambdas, x, y, n_it, wt)
  matplot(
    log(lambdas), t(do.call(cbind, betas)),
    type = "l", lty = 1,
    xlab = expression(paste(log(lambda))), ylab = expression(paste(beta)),
    main = "Co-ordinate Descent Method")
}

#' Example
#'
#' Example 1 of LASSO using different optimization techniques
#' @param wt weights for the coefficients of weighted LASSO. Defaults to NULL
#' @export

example_lasso_1 = function(wt = NULL) 
{
  if ((is.null(wt) == T)|(length(wt) != 3))
    wt = rep(1, 3)
  else
    wt = 3 * wt/sum(wt)

  x = matrix(data = rnorm(9000), ncol = 3)
  b = as.matrix(c(-3,0,3))
  er = as.matrix(rnorm(3000))
  y = x %*% b + er
  lambdas = exp(seq(-5,2,0.2))
  beta0 = as.matrix(rep(0,3))
  ts = opt_ts(0.1, 1000, 1000)
  lasso_sg_plot(lambdas, x, y, beta0, ts, wt)
  lasso_pg_plot(lambdas, x, y, beta0, ts, wt)
  lasso_cd_plot(lambdas, x, y, 100, wt)
}

#' Example
#'
#' Example 2 of LASSO using different optimization techniques
#' @param wt weights for the coefficients of weighted LASSO. Defaults to NULL
#' @export

example_lasso_2 = function(wt = NULL)
{
  if ((is.null(wt) == T)|(length(wt) != 24))
    wt = rep(1, 24)
  else
    wt = 24 * wt/sum(wt)
  
  x = matrix(data = rnorm(1200), nrow = 50, ncol = 24)
  b = as.matrix(rep(c(-3,-2,-1,1,2,3), 4))
  er = as.matrix(rnorm(50))
  y = x %*% b + er
  lambdas = exp(seq(-5,2,0.2))
  beta0 = b
  ts = opt_ts(0.1, 1000, 1000)
  lasso_sg_plot(lambdas, x, y, beta0, ts, wt)
  lasso_pg_plot(lambdas, x, y, beta0, ts, wt)
  lasso_cd_plot(lambdas, x, y, 100, wt)
}
