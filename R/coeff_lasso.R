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
#' @export
square_f = function(x, y, beta)
  sum((y - x %*% beta)^2) / (2 * nrow(x))

#' OLS differential
#' @export
square_df = function(x, y, beta)
  -t(x) %*% (y - x %*% beta) / nrow(x)

#' LASSO penalty term
#' @export
lasso_f = function(lambda, beta, wt) lambda * sum(abs(beta * wt))


#' LASSO penalty subgradient (Note: not actual gradient.)
#' @export
lasso_df = function(lambda, beta, wt) lambda * sign(beta) * wt

#' Proximal operator for lasso_f.
#' @export
lasso_p = function(lambda, wt) function(t, x) sign(x) * pmax(0, abs(x) - lambda * t * wt)

#' Soft threshold operator
#' @export
soft  = function(lambda, wt)function(x, i) sign(x) * max(0, abs(x) - lambda * wt[i])

#' CD LASSO Soft Threshold term
#' @export
st_f = function(i, x, y, beta)
  t(x[,i]) %*% (y - x[,-i] %*% beta[-i]) / (t(x[,i]) %*% x[,i])

#' LASSO objective
#' @export
square_lasso_f = function(lambda, x, y, beta, wt)
  square_f(x, y, beta) + lasso_f(lambda, beta, wt)

#' LASSO sub-gradient
#' @export
square_lasso_df = function(lambda, x, y, beta, wt)
  square_df(x, y, beta) + lasso_df(lambda, beta, wt)

#' LASSO optimization (sub-gradient method)
#' @param lambda Penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @export
lasso_optim_sg = function(lambda, x, y, beta0, ts, wt) {
  f = function(beta) square_lasso_f(lambda, x, y, beta, wt)
  df = function(beta) square_lasso_df(lambda, x, y, beta, wt)
  sg_optim(x = beta0, f = f, df = df, ts = ts)
}

#' LASSO optimization (proximal gradient)
#' @param lambda Penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @export
lasso_optim_pg = function(lambda, x, y, beta0, ts, wt) {
  f = function(beta) square_f(x, y, beta)
  df = function(beta) square_df(x, y, beta)
  pg = lasso_p(lambda, wt)
  pg_optim(x = beta0, f = f, df = df, pg = pg, ts = ts)
}

#' LASSO optimization using coordinate descent
#' @param lambda Penalty term
#' @param x Predictors
#' @param y Response
#' @param beta0 Initial guess of the regression coefficients
#' @param n_it Number of iterations. Default value 100
#' @param wt weights for the coefficients of weighted LASSO.
#' @export
lasso_optim_cd = function(lambda, x, y, beta0, n_it = 100, wt){
  s = soft(lambda, wt)
  f = function(beta) square_lasso_f(lambda, x, y, beta, wt)
  v = function(i, beta) st_f(i, x, y, beta)
  cd_optim(x = beta0, f = f, v = v, s = s, n_it = n_it)
}

#' Coefficient path (sg)
#' @param x Predictors
#' @param y Response
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @export
lasso_sg = function(x, y, ts, wt, ...) {
  
  beta0 = as.matrix(rep(1, ncol(x)))
  betalist = list(solve(t(x) %*% x) %*% (t(x) %*% y))
  lambdas = c(0)
  j = 1
  h = max(abs(solve(t(x) %*% x) %*% (t(x) %*% y)) + 5) / 50
  while (sum(abs(betalist[[j]])) > 0.0001) {
    j = j + 1
    lambdas[j] = exp(-5 + h * j)
    betalist[[j]] = lasso_optim_sg(lambdas[j], x, y, beta0, ts, wt)
    
  }
  
  betas = do.call(cbind, betalist)
  output = list("lambdas" = lambdas, "betas" = betas)
  return(output)
}

#' Coefficient path (pg)
#' @param x Predictors
#' @param y Response
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @export
lasso_pg = function(x, y, ts, wt, ...) {
  
  beta0 = as.matrix(rep(1, ncol(x)))
  betalist = list(solve(t(x) %*% x) %*% (t(x) %*% y))
  lambdas = c(0)
  j = 1
  h = max(abs(solve(t(x) %*% x) %*% (t(x) %*% y)) + 5) / 50
  while (sum(abs(betalist[[j]])) > 0.0001) {
    j = j + 1
    lambdas[j] = exp(-5 + h * j)
    betalist[[j]] = lasso_optim_pg(lambdas[j], x, y, beta0, ts, wt)
    
  }
  
  betas = do.call(cbind, betalist)
  output = list("lambdas" = lambdas, "betas" = betas)
  return(output)
}

#' Coefficient path (cd)
#' @param x Predictors
#' @param y Response
#' @param n_it Number of iterations. Default value 100.
#' @param wt weights for the coefficients of weighted LASSO.
#' @export

lasso_cd = function(x, y, n_it=100, wt, ...){
  
  beta0 = as.matrix(rep(1, ncol(x)))
  betalist = list(solve(t(x) %*% x) %*% (t(x) %*% y))
  lambdas = c(0)
  j = 1
  h = max(abs(solve(t(x) %*% x) %*% (t(x) %*% y)) + 5) / 50
  while (sum(abs(betalist[[j]])) > 0.0001) {
    j = j + 1
    lambdas[j] = exp(-5 + h * j)
    betalist[[j]] = lasso_optim_cd(lambdas[j], x, y, beta0, n_it = 100, wt)
    
  }
  
  betas = do.call(cbind, betalist)
  output = list("lambdas" = lambdas, "betas" = betas)
  return(output)
}

#' sg Plot
#' @param x Predictors
#' @param y Response
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @export
lasso_sg_plot = function(x, y, ts, wt) {
  out = lasso_sg(x, y, ts, wt)
  matplot(
    log(out$lambdas), t(out$betas),
    type = "l", lty = 1,
    xlab = expression(paste(log(lambda))), ylab = expression(paste(beta)),
    main = "Sub-gradient Method")
}

#' pg plot
#' @param x Predictors
#' @param y Response
#' @param ts Stepsize for optimization method
#' @param wt weights for the coefficients of weighted LASSO.
#' @export
lasso_pg_plot = function(x, y, ts, wt) {
  out = lasso_pg(x, y, ts, wt)
  matplot(
    log(out$lambdas), t(out$betas),#t(do.call(cbind, betas)),
    type = "l", lty = 1,
    xlab = expression(paste(log(lambda))), ylab = expression(paste(beta)),
    main = "Proximal Gradient Method")
}

#' cd plot
#' @param x Predictors
#' @param y Response
#' @param n_it Number of iterations. Default value 100
#' @param wt weights for the coefficients of weighted LASSO.
#' @export
lasso_cd_plot = function(x, y, n_it = 100, wt) {
  out = lasso_cd(x, y, n_it, wt)
  matplot(
    log(out$lambdas), t(out$betas),
    type = "l", lty = 1,
    xlab = expression(paste(log(lambda))), ylab = expression(paste(beta)),
    main = "Co-ordinate Descent Method")
}

#' Examples 1
#' @import glmnet
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
  ts = opt_ts(0.1, 1000, 1000)
  lasso_sg_plot(x, y, ts, wt)
  lasso_pg_plot(x, y, ts, wt)
  lasso_cd_plot(x, y, 100, wt)
}
#' Examples 2
#' @import glmnet
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
  ts = opt_ts(0.1, 1000, 1000)
  lasso_sg_plot(x, y, ts, wt)
  lasso_pg_plot(x, y, ts, wt)
  lasso_cd_plot(x, y, 100, wt)
}
