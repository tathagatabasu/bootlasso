###################################################################################################
#
# Lower and upper estimates of LASSO coefficients
#
# Tathagata Basu
# 15 Dec 2018
#
###################################################################################################

#' Lower and upper estimates of LASSO coefficients
#'
#' Lower and upper estimates of LASSO coefficients using optimization over weights
#' @param lambda value of the penalty parameter
#' @param x predictors
#' @param y response
#' @param wtl lower bound of the weights for the coefficients of weighted LASSO
#' @param wtu upper bound of the weights for the coefficients of weighted LASSO
#' @param ts stepsize for proximal gradient and sub-gradient method (use opt_ts() to generate stepsize). Defaults to NULL
#' @param method lasso optimization function.  Three different methods are available to use. method = c(lasso_optim_cd, lasso_optim_sg, lasso_optim_pg). Defaults to lasso_optim_cd
#' @param n_it number of iteration for lasso_optim_cd method. Default value is 10.
#' @return The function returns the lower and upper estimates of LASSO coefficients.
#' @export

imp.lasso = function(lambda, x, y, wtl, wtu, ts = NULL, method = lasso_optim_cd, n_it = 10)
{
	theta = (wtl + wtu)/2

	f = function(wt)
	{
		if ((is.null(wt) == T) | (length(wt) != ncol(x))) 
			wt = rep(1, ncol(x))
		else wt = ncol(x) * wt/sum(wt)
    
		beta = method(lambda = lambda,
		x, y, n_it = n_it, wt)
		return (beta[i])
	}

	ui = rbind(-diag(ncol(x)), diag(ncol(x)))

	ci = c(-wtu, wtl)


	betalow = matrix(nrow = ncol(x), ncol = 1)
	betaup = matrix(nrow = ncol(x), ncol = 1)

	for (i in 1:ncol(x))
	{
		low = constrOptim(theta, f, grad = NULL, ui, ci)
		betalow[i] = low$value
		up = constrOptim(theta, f, grad = NULL, ui, ci, control = list(fnscale = -1))
		betaup[i] = up$value
	}

	betabound = cbind(betalow, betaup)
	
	rownames(betabound)[1:nrow(betabound)] = sprintf("var %d", 1:nrow(betabound))
	colnames(betabound) = c("low", "up")

	return(betabound)
}

###################################################################################################

#' Tuning of LASSO coefficients bounds
#'
#' Tuning of the model over lambda using lower and upper estimates of LASSO coefficients using optimization over weights
#' @param lambdas values of the penalty parameter
#' @param x predictors
#' @param y response
#' @param wtl lower bound of the weights for the coefficients of weighted LASSO
#' @param wtu upper bound of the weights for the coefficients of weighted LASSO
#' @param ts stepsize for proximal gradient and sub-gradient method (use opt_ts() to generate stepsize). Defaults to NULL
#' @param method lasso optimization function.  Three different methods are available to use. method = c(lasso_optim_cd, lasso_optim_sg, lasso_optim_pg). Defaults to lasso_optim_cd
#' @param n_it number of iteration for lasso_optim_cd method. Default value is 10.
#' @return The function returns the lower and upper estimates of LASSO coefficients.
#' @export

imp.lasso.tuning = function(x, y, wtl, wtu, ts = NULL, method = lasso_optim_cd, n_it = 10)
{
	lmax = 1/nrow(x)*max(abs(t(x)%*%y))
	lambdas = exp(seq(-5, log(lmax), length.out = 51))
	betabounds = lapply(1:length(lambdas), function(i) imp.lasso(lambdas[i], x, y, wtl, wtu, ts = ts, method = method, n_it = n_it))
	error = matrix(ncol = 1, nrow = length(lambdas))
	
	f = function(beta)cv.mse(y, x%*%beta)
	
	for (i in 1:length(lambdas))
	{
		theta = (betabounds[[i]][,1] + betabounds[[i]][,2] + 10^-8) / 2
		ui = rbind(-diag(ncol(x)), diag(ncol(x)))
		ci = c(-(betabounds[[i]][,2]+10^-8), betabounds[[i]][,1])
		
		tune = constrOptim(theta, f, grad = NULL, ui, ci, control = list(fnscale = -1))
		error[i] = tune$value
	}
	
	output = list("betabounds" = betabounds, "maxerror" = error)
	return(output)
}

###################################################################################################

#' Example
#'
#' Example to check lower and upper estimates for LASSO.
#' @export

example.imp = function()
{
	x = matrix(rnorm(1200), nrow = 100, ncol = 12)
	b = rep(c(1, 2, 3), 4)
	er = rnorm(100)
	y = x%*%b + er

	wtl = rep(1/2, ncol(x))
	wtu = rep(1, ncol(x))

	ex = imp.lasso.tuning(x, y, wtl, wtu)
	
	print(ex$maxerror)
}
