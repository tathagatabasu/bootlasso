###################################################################################################
#
# Sensitivity using weights for LASSO
#
# Tathagata Basu
# 1 Dec 2018
#
###################################################################################################

#' Perturbation of weights
#'
#' Sensitivity analysis of predictors using perturbed weights
#' @param lambdas values of the penalty parameter
#' @param x predictors
#' @param y response
#' @param wt weights for the coefficients of weighted LASSO. Defaults to NULL
#' @param ts stepsize for proximal gradient and sub-gradient method (use opt_ts() to generate stepsize). Defaults to NULL
#' @param method lasso optimization function.  Three different methods are available to use. method = c(lasso_cd, lasso_sg, lasso_pg). Defaults to lasso_cd
#' @param k number of fold. Default value is 5.
#' @param n_it number of iteration for lasso_cd method. Default value is 10.
#' @param df Degree of freedom. Number of desired variables to be zero. Defaults to NULL
#' @param nsim number of perturbed weights. Default value is 50
#' @return The function returns the summary of sensitivity analysis.
#' @export

sensitivity_lasso = function(lambdas, x, y, wt = NULL, ts = NULL, method = lasso_cd, k = 5, n_it = 10, df = NULL, nsim = 50)
	{
	wts = t(matrix(rep(wt, nsim), ncol = nsim) + wt / 100 * matrix(rnorm(nsim * ncol(x)), ncol = nsim))
	
	cv = cv.lasso(lambdas = lambdas, x = x, y = y, wt = wt, ts = ts, method = method, k = k, n_it = n_it, df = df)
	out = cv$coeff
	
	beta = out[2:(ncol(x)+1)]
	lambda = out[1]
	
	coef = sapply(1:nsim, function(i)lasso_optim_cd(lambda = lambda, x, y, beta0 = beta, n_it = 100, wts[i,]))
	
	rownames(coef)[1:nrow(coef)] = sprintf("var %d", 1:nrow(coef))
	
	#coefficient distribution plot
  
    boxplot(t(coef),
          main = "Boxplot of Regression Coefficients",
          ylab = "value of coeffs")
    points(rowMeans(coef), col = "red")
    abline(h = 0, lty = 2)
    legend("topright", legend = c("mean"), col = c("red"), pch = c(19), lty = 1, cex = 0.8)
	
	return(summary(t(coef)))
	}

###################################################################################################

#' Example
#'
#' Example to check cross-validation for LASSO.
#' @param wt weights for the coefficients of weighted LASSO. Defaults to NULL
#' @param nsim number simulations with perturbed weights. Default value is 200.
#' @export

example.sensitivity = function(wt=NULL, nsim = 200)
{
  cat(sprintf("\nSimulated Dataset: 24 predictors and 50 observations, no of fold is 5.\n\n"))
  x = matrix(data = rnorm(1200), nrow = 50, ncol = 24)
  b = as.matrix(rep(c(-3,-2,-1,1,2,3), 4))
  er = as.matrix(rnorm(50))
  y = x %*% b + er
  if ((is.null(wt) == T)|(length(wt) != length(b)))
    wt = rep(1, length(b))
  else
    wt = length(b) * wt/sum(wt)
  
  lambdas = exp(seq(-5,3,0.1))
  
  sensitivity_lasso(lambdas, x, y, wt = wt, ts = NULL, method = lasso_cd, k = 5, n_it = 10, df = NULL, nsim = nsim)
  
}
