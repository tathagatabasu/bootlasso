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
	
	cv = lapply(1:nsim, function(i) cv.lasso(lambdas = lambdas, x = x, y = y, wt = wts[i,], ts = ts, method = method, k = k, n_it = n_it, df = df))
	coeff = lapply(1:nsim, function(i) cv[[i]]$coeff)
				   
	out = matrix(unlist(coeff), nrow = ncol(x)+2, ncol = nsim)
	coef = out[2:(ncol(x)+1),]
	
	rownames(coef)[1:nrow(coef)] = sprintf("var %d", 1:nrow(coef))
	
	#coefficient distribution plot
  
    boxplot(t(coef),
          main = "Boxplot of Regression Coefficients",
          ylab = "value of coeffs")
    points(rowMeans(coef), col = "red")
    abline(h = 0, lty = 2)
    legend("topright", legend = c("mean"), col = c("red"), pch = c(19), lty = 1, cex = 0.8)
	
	return(summary(coef))
	}