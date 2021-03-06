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
#' @param x predictors
#' @param y response
#' @param wt weights for the coefficients of weighted LASSO. Defaults to NULL
#' @param ts stepsize for proximal gradient and sub-gradient method (use opt_ts() to generate stepsize). Defaults to NULL
#' @param method lasso optimization function.  Three different methods are available to use. method = c(lasso_cd, lasso_sg, lasso_pg). Defaults to lasso_cd
#' @param k number of fold. Default value is 5.
#' @param n_it number of iteration for lasso_cd method. Default value is 10.
#' @param df Degree of freedom. Number of desired variables to be zero. Defaults to NULL
#' @param nsim number of perturbed weights. Default value is 50
#' @param rel_err Relative purturbed error. Default value is 1
#' @return The function returns the summary of sensitivity analysis.
#' @export

sensitivity.lasso = function(x, y, wt = NULL, ts = NULL, method = lasso_cd, k = 5, n_it = 10, df = NULL, nsim = 50, rel_err = 1)
	{
	
	if ((is.null(wt) == T)|(length(wt) != ncol(x)))
		wt = rep(1, ncol(x))
	else
		wt = ncol(x) * wt / sum(wt)
	
	lmax = 1/nrow(x)*max(abs(t(x)%*%y))
	lambdas = exp(seq(-5, log(lmax), length.out = 101))
	wts = t(matrix(rep(wt, nsim), ncol = nsim) + wt * (rel_err / 100) * matrix(rnorm(nsim * ncol(x)), ncol = nsim))
	
	sensty= function(x, y, wt, ts, method, k, n_it, df)
	{
	  cv = cv.lasso(x = x, y = y, wt = wt, ts = ts, method = method, k = k, n_it = n_it, df = df)
	  out = cv$coeff
	  
	  beta = out[2:(ncol(x)+1)]
	  return(beta)
	}
	
	coef = sapply(1:nsim, function(i)unlist(sensty(x = x, y = y, wt = wts[i,], 
	                                               ts = ts, method = method, k = k, n_it = n_it, df = df)))
	
	rownames(coef)[1:nrow(coef)] = sprintf("var %d", 1:nrow(coef))
	
	#-----------------------------------------------------------------------------
    # Means, medians and SDs of the bootstrapped statistics
    #-----------------------------------------------------------------------------
  
    sensty.means = rowMeans(coef, na.rm=T)
  
    sensty.medians = apply(coef, 1, median, na.rm=T)
  
    sensty.sds = apply(coef, 1, sd, na.rm=T)
  
    #-----------------------------------------------------------------------------
    # Basic bootstrap CIs
    #-----------------------------------------------------------------------------
  
    conf.mat = matrix(apply(coef, 1 ,quantile, c(0.025, 0.975), na.rm=T),
                     ncol=2, byrow=TRUE)
    colnames(conf.mat) = c("95%-CI Lower", "95%-CI Upper")
  
    #-----------------------------------------------------------------------------
    # Summary
    #-----------------------------------------------------------------------------
  
    summary.frame = data.frame(mean=sensty.means, median=sensty.medians, sd=sensty.sds, 
                              "CI_lower" = conf.mat[,1], "CI_upper" = conf.mat[,2])
    sensty_summary = summary.frame
  
  
	#coefficient distribution plot
  
    boxplot(t(coef),
          main = "Boxplot of Regression Coefficients",
          ylab = "value of coeffs")
    points(rowMeans(coef), col = "red")
    abline(h = 0, lty = 2)
    legend("topright", legend = c("mean"), col = c("red"), pch = c(19), lty = 1, cex = 0.8)
	
	return(sensty_summary)
	}

###################################################################################################

#' Example
#'
#' Example to check cross-validation for LASSO.
#' @param wt weights for the coefficients of weighted LASSO. Defaults to NULL
#' @param nsim number simulations with perturbed weights. Default value is 200.
#' @param rel_err Relative purturbed error. Default value is 1
#' @export

example.sensitivity = function(wt=NULL, nsim = 200, rel_err = 1)
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
  
  sensitivity.lasso(x, y, wt = wt, ts = NULL, method = lasso_cd, k = 5, n_it = 10, df = NULL, nsim = nsim, rel_err = rel_err)
  
}
