###################################################################################################
#
# Bootstrap for LASSO
#
# Tathagata Basu
# 11 Oct 2018
#
###################################################################################################

#' Bootstrap for LASSO.
#'
#' Function for getting bootstrap estimates for LASSO.
#' @param lambdas values of the penalty parameter
#' @param x predictors 
#' @param y response
#' @param wt weights for the coefficients of weighted LASSO. Defaults to NULL
#' @param ts stepsize for proximal gradient and sub-gradient method (use opt_ts() to generate stepsize). Defaults to NULL
#' @param method optimization method. Three different methods are available to use. method = c(lasso_cd, lasso_sg, lasso_pg). Defaults to lasso_cd
#' @param k no. of folds for cross-validation. Default value is 5.
#' @param n_it no. of iteration for lasso_cd optimization. Default value is 10.
#' @param df degree of freedom. Number of desired variables to be zero. Defaults to NULL
#' @param n.sim no. of bootstrap replicates. Default value is 500.
#' @return The summary frame of the bootstrap comprising mean, median, bias, standard deviation and confidence intervals.
#' @export

boot.lasso = function(lambdas, x, y, wt = NULL, ts = NULL, 
                      method = lasso_cd, k = 5, n_it = 10, df = NULL, n.sim = 500)
  {

  #-----------------------------------------------------------------------------
  # Main funtion
  #-----------------------------------------------------------------------------
  
  cv1 = function(x, y)
    {
    data.cv.temp = cv.lasso(lambdas, x, y, wt = wt, ts = ts, 
                            method = method, n_it = n_it, k = k, df =df)
    data.fitlasso = data.cv.temp$coeff[2:(ncol(x)+1)]
    
    lc = as.matrix(data.fitlasso)
    return(lc)
    }
  
  #-----------------------------------------------------------------------------
  # Original estimates
  #-----------------------------------------------------------------------------
  
  original.estimates = t(cv1(x, y))
  
  #-----------------------------------------------------------------------------
  # Number of replications
  #-----------------------------------------------------------------------------
  
  n.sim = n.sim
  
  #-----------------------------------------------------------------------------
  # The loop
  #-----------------------------------------------------------------------------
  
  set.seed(123)
  
  boot = function(i) 
    {
    
    #-----------------------------------------------------------------------------
    # Sampling with replacement
    #-----------------------------------------------------------------------------
    
    mydata = cbind(x, y)
    mydata.new = mydata[sample(1:dim(mydata)[1], dim(mydata)[1], replace=TRUE),]
    
    X.new = as.matrix(mydata.new[,1:ncol(x)])
    Y.new = as.matrix(mydata.new[,ncol(mydata)])
    
    #-----------------------------------------------------------------------------
    # Bootstrap function
    #-----------------------------------------------------------------------------
    
    cv <- function(x, y)
      {
      data.cv.temp = cv.lasso(lambdas, x, y, wt = wt, ts = ts, 
                            method = method, n_it = n_it, k = k, df =df)
      data.fitlasso = data.cv.temp$coeff[2:(ncol(x)+1)]
      
      lc = as.matrix(data.fitlasso)
      return(lc)
      }
    
    #-----------------------------------------------------------------------------
    # Results Matrix
    #-----------------------------------------------------------------------------
    
    store.matrix = t(cv(X.new, Y.new))
    return(store.matrix)
  }
  
  store = sapply(1:n.sim, boot, simplify = "array")
  
  store.matrix = t(as.matrix(store[1,,]))
  colnames(store.matrix)[1:ncol(store.matrix)] = sprintf("var %d", 1:ncol(store.matrix))
  
  #-----------------------------------------------------------------------------
  # Means, medians and SDs of the bootstrapped statistics
  #-----------------------------------------------------------------------------
  
  boot.means = colMeans(store.matrix, na.rm=T)
  
  boot.medians = apply(store.matrix, 2, median, na.rm=T)
  
  boot.sds = apply(store.matrix, 2, sd, na.rm=T)
  
  #-----------------------------------------------------------------------------
  # Bootstrap bias
  #-----------------------------------------------------------------------------
  
  boot.bias = colMeans(store.matrix, na.rm=T) - original.estimates
  
  #-----------------------------------------------------------------------------
  # Basic bootstrap CIs
  #-----------------------------------------------------------------------------
  
  conf.mat = matrix(apply(store.matrix, 2 ,quantile, c(0.025, 0.975), na.rm=T),
                     ncol=2, byrow=TRUE)
  colnames(conf.mat) = c("95%-CI Lower", "95%-CI Upper")
  
  #-----------------------------------------------------------------------------
  # Summary
  #-----------------------------------------------------------------------------
  
  summary.frame = data.frame(mean=boot.means, median=boot.medians, bias=t(boot.bias),
                              sd=boot.sds, "CI_lower" = conf.mat[,1], "CI_upper" = conf.mat[,2])
  boot = summary.frame
  
  #-----------------------------------------------------------------------------
  # Graphs
  #-----------------------------------------------------------------------------
  
  #coefficient distribution
  
  x11()
  boxplot(store.matrix,
          main = "Boxplot of Regression Coefficients",
          ylab = "value of coeffs")
  points(summary.frame$mean, col = "red")
  abline(h = 0, lty = 2)
  legend("topright", legend = c("mean"), col = c("red"), pch = c(19), lty = 2, cex = 0.8)
  
  return(boot)
}

#' Example
#'
#' Example to check boot-strap LASSO
#' @param n.sim No. of bootstrap replicates. Default value is 200.
#' @export

example.boot = function(n.sim = 200)
{
  cat(sprintf("\nSimulated Dataset: 24 predictors and 50 observations, no of fold is 5.\n\n"))
  x = matrix(data = rnorm(1200), nrow = 50, ncol = 24)
  b = as.matrix(rep(c(-3,-2,-1,1,2,3), 4))
  er = as.matrix(rnorm(50))
  y = x %*% b + er
  
  lambdas = exp(seq(-5,2,0.1))
  
  boot.test = boot.lasso(lambdas, x, y, n.sim = n.sim)
  
  cat(sprintf("Least Square Model \n"))
  print(lm(y ~ x - 1))
  
  cat(sprintf("boot-strap LASSO estimates \n"))
  print(boot.test$median)
}
