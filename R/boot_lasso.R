#' Bootstrap for LASSO.
#'
#' @param lambdas values of the penalty parameter
#' @param x Predictors 
#' @param y Response
#' @param beta0 Initial guess of regression co-efficients
#' @param ts Stepsize for proximal gradient and sub-gradient method. Use opt_ts() to create your own.
#' @param method Optimization method. Three different methods are available to use. method = c(lasso_cd, lasso_sg, lasso_pg)
#' @param k No. of folds for cross-validation. Default value is 5.
#' @param n.sim No. of bootstrap replicates. Default is 500.
#' @export

boot.lasso<-function(lambdas, x, y, beta0=rep(0, ncol(x)), ts = opt_ts(0.1, 1000, 1000), method = lasso_cd, k=5, n.sim=500){

  #-----------------------------------------------------------------------------
  # Main funtion
  #-----------------------------------------------------------------------------
  
  cv1 <- function(x, y){
    data.cv.temp = cv.lasso(lambdas, x, y, beta0=beta0, ts=ts, method=method, k=k)
    data.fitlasso <- data.cv.temp$cv.model[-1,max(which(colMeans(data.cv.temp$cv.error)==min(colMeans(data.cv.temp$cv.error))))]
    lc=as.matrix(data.fitlasso)
    return(lc)}
  
  #-----------------------------------------------------------------------------
  # Original estimates
  #-----------------------------------------------------------------------------
  
  original.estimates <- t(cv1(x, y))
  
  #-----------------------------------------------------------------------------
  # Number of replications
  #-----------------------------------------------------------------------------
  
  n.sim <- n.sim
  
  #-----------------------------------------------------------------------------
  # The loop
  #-----------------------------------------------------------------------------
  
  set.seed(123)
  
  boot = function(i) {
    
    #-----------------------------------------------------------------------------
    # Sampling with replacement
    #-----------------------------------------------------------------------------
    
    mydata <- cbind(x, y)
    mydata.new <- mydata[sample(1:dim(mydata)[1], dim(mydata)[1], replace=TRUE),]
    X.new <- as.matrix(mydata.new[,1:ncol(x)])
    Y.new <- as.matrix(mydata.new[,ncol(mydata)])
    
    #-----------------------------------------------------------------------------
    # Bootstrap function
    #-----------------------------------------------------------------------------
    
    cv <- function(X.new, Y.new){
      data.cv.temp =cv.lasso(lambdas, x, y, beta0=beta0, ts=ts, method=method, k=k )
      data.fitlasso <- data.cv.temp$cv.model[-1,max(which(colMeans(data.cv.temp$cv.error)==min(colMeans(data.cv.temp$cv.error))))]
      lc=as.matrix(data.fitlasso)
      return(lc)}
    
    #-----------------------------------------------------------------------------
    # Results Matrix
    #-----------------------------------------------------------------------------
    
    store.matrix <- t(cv(X.new, Y.new))
    return(store.matrix)
  }
  
  store <- sapply(1:n.sim, boot, simplify = "array")
  
  store.matrix<-t(as.matrix(store[1,,]))
  colnames(store.matrix)[1:ncol(store.matrix)] <- sprintf("var %d", 1:ncol(store.matrix))
  
  #-----------------------------------------------------------------------------
  # Means, medians and SDs of the bootstrapped statistics
  #-----------------------------------------------------------------------------
  
  boot.means <- colMeans(store.matrix, na.rm=T)
  
  boot.medians <- apply(store.matrix,2,median, na.rm=T)
  
  boot.sds <- apply(store.matrix,2,sd, na.rm=T)
  
  #-----------------------------------------------------------------------------
  # Bootstrap bias
  #-----------------------------------------------------------------------------
  
  boot.bias <- colMeans(store.matrix, na.rm=T) - original.estimates
  
  #-----------------------------------------------------------------------------
  # Basic bootstrap CIs
  #-----------------------------------------------------------------------------
  
  conf.mat <- matrix(apply(store.matrix, 2 ,quantile, c(0.025, 0.975), na.rm=T),
                     ncol=2, byrow=TRUE)
  colnames(conf.mat) <- c("95%-CI Lower", "95%-CI Upper")
  
  #-----------------------------------------------------------------------------
  # Summary
  #-----------------------------------------------------------------------------
  
  summary.frame <- data.frame(mean=boot.means, median=boot.medians, bias=t(boot.bias),
                              sd=boot.sds, "CI_lower"=conf.mat[,1], "CI_upper"=conf.mat[,2])
  boot = summary.frame
  
  #-----------------------------------------------------------------------------
  # Graphs
  #-----------------------------------------------------------------------------
  
  #coefficient distribution
  
  
  boxplot(store.matrix,
          main="Boxplot of Regression Coefficients",
          ylab="value of coeffs")
  points(summary.frame$mean, col="red")
  abline(h=0, lty=2)
  legend("topright", legend = c("mean"), col=c("red"), pch = c(19), lty=1:2, cex=0.8)
  
  return(boot)
}
