###################################################################################################
#
# Cross-validation for LASSO
#
# Tathagata Basu
# 11 Oct 2018
#
###################################################################################################

#' Cross-validation Tools

#' cross-validation
#' 
#' @param lambda penalty parameter
#' @param x predictors
#' @param y response
#' @param beta0 initial guess of beta
#' @param ts sequence of stepsize (for lasso_sg, lasso_pg)
#' @param method lasso optimization function
#' @param k number of fold
#' @param n_it number of iteration for lasso_cd method
#' @param relax relaxation on the accepted error should be between 0 and 1.
#' @export

cv.lasso<-function(lambda, x, y, beta0, ts, method, k=5, n_it=100, relax)
{
  
  data.partition<-cv.random.partition(x, y, k=k)
  
  # function for training
  
  cv.train<-function (i)
  {
    traindata<-as.matrix(apply(data.partition[,,-i], 2, rbind))
    x<-traindata[,2:ncol(traindata)]
    y<-traindata[,1]
    model<-method(lambda, x, y, beta0, ts, n_it=n_it)
  }
  model<-sapply(1:k,cv.train, simplify = "array")
  model<-array(unlist(model), dim = c(ncol(x), length(lambda), k))
  
  # function for testing
  
  cv.test<-function(i)
  {
    testdata<-as.matrix(data.partition[,,i])
    x<-testdata[,2:ncol(testdata)]
    y<-testdata[,1]
    estimate<-matrix(data = NA, nrow = nrow(testdata), ncol = length(lambda))
    model<-model[,,i]
    cv.predict<-function(beta, x) x%*%beta
    estimate<-as.matrix(apply(model, 2, cv.predict, x = x))
    error<-lapply(1:length(lambda),function(j){cv.mse(y, estimate[,j])})
    cv.error<-t(as.matrix(error))
    return(cv.error)
  }
  cv.error<-sapply(1:k,cv.test, simplify = "array")
  cv.error<-matrix(unlist(cv.error), nrow = k, ncol = length(lambda),  byrow = T)
  beta<-apply(model, 2, rowMeans)
  cv.model<-rbind(t(lambda), beta)
  rownames(cv.model)[2:nrow(cv.model)]<-sprintf("var %d",1:(nrow(cv.model)-1))
  m<-cv.model
  e<-cv.error
  me=colMeans(e)
  cv.error.index<-max(which(me<=(min(me)+(max(me)-min(me))*relax)))
  cv.coef<-m[,cv.error.index]
  output<-list("cv.model"=cv.model, "cv.error"=cv.error, "cv.coeff"=cv.coef, "index"=cv.error.index)
  return(output)
}

#' random partitioning for cross-validation
#' 
#' @param x predictors
#' @param y response
#' @param k number of fold
#' @export

cv.random.partition<-function(x, y, k)
{
  
  cv.data<-cbind(y,x)
  cv.data<-cv.data[1:((nrow(cv.data)%/%k)*k),]
  index.sample<-sample(1:nrow(cv.data), nrow(cv.data))
  cv.data<-cv.data[index.sample,]
  div<-(nrow(cv.data)/k)
  split.data<-sapply(1:k,function(i){cv.data[((i-1)*div+1):(i*div),]},simplify="array")
  
  return(split.data)
}

#' root mean-squared error
#' 
#' @param original original value
#' @param estimate estimated value
#' @export

cv.mse<-function(original, estimate)
{
  difference<-original-estimate
  difference.square<-difference^2
  total.error<-sum(difference.square)
  error<-total.error/length(estimate)
  
  return(error)
}

###################################################################################################

#' function for normalizing the data
#' 
#' @param x predictor
#' @param y response
#' @export

normalize<-function(x, y)
{
  # predictor
  x<-as.matrix(x)
  design<-as.matrix(rep(1, nrow(x)))
  mean.x<-colMeans(x)
  x<-x-design%*%mean.x
  sd.x<-apply(x,2,sd)
  x<-x/(design%*%sd.x)
  colnames(x)[1:ncol(x)]<-sprintf("var %d",1:ncol(x))
  
  # response
  y<-as.matrix(y)
  design<-as.matrix(rep(1, nrow(y)))
  mean.y<-colMeans(y)
  y<-y-design%*%mean.y
  sd.y<-sd(y)
  y<-y/(design%*%sd.y)
  colnames(y)<-"y"
  
  output<-list("x"=x,"y"=y)
  return(output)
}

###################################################################################################

#' lasso plots
#' @export

lasso.plot<-function(output, error, index, main=NULL)
{
  y.lim.up<-max(abs(output[2:nrow(output),]))
  y.lim.bel<-min(output[2:nrow(output),])
  cut<-colMeans(error)
  x<-output[1,]
  matplot(log(output[1,]),t(output[2:nrow(output),]),type = "l", 
          xlab = expression(paste(log(lambda))),
          ylab = "Coefficients", ylim = c(y.lim.bel,y.lim.up),
          main = main,
          lty = 1, col = 1:6)
  abline(h=0, col="black", lty=2)
  abline(v=log(x[max(index)]), lty=2)
}

#' cv curve plot
#' @export
cv.plot<-function(error,lambda, index, main=NULL)
{
  error<-as.matrix(error)
  y<-colMeans(error)
  ysd<-apply(error, 2, sd)
  x<-log(lambda)
  xHigh <- x
  yHigh <- y+ysd/2
  xLow <- x
  yLow <- y-ysd/2
  plot(x,y,xlab=expression(log(lambda)),ylab="mean-squared error",
       main=main, 
       ylim = c(min(yLow),max(yHigh)), col="red", pch=20)
  arrows(xHigh,yHigh,xLow,yLow,col="grey",angle=90,length=0.03,code=3)
  lines(x,y,col="red", lty=2)
  abline(v=x[max(index)], lty=2)
  abline(h=y[which(y==min(y))], lty=2)
}

###################################################################################################

# Examples
#' @export
cv.ex.1 = function(relax)
{
  cat(sprintf("\nSimulated Dataset: 24 predictors and 50 observations, no of fold is 5.\n\n"))
  x = matrix(data = rnorm(1200), nrow = 50, ncol = 24)
  b = as.matrix(rep(c(-3,-2,-1,1,2,3), 4))
  er = as.matrix(rnorm(50))
  y = x%*%b + er
  lambdas = exp(seq(-5,2,0.1))
  beta0 = b
  ts = opt_ts(0.1, 1000, 1000)
  third.test.cv<-cv.lasso(lambdas, x, y, beta0, ts, method = lasso_cd, relax=relax)
  m<-third.test.cv$cv.model
  e<-third.test.cv$cv.error
  cv.plot(e,m[1,], third.test.cv$index, main = "LASSO using co-ordinate descent")
  lasso.plot(m, e, third.test.cv$index)
  
  cat(sprintf("Least Square Model \n"))
  print(lm(y~x-1))
  
  cat(sprintf("LASSO estimates \n"))
  print(third.test.cv$cv.coeff)
}
