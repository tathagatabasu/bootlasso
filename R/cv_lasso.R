###################################################################################################
#
# Cross-validation for LASSO
#
# Tathagata Basu
# 11 Oct 2018
#
###################################################################################################

#' Cross-validation
#' 
#' Cross-validation for lasso
#' @param x predictors
#' @param y response
#' @param wt weights for the coefficients of weighted LASSO. Defaults to NULL
#' @param ts stepsize for proximal gradient and sub-gradient method (use opt_ts() to generate stepsize). Defaults to NULL
#' @param method lasso optimization function.  Three different methods are available to use. method = c(lasso_cd, lasso_sg, lasso_pg). Defaults to lasso_cd
#' @param k number of fold. Default value is 5.
#' @param n_it number of iteration for lasso_cd method. Default value is 10.
#' @param df Degree of freedom. Number of desired variables to be zero. Defaults to NULL
#' @return The function returns the  following list of outputs
#' \item{model}{The values of coefficients corresponding to each lambda}
#' \item{error}{Mean-squared error of the cross-validated model corresponding to each lambda}
#' \item{coeff}{The values of coefficients corresponding to minimum mean-squared error}
#' \item{index}{Index of the minimum mean-squared error (for internal use)}
#' @export

cv.lasso = function(x, y, wt = NULL, ts = NULL, method = lasso_cd, k = 5, n_it = 10, df = NULL)
{
  data.partition = cv.random.partition(x, y, k = k)
  x = scale(x, scale = F)
  y = scale(y, scale = F)
  lmax = max(abs(t(x) %*% y / diag(t(x)%*%x)))
  lambdas = as.matrix(exp(seq(-5, log(lmax), length.out = 51)))
  colnames(lambdas) = "lambda"
  

  if ((is.null(wt) == T)|(length(wt) != ncol(x)))
    wt = rep(1, ncol(x))
  else
    wt = ncol(x) * wt / sum(wt)
  
  # function for training
  
  cv.train = function (i)
  {
    traindata = as.matrix(apply(data.partition[,,-i], 2, rbind))
    x = traindata[,2:ncol(traindata)]
    y = traindata[,1]
    model = method(lambdas = lambdas, x = x, y = y, ts = ts, n_it = n_it, wt = wt)
  }
  model = sapply(1:k, cv.train, simplify = "array")
  model = array(unlist(model), dim = c(ncol(x), length(lambdas), k))
  
  # function for testing
  
  cv.test = function(i)
  {
    testdata = as.matrix(data.partition[,,i])
    x = testdata[,2:ncol(testdata)]
    y = testdata[,1]
    estimate = matrix(data = NA, nrow = nrow(testdata), ncol = length(lambdas))
    model = model[,,i]
    cv.predict = function(beta, x) x %*% beta
    estimate = as.matrix(apply(model, 2, cv.predict, x = x))
    error = lapply(1:length(lambdas), function(j) cv.mse(y, estimate[,j]))
    cv.error = t(as.matrix(error))
    return(cv.error)
  }
  
  cv.error = sapply(1:k, cv.test, simplify = "array")
  cv.error = matrix(unlist(cv.error), nrow = k, ncol = length(lambdas),  byrow = T)
  
  beta = apply(model, 2, rowMeans)
  cv.model = rbind(t(lambdas), beta)
  rownames(cv.model)[2:nrow(cv.model)] = sprintf("var %d", 1:(nrow(cv.model) - 1))
  
  m = cv.model
  e = cv.error
  me = colMeans(e)
  cv.error.index = max(which(me <= (min(me))))
  
  vdf = as.matrix(colSums(beta != 0))
  er = as.matrix(me); colnames(er) = "error"
  cv.df = cbind(lambdas, vdf, er)
  colnames(cv.df)[2] = "df"

  mse = mean(e[,cv.error.index])
  cv.coef = c(m[,cv.error.index], mse)
  names(cv.coef)[length(cv.coef)] = "mse"
  
  if(is.null(df)==F)
  {
    cv.error.index = which(colMeans(e)==min(colMeans(e)[which(vdf <= df)]))
    mse = mean(e[,cv.error.index])
    cv.coef = c(m[,cv.error.index], mse)
    names(cv.coef)[length(cv.coef)] = "mse"
  }
  
  output = list("model" = cv.model, "error" = cv.error, "coeff" = cv.coef, "index" = cv.error.index,
                "df" = cv.df)
  return(output)
}

#' Random partition
#' 
#' For internal use
#' @param x predictors
#' @param y response
#' @param k number of fold
#' @return Partitioned dataset
#' @export

cv.random.partition = function(x, y, k)
{
  cv.data = cbind(y, x)
  cv.data = cv.data[1:((nrow(cv.data) %/% k) * k),]
  
  index.sample = sample(1:nrow(cv.data), nrow(cv.data))
  cv.data = cv.data[index.sample,]
  div = (nrow(cv.data) / k)
  
  split.data = sapply(1:k, function(i) cv.data[((i - 1) * div + 1):(i * div),], simplify="array")
  
  return(split.data)
}

#' Root mean-squared error
#' 
#' For internal use
#' @param original original value
#' @param estimate estimated value
#' @return The root-mean squared error
#' @export

cv.mse = function(original, estimate)
{
  difference = original - estimate
  difference.square = difference ^ 2
  total.error = sum(difference.square)
  error = total.error / length(estimate)
  
  return(error)
}


###################################################################################################

#' Cross-validation plot
#'
#' Plots the cross-valiadtion curve
#' @param cv model obtained by doing cross-validation.
#' @param main Title of the plot.
#' @export

cv.plot = function(cv, main=NULL)
{
  error = as.matrix(cv$error)
  y = colMeans(cv$error)
  ysd = apply(cv$error, 2, sd)
  yHigh = y + ysd / 2
  yLow = y - ysd / 2
  
  x = log(cv$model[1,])
  xHigh = x
  xLow = x
  
  plot(x, y, xlab = expression(paste(log(lambda))),
       ylab="mean-squared error",
       main=main, 
       ylim = c(min(yLow), max(yHigh)), col = "red", pch = 20)
  
  arrows(xHigh, yHigh, xLow, yLow, col = "grey", angle = 90, length = 0.03, code = 3)
  lines(x, y, col = "red", lty = 2)
  abline(v = x[max(cv$index)], lty = 2)
}

###################################################################################################

#' Example
#'
#' Example to check cross-validation for LASSO.
#' @param wt weights for the coefficients of weighted LASSO. Defaults to NULL
#' @export

example.cv = function(wt=NULL)
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
  
  
  test.cv = cv.lasso(x, y, wt = wt)
  
  cv.plot(test.cv, main = "Cross-validation error (LASSO using co-ordinate descent)")

  lasso_cd_plot(x, y, wt = wt)
    
  cat(sprintf("Least Square Model \n"))
  print(lm(y ~ x - 1))
  
  cat(sprintf("LASSO estimates \n"))
  print(test.cv$coeff)
  
  cat(sprintf("Degrees of freedom w.r.t. lambda \n"))
  print(test.cv$df)
}
