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
#' @param lambdas values of the penalty parameter
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

cv.lasso = function(lambdas, x, y, wt = NULL, ts = NULL, method = lasso_cd, k = 5, n_it = 10, df = NULL)
{
  data.partition = cv.random.partition(x, y, k = k)
  lambdas = as.matrix(lambdas)
  colnames(lambdas) = "lambda"
  beta0 = rep(0, ncol(x))
  
  if ((is.null(wt) == T)|(length(wt) != length(beta0)))
    wt = rep(1, length(beta0))
  else
    wt = length(beta0) * wt / sum(wt)
  
  # function for training
  
  cv.train = function (i)
  {
    traindata = as.matrix(apply(data.partition[,,-i], 2, rbind))
    x = traindata[,2:ncol(traindata)]
    y = traindata[,1]
    model = method(lambdas = lambdas, x = x, y = y, beta0 = beta0, ts = ts, n_it = n_it, wt = wt)
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
  
  nvar = as.matrix(colSums(beta != 0))
  vdf = ncol(x) - nvar; colnames(vdf) = "df"
  er = as.matrix(me); colnames(er) = "error"
  cv.df = cbind(lambdas, vdf, er)
  
  s = min(which(nvar == 0))
  
  if (s == Inf)
    s = length(lambdas)
	
  if(cv.error.index > s)
    cv.error.index = s
  
  mse = mean(e[,cv.error.index])
  cv.coef = c(m[,cv.error.index], mse)
  names(cv.coef)[length(cv.coef)] = "mse"
  
  if(is.null(df)==F)
  {
    cv.error.index = which(colMeans(e)==min(colMeans(e)[which(nvar <= (ncol(x) - df))]))
    mse = mean(e[,cv.error.index])
    cv.coef = c(m[,cv.error.index], mse)
    names(cv.coef)[length(cv.coef)] = "mse"
  }
  
  output = list("model" = cv.model[,1:s], "error" = cv.error[,1:s], "coeff" = cv.coef, "index" = cv.error.index,
                "df" = cv.df[1:s,])
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

#' Normalizing the data
#' 
#' Function for normalizing the dataset
#' @param x predictor
#' @param y response
#' @return The function returns following list of outputs
#' \item{x}{normalized predictor}
#' \item{y}{normalized response}
#' @export

normalize = function(x, y)
{
  # predictor
  x = as.matrix(x)
  design = as.matrix(rep(1, nrow(x)))
  mean.x = colMeans(x)
  x = x - design %*% mean.x
  sd.x = apply(x, 2, sd)
  x = x / (design %*% sd.x)
  colnames(x)[1:ncol(x)] = sprintf("var %d", 1:ncol(x))
  
  # response
  y = as.matrix(y)
  design = as.matrix(rep(1, nrow(y)))
  mean.y = colMeans(y)
  y = y - design %*% mean.y
  sd.y = sd(y)
  y = y / (design %*% sd.y)
  colnames(y) = "y"
  
  output = list("x" = x, "y" = y)
  return(output)
}

###################################################################################################

#' LASSO plot
#'
#' Plots the LASSO coefficient path obtained from cross-validation.
#' @param cv model obtained by doing cross-validation.
#' @param main Title of the plot.
#' @export

lasso.plot = function(cv, main=NULL)
{
  y.lim.up = max(abs(cv$model[2:nrow(cv$model),]))
  y.lim.bel = min(cv$model[2:nrow(cv$model),])
  x = cv$model[1,]
  
  matplot(log(cv$model[1,]), t(cv$model[2:nrow(cv$model),]), type = "l", 
          xlab = expression(paste(log(lambdas))),
          ylab = "Coefficients", ylim = c(y.lim.bel, y.lim.up),
          main = main,
          lty = 1, col = 1:6)
  abline(h = 0, col = "black", lty = 2)
  abline(v = log(x[max(cv$index)]), lty = 2)
}

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
  
  plot(x, y, xlab = expression(paste(log(lambdas))),
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
  
  lambdas = exp(seq(-5,3,0.1))
  
  test.cv = cv.lasso(lambdas, x, y, wt = wt)
  
  cv.plot(test.cv, main = "Cross-validation error (LASSO using co-ordinate descent)")
  
  lasso.plot(test.cv, main = "LASSO using co-ordinate descent")
  
  cat(sprintf("Least Square Model \n"))
  print(lm(y ~ x - 1))
  
  cat(sprintf("LASSO estimates \n"))
  print(test.cv$coeff)
  
  cat(sprintf("Degrees of freedom w.r.t. lambda \n"))
  print(test.cv$df)
}
