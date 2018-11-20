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
#' @param lambdas penalty parameter
#' @param x predictors
#' @param y response
#' @param beta0 initial guess of beta
#' @param wt weights for the coefficients of weighted LASSO.
#' @param ts sequence of stepsize (for lasso_sg, lasso_pg)
#' @param method lasso optimization function
#' @param k number of fold
#' @param n_it number of iteration for lasso_cd method
#' @param df Degree of freedom. Number of desired variables to be zero.
#' @export

cv.lasso = function(lambdas, x, y, beta0, wt = NULL, ts, method, k = 5, n_it = 10, df = NULL)
{
  data.partition = cv.random.partition(x, y, k = k)
  lambdas = as.matrix(lambdas)
  colnames(lambdas) = "lambda"
  
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

#' random partitioning for cross-validation
#' 
#' @param x predictors
#' @param y response
#' @param k number of fold
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

#' root mean-squared error
#' 
#' @param original original value
#' @param estimate estimated value
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

#' function for normalizing the data
#' 
#' @param x predictor
#' @param y response
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

#' lasso plots
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

#' cv curve plot
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

#' Examples
#' @param wt weights for the coefficients of weighted LASSO.
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
  beta0 = b
  ts = opt_ts(0.1, 1000, 1000)
  
  test.cv = cv.lasso(lambdas, x, y, beta0, wt = wt, ts, method = lasso_cd)
  
  cv.plot(test.cv, main = "Cross-validation error (LASSO using co-ordinate descent)")
  x11()
  lasso.plot(test.cv, main = "LASSO using co-ordinate descent")
  
  cat(sprintf("Least Square Model \n"))
  print(lm(y ~ x - 1))
  
  cat(sprintf("LASSO estimates \n"))
  print(test.cv$coeff)
  
  cat(sprintf("Degrees of freedom w.r.t. lambda \n"))
  print(test.cv$df)
}
