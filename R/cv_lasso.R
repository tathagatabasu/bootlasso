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
#' @param ts sequence of stepsize (for lasso_sg, lasso_pg)
#' @param method lasso optimization function
#' @param k number of fold
#' @param n_it number of iteration for lasso_cd method
#' @param rel_er Relative error. Should be between 0 (minimum error) and 1 (maximum error).
#' @param df Maximum number of desired (non-zero) variable in the model.
#' @export

cv.lasso = function(lambdas, x, y, beta0, ts, method, k = 5, n_it = 10, rel_er = 0, df = NULL)
{
  data.partition = cv.random.partition(x, y, k = k)
  lambdas = as.matrix(lambdas)
  colnames(lambdas) = "lambda"
  
  # function for training
  
  cv.train = function (i)
  {
    traindata = as.matrix(apply(data.partition[,,-i], 2, rbind))
    x = traindata[,2:ncol(traindata)]
    y = traindata[,1]
    model = method(lambdas, x, y, beta0, ts, n_it = n_it)
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
  cv.error.index = max(which(me <= (min(me) + (max(me) - min(me)) * rel_er)))
  
  nvar = as.matrix(colSums(beta != 0)); colnames(nvar) = "nvar"
  df = ncol(x) - nvar
  er = as.matrix(me); colnames(er) = "error"
  cv.df = cbind(lambdas, nvar, er)
  
  s = min(which(nvar == 0))
  
  if(cv.error.index > s)
    cv.error.index = s
  
  cv.coef = m[,cv.error.index]
  
  if(is.null(df)==F)
  {
    colSums(m[-1,]!=0)
    cv.coef = m[-1,which(colMeans(e)==min(colMeans(e)[which(nvar <= (ncol(x) - df))]))]
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
#' @param output model obtained by doing cross-validation.
#' @param index index obtained from cross-validation for checking desired lambdas intersect.
#' @param main Title of the plot.
#' @export

lasso.plot = function(output, index, main=NULL)
{
  y.lim.up = max(abs(output[2:nrow(output),]))
  y.lim.bel = min(output[2:nrow(output),])
  x = output[1,]
  
  matplot(log(output[1,]), t(output[2:nrow(output),]), type = "l", 
          xlab = expression(paste(log(lambdas))),
          ylab = "Coefficients", ylim = c(y.lim.bel, y.lim.up),
          main = main,
          lty = 1, col = 1:6)
  abline(h = 0, col = "black", lty = 2)
  abline(v = log(x[max(index)]), lty = 2)
}

#' cv curve plot
#' @param error error obtained from cross-validation
#' @param lambdas sequence of lambdas
#' @param index index obtained from cross-validation for checking desired lambdas intersect.
#' @param main Title of the plot.
#' @export

cv.plot = function(error, lambdas, index, main=NULL)
{
  error = as.matrix(error)
  y = colMeans(error)
  ysd = apply(error, 2, sd)
  yHigh = y + ysd / 2
  yLow = y - ysd / 2
  
  x = log(lambdas)
  xHigh = x
  xLow = x
  
  plot(x, y, xlab = expression(log(lambdas)),
       ylab="mean-squared error",
       main=main, 
       ylim = c(min(yLow), max(yHigh)), col = "red", pch = 20)
  
  arrows(xHigh, yHigh, xLow, yLow, col = "grey", angle = 90, length = 0.03, code = 3)
  lines(x, y, col = "red", lty = 2)
  abline(v = x[max(index)], lty = 2)
}

###################################################################################################

#' Examples
#' @param rel_er Relative error. Should be between 0 (minimum error) and 1 (maximum error).
#' @export
example.cv = function(rel_er)
{
  cat(sprintf("\nSimulated Dataset: 24 predictors and 50 observations, no of fold is 5.\n\n"))
  x = matrix(data = rnorm(1200), nrow = 50, ncol = 24)
  b = as.matrix(rep(c(-3,-2,-1,1,2,3), 4))
  er = as.matrix(rnorm(50))
  y = x %*% b + er
  
  lambdas = exp(seq(-5,3,0.1))
  beta0 = b
  ts = opt_ts(0.1, 1000, 1000)
  
  test.cv = cv.lasso(lambdas, x, y, beta0, ts, method = lasso_cd, rel_er = rel_er)
  m = test.cv$model
  e = test.cv$error
  i = test.cv$index
  cv.plot(e, m[1,], i, main = "Cross-validation error (LASSO using co-ordinate descent)")
  lasso.plot(m, i, main = "LASSO using co-ordinate descent")
  
  cat(sprintf("Least Square Model \n"))
  print(lm(y ~ x - 1))
  
  cat(sprintf("LASSO estimates \n"))
  print(test.cv$coeff)
  
  cat(sprintf("No of variables w.r.t. lambda \n"))
  print(test.cv$df)
}
