##################################################################################
## ADMM.R
## Alternative Direction Methods of Multiplier approach to solve generalized lasso.
## argmin_{\beta} 1/2||Y - X \beta ||_2^2 + \lambda ||\alpha||_1 
## subject to     \alpha = A \beta
## Author: Meilei
##################################################################################
require(Matrix)
require(MASS)

source('R/ADMM/Soft_Thresholding_Operator.R')
admm_gen_lasso = function(Y, X, A, epsilon = 10^(-6), lambda, rho, maxstep, minstep = 100, 
                          beta0 = NULL, u0 = NULL, Print = FALSE){
  ## initialization
  if(length(Y) != dim(X)[1]) stop("dimensions of Y and X do not match!")
  p = dim(X)[2]
  if(p != dim(A)[2]) stop("dimensions of X and A do not match!")
  if(is.null(beta0)) beta0 = rep(0, p)
  if(is.null(u0)) u0 = 0
  alpha0 = A %*% beta0
  beta = matrix(nrow = maxstep, ncol = length(beta0))
  alpha = matrix(nrow = maxstep, ncol = length(alpha0))
  u = matrix(nrow = maxstep, ncol = length(alpha0))
  beta[1,] = beta0
  alpha[1,] = alpha0
  u[1,] = u0
  ## ADMM
  for (k in 2:maxstep) {
    beta[k,] = ginv(t(X)%*% X + rho * t(A) %*% A) %*% (t(X) %*% Y + rho*t(A)%*% (alpha[k-1,] + u[k-1,]))
    alpha[k,] = soft_threshold(A %*% beta[k-1,] - u[k-1,], lambda/rho)
    u[k,] = u[k-1,] + rho * (alpha[k-1,] - A %*% beta[k-1,])
    if(k > minstep & sum(abs(beta[k,] - beta[k-1,]))/p < epsilon & sum(abs(u[k,] - u[k-1,]))/p < epsilon ) break
  }
  if(Print)print(paste0("Stop at: ",k))
  return(list(matrix = beta, result = beta[k,], u = u[k,]))
}

# Example

# X = cbind(rep(1,1000), rbind(diag(1, nrow  = 200), diag(2, nrow = 200), 
#                              diag(3, nrow  = 200), diag(4, nrow = 200), diag(-1, nrow  = 200)))
# Beta = c(0 ,rep(c(1,-1,2,-2,0), rep(40,5)))
# 
# Y = X %*% Beta + rnorm(1000)
# 
# 
# A = rbind(10 * getDtf(201, 1),  diag(1, 201))
# 
# hBeta = admm_gen_lasso(Y, X, A, epsilon = 10^(-5), lambda = 5, rho = 1, maxstep = 2000)
# plot(hBeta[,10], type = "l")
# plot(hBeta[,50], type = "l")
# plot(hBeta[1717,], type = "l")
