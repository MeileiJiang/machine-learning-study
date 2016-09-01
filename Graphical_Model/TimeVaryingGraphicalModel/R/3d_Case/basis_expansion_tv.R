#############################################################################################
## basis_expansion_tv.R
## author: Meilei Jiang
## Date: July 16th.
#############################################################################################
library(splines)
library(dplyr)
library(reshape2)
library(genlasso)
library(Matrix)
source("R/ADMM/ADMM.R")
source("R/3d_Case/basis_expansion_penalty.R")
# Assume that time points are evenly spaced and each time point contains the same number of samples.
# gamma0 is the initial value
basis_expansion_tv = function(X, Y, time, grid, order = 2, w0 = 1, w = 0.05, epsilon = 10^(-4), lambda = 10, rho = 1, 
                              set = seq(0, 1, by = 0.001), gamma0 = NULL, u0 = NULL, maxstep = 10000, minstep = 100, Print = F){
  # bspline basis matrix
  B  = splineDesign(knots = grid, x = time, ord = order, outer.ok = TRUE)
  
  # basis expansion
  n = dim(X)[1]/length(time) # number of samples on each time point
  p = dim(X)[2] # number of nodes
  
  U = matrix(nrow = length(time)*n, ncol = dim(X1)[2]*dim(B)[2])
  for(t in 1:length(time)){
    diagB = as.matrix(bdiag( replicate(p, B[t,], simplify = FALSE)))
    U[ (1+(t-1)*n) : (t*n),] = X[(1+(t-1)*n) : (t*n),]%*%t(diagB)
  }
  
  
  # penalty matrix
  A  = basis_expansion_penalty(grid, time, order, p, w0, w)$mat
  
  # solve the optimization problem by ADMM algorithm
  GAMMA = admm_gen_lasso(Y, U, A, epsilon = epsilon, lambda = lambda, rho = rho, maxstep = maxstep, 
                         minstep = minstep, beta0 = gamma0, u0 = u0, Print = Print)
  gamma = GAMMA$result
  thread = (abs(gamma) > 50*epsilon)
  gamma_mat = matrix(gamma * thread, ncol = p, byrow = FALSE)
  
  # value of beta
  basis = splineDesign(knots = grid, x = set, ord = order, outer.ok = TRUE)
  beta = data.frame(basis%*% gamma_mat)
  
  # other values for initializing
  u = GAMMA$u
  return(list(gamma = gamma_mat, beta = beta,  u = u))
}
