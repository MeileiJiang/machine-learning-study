##########################################################################
## Genlasso_BSpline.R
## BSpline Across different time on the coefficient for GenlassoExample1.
## Author: Meilei
## Date: April 2016
#########################################################################
library(splines)
library(dplyr)
library(reshape2)
library(ggplot2)
library(genlasso)
library(glmnet)
library(Matrix)
library(RUnit)
require(stats); require(graphics)

load("R/genlasso/genlassoExample.RData")


time = c.df$grid

grid = seq(0, 1, by = 0.05)

# grid = time

# bspline basis matrix
B = cbind(splineDesign(knots = grid, x = time, ord = 1, outer.ok = TRUE), 
          splineDesign(knots = grid, x = time, ord = 2, outer.ok = TRUE), 
          splineDesign(knots = grid, x = time, ord = 3, outer.ok = TRUE), 
          splineDesign(knots = grid, x = time, ord = 4, outer.ok = TRUE))

# dB = cbind(splineDesign(knots = grid, x = time, ord = 1, outer.ok = TRUE), 
#            splineDesign(knots = grid, x = time, ord = 2, outer.ok = TRUE), 
#            splineDesign(knots = grid, x = time, ord = 3, outer.ok = TRUE))

# B = cbind(splineDesign(knots = grid, x = time, ord = 1, outer.ok = TRUE)) 

dim(B)

# First order difference operator 
DifferenceOperator = function(B, time){
  DB = matrix(nrow = dim(B)[1]-1, ncol = dim(B)[2])
  for(t in 1:dim(DB)[1]){
    DB[t,] = (B[t+1,] - B[t])#/(time[t+1] - time[t])
  }
  return(DB)
} 
D1B = DifferenceOperator(B, time)
D2B = DifferenceOperator(D1B, time)
D3B = DifferenceOperator(D2B, time)

Bp = solve(t(B)%*% B ) %*% t(B)


w0 = 0.1
w1 = 20
w2 = 20
A = rbind(w0*B, w2*D2B)


## Set the parameter for plots

set = seq(0, 1, by = 0.005)
# basis = cbind(bs(set, knots = grid, intercept = F, degree = 1), bs(set, knots = grid, intercept = F, degree = 2),
#               bs(set, knots = grid, intercept = F, degree = 3))

basis = cbind(splineDesign(knots = grid, x = set, ord = 1, outer.ok = TRUE), splineDesign(knots = grid, x = set, ord = 2, outer.ok = TRUE), 
              splineDesign(knots = grid, x = set, ord = 3, outer.ok = TRUE), splineDesign(knots = grid, x = set, ord = 4, outer.ok = TRUE))

# basis = cbind(splineDesign(knots = grid, x = set, ord = 1, outer.ok = TRUE))

# node 1 ------------------------------------------------------------------

X1 = diag(dim(B)[2])
Y1 = c.df$y


# solve by fused lasso ---------------------------------------------------------

out1 = fusedlasso1d(y = Y1, X = B, gamma =  1)

plot(out1)

gamma1 = coef.genlasso(out1, lambda = 0.32)$beta
eta1 = basis %*% gamma1
plot( eta1, type = "l")


# solve by trend filter ---------------------------------------------------


out2 = trendfilter(y = Y1, X = B, ord = 1)
plot(out2, nlam = 5)

gamma2 = coef(out2, lambda = 2)$beta
eta2 = basis %*% gamma2
plot( eta2, type = "l")



# solve by generalized lasso ---------------------------------------------------------

out3= genlasso(y = Y1, X = B,  D = A)
plot(out3)


gamma3 = coef.genlasso(out1, lambda = 0.18)$beta
eta3 = basis %*% gamma3
plot( eta3, type = "l")


