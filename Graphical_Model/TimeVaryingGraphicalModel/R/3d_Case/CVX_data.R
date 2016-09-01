##################################################################################
## CVX_data.R
## Prepare .mat for solving generalized lasso by CVX in Matlab.
## Author: Meilei Jiang
##################################################################################
library(R.matlab)
library(dplyr)
library(splines)

source("R/3d_Case/basis_expansion_penalty.R")

load("rdata/3dtoyexample1.Rdata")

Data = Toydata[,-4]
time = unique(Toydata$Time)
n = length(Toydata$Time)/length(time)

delta = 0.05
grid = seq(0, 1, by = delta)
order = 2

# node 1
X = as.matrix(Data[,-1])
Y = as.matrix(Data[,1])

# bspline basis matrix
B  = splineDesign(knots = grid, x = time, ord = order, outer.ok = TRUE)

# basis expansion
n = dim(X)[1]/length(time) # number of samples on each time point
p = dim(X)[2] # number of nodes

U = matrix(nrow = length(time)*n, ncol = dim(X)[2]*dim(B)[2])
for(t in 1:length(time)){
  diagB = as.matrix(bdiag( replicate(p, B[t,], simplify = FALSE)))
  U[ (1+(t-1)*n) : (t*n),] = X[(1+(t-1)*n) : (t*n),]%*%t(diagB)
}


# penalty matrix
A  = basis_expansion_penalty(grid, time, order, p, w0, w)$mat

writeMat("R/3d_Case/temp_rdata/Genlasso_cvx.mat", Y = Y, X = X, U = U, A = A)
writeMat("matlab/mat_data/3dtoyexample1.mat", Toydata = Toydata)
