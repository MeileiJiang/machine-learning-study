##################################################################################
## basis_expansion_penalty.R
## Generate the penalty matrix for basis expansion method of varying coefficient
## model.
## Author: Meilei
##################################################################################
library(splines)
library(dplyr)
library(genlasso)
library(Matrix)

basis_expansion_penalty = function(grid, time, order, p, w0 = 1, w = 0.05){
  # bspline basis matrix
  B  = splineDesign(knots = grid, x = time, ord = order, outer.ok = TRUE)
  # total variation penalty matrix on order-1
  if(order > 1){
    delta = grid[2] - grid[1];
    b1 = splineDesign(knots = grid, x = time, ord = 1, outer.ok = TRUE)
    db = b1 %*% t(getDtf(dim(b1)[2], ord = order-2))/delta^{order-1}
    DB = getD1d(dim(db)[1]) %*% db
  }else{
    DB = getD1d(dim(B)[1]) %*% B
  }
  # penalty matrix
  A0 = rbind(w0*B, w*DB)
  A  = as.matrix(bdiag( replicate(p, A0, simplify = FALSE)))
  return(list("mat" = A, "mat0" = A0))
}