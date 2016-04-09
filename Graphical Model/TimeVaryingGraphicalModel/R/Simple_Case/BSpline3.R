##########################################################################
## BSpline3.R
## BSpline Across different time on the coefficient. Refine the direvatve matrix
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

load("rdata/2dtoyexample4.Rdata")
load("R/Simple_Case/2dexample2.RData")

Data = Toydata[,-3]
time = unique(Toydata$Time)
n = length(Toydata$Time)/length(time)

grid = seq(0, 1, by = 0.05)
delta = 0.05

S = 1
L = 65
# grid = time

# bspline basis matrix
b1 = splineDesign(knots = grid, x = time, ord = 1, outer.ok = TRUE)
d1b1 = matrix(0, nrow = dim(b1)[1], ncol = dim(b1)[2])
d2b1 = d1b1
d3b1 = d2b1

b2 = splineDesign(knots = grid, x = time, ord = 2, outer.ok = TRUE)
d1b2 = (b1[,-dim(b1)[2]] - b1[,-1])/delta
d2b2 = matrix(0, nrow = dim(b2)[1], ncol = dim(b2)[2])
d3b2 = d2b2

b3 = splineDesign(knots = grid, x = time, ord = 3, outer.ok = TRUE)
d1b3 = (b2[,-dim(b2)[2]] - b2[,-1])/delta
d2b3 = (d1b2[,-dim(d1b2)[2]] - d1b2[,-1])/delta
d3b3 = matrix(0, nrow = dim(b3)[1], ncol = dim(b3)[2])

b4 = splineDesign(knots = grid, x = time, ord = 4, outer.ok = TRUE)
d1b4 = (b3[,-dim(b3)[2]] - b3[,-1])/delta
d2b4 = (d1b3[,-dim(d1b3)[2]] - d1b3[,-1])/delta
d3b4 = (d2b3[,-dim(d2b3)[2]] - d2b3[,-1])/delta


B = cbind(b1, b2, b3, b4)

d1B = cbind(d1b1, d1b2, d1b3, d1b4)
d2B = cbind(d2b1, d2b2, d2b3, d2b4)
d3B = unique(cbind(d3b1, d3b2, d3b3, d3b4))


# B = cbind(splineDesign(knots = grid, x = time, ord = 1, outer.ok = TRUE)) 

checkEquals(dim(B),dim(d1B)) 
checkEquals(dim(B),dim(d2B)) 
checkEquals(dim(B),dim(d3B)) 

# # First order difference operator 
# DifferenceOperator = function(B, time){
#   DB = matrix(nrow = dim(B)[1]-1, ncol = dim(B)[2])
#   for(t in 1:dim(DB)[1]){
#     DB[t,] = (B[t+1,] - B[t,])#/(time[t+1] - time[t])
#   }
#   return(DB)
# } 
# D1B = DifferenceOperator(B, time)
# D2B = DifferenceOperator(D1B, time)
# D3B = DifferenceOperator(D2B, time)


## Set the parameter for plots

set = seq(0, 1, by = 0.005)
# basis = cbind(bs(set, knots = grid, intercept = F, degree = 1), bs(set, knots = grid, intercept = F, degree = 2),
#               bs(set, knots = grid, intercept = F, degree = 3))

basis = cbind(splineDesign(knots = grid, x = set, ord = 1, outer.ok = TRUE), splineDesign(knots = grid, x = set, ord = 2, outer.ok = TRUE), 
              splineDesign(knots = grid, x = set, ord = 3, outer.ok = TRUE), splineDesign(knots = grid, x = set, ord = 4, outer.ok = TRUE))

# basis = cbind(splineDesign(knots = grid, x = set, ord = 1, outer.ok = TRUE))

# node 1 ------------------------------------------------------------------

X1 = as.matrix(Data[,-1] )
Y1 = as.matrix(Data[,1])


checkEquals( length(time)*n, dim(X1)[1])
U1 = matrix(nrow = length(time)*n, ncol = dim(X1)[2]*dim(B)[2])
for(t in 1:length(time)){
  U1[ (1+(t-1)*n) : (t*n),] = X1[(1+(t-1)*n) : (t*n),]%*%t(B[t,])
}




# solve by generalized lasso ---------------------------------------------------------
w0 = 0.5
w1 = 0.01
w2 = 1
w3 = 0.1

A = rbind(w0*B, w1*d1B, w2*d2B)
A = rbind(w0*B, w1*d1B)
A = rbind(w0*B, w3*d3B)

# Ap = solve(t(A)%*% A ) %*% t(A)


out1 = genlasso(y = Y1, X = U1, D = A)

plot(out1)


gamma1 = coef.genlasso(out1, lambda = 3.3)$beta
eta1 = basis %*% gamma1
plot( eta1, type = "l")

eta1.0 = B%*%gamma1
plot( eta1.0, type = "l")
eta1.1 = d1B%*%gamma1
plot( eta1.1, type = "l")
eta1.2 = d2B%*%gamma1
plot( eta1.2, type = "l")
eta1.3 = d3B%*%gamma1
plot(eta1.3, type = "l")

beta12 =  basis %*% gamma1[1:dim(B)[2],1]


beta1 = data.frame(beta12, set)
mbeta1 = melt(beta1, id.vars = "set")

g1 = ggplot(data = mbeta1, aes(x = set, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Functional Coefficient Estimation of Node 1 \n
       B-spline basis by generalized lasso")

pdf(file = "figures/2dexample_genlasso.pdf")
g1
dev.off()

