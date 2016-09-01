###############################################################################################
## Total_variation_penalty_order1.R
## Apply total variation penalty on splines.
## Author: Meilei
## Date: May 2016
################################################################################################

library(splines)
library(dplyr)
library(reshape2)
library(ggplot2)
library(genlasso)
library(glmnet)
library(Matrix)
library(RUnit)
require(stats); require(graphics)
source("R/ADMM/ADMM.R")
load("R/Simple_Case/2dexample2.RData")
# visualize the function --------------------------------------------------


ggplot(data = c0.df, aes(x = grid, y = y0)) +
  geom_line(col = "blue") +
  scale_y_continuous(limits = c(-1,1)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "time", y = "partical correlation", title = "Sparse Piecewise Constant")




# load the simulation data ------------------------------------------------

load("rdata/2dtoyexample4.Rdata")
Toydata = Toydata0

Data = Toydata[,-3]
time = unique(Toydata$Time)
time
# [1] 0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30
# [17] 0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50 0.52 0.54 0.56 0.58 0.60 0.62
# [33] 0.64 0.66 0.68 0.70 0.72 0.74 0.76 0.78 0.80 0.82 0.84 0.86 0.88 0.90 0.92 0.94
# [49] 0.96 0.98 1.00

n = length(Toydata$Time)/length(time)
n
# [1] 100

delta = 1/30
grid = seq(0, 1, by = delta)


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


# upto order 3
B = cbind(b1, b2, b3)

d1B = cbind(d1b1, d1b2, d1b3)
d2B = cbind(d2b1, d2b2, d2b3)
d3B = cbind(d3b1, d3b2, d3b3)


# order 3
B = b3;
d1B = d1b3;
d2B = d2b3;
d3B = d3b3;

# order 2
B = b2
d1B = d1b2;
d2B = d2b2;
d3B = d3b2;

#  order 1
B = b1
d1B = d1b1;
d2B = d2b1;
d3B = d3b1;


# set up penalty matrix ---------------------------------------------------

# Total variation penalty
D1 = getD1d(dim(B)[1])
D1B = D1 %*% B
dim(D1B)

Dd1 = getD1d(dim(d1B)[1])
Dd1B = Dd1 %*% d1B

Dd2 = getD1d(dim(d2B)[1])
Dd2B = Dd2 %*% d2B

Dd3 = getD1d(dim(d3B)[1])
Dd3B = Dd3 %*% d3B


# # FLirTi
# 
# D2 = getD1d(dim(D1B)[1])
# D2B = D2 %*% D1B
# dim(D2B)
# # [1] 100  74 
# D3 = getD1d(dim(D2B)[1])
# D3B = D3 %*% D2B
# dim(D3B)

checkEquals(dim(B),dim(d1B)) 
# [1] TRUE
checkEquals(dim(B),dim(d2B)) 
# [1] TRUE

## Set the parameter for plots

set = seq(0, 1, by = 0.005)

# Up to order 3
basis = cbind(splineDesign(knots = grid, x = set, ord = 1, outer.ok = TRUE), splineDesign(knots = grid, x = set, ord = 2, outer.ok = TRUE), 
              splineDesign(knots = grid, x = set, ord = 3, outer.ok = TRUE))

# Order 3
basis = splineDesign(knots = grid, x = set, ord = 3, outer.ok = TRUE)

# order 2
basis = cbind(splineDesign(knots = grid, x = set, ord = 2, outer.ok = TRUE))

# order 1
basis = cbind(splineDesign(knots = grid, x = set, ord = 1, outer.ok = TRUE))

# node 1 ------------------------------------------------------------------

X1 = as.matrix(Data[,-1] )
Y1 = as.matrix(Data[,1])


checkEquals( length(time)*n, dim(X1)[1])
U1 = matrix(nrow = length(time)*n, ncol = dim(X1)[2]*dim(B)[2])
for(t in 1:length(time)){
  U1[ (1+(t-1)*n) : (t*n),] = X1[(1+(t-1)*n) : (t*n),]%*%t(B[t,])
}




# solve by generalized lasso ---------------------------------------------------------
w0 = 1
w1 = 100



# total variation penalty
A = rbind(w0*B, w1*D1B)


# # Flirti
# A = rbind(w0*B, w1*D1B)
# A = rbind(w0*B, w2*D2B)
# A = rbind(w0*B, w3*D3B)
# 
# 

# genlasso
A1 = unique(A)
out1 = genlasso(y = Y1, X = U1, D = A1)
plot(out1)
gamma1 = coef.genlasso(out1, lambda = 10)$beta

# admm_genlasso
GAMMA1 = admm_gen_lasso(Y1, U1, A, epsilon = 10^(-6), lambda = 20, rho = 1, maxstep = 15000, minstep = 200)
gamma_mat = GAMMA1$matrix
plot(gamma_mat[,3], type = "l")

gamma1 = GAMMA1$result


# lasso
wcv = cv.glmnet(x = U1, y = Y1, intercept = F)
lamb = wcv$lambda.1se

wfit = glmnet(x =U1, y = Y1, lambda = 3*lamb, intercept = F)
gamma1 = coef(wfit)[-1,]
# 
# 
# # Fused lasso
# out1.2 = fusedlasso1d(y = Y1, X = U1)
# plot(out1.2)
# 
# gamma1 = coef(out1.2, lambda = 20)$beta


# visualize the result ----------------------------------------------------

plot(gamma1, type = "l")

eta1 = B %*% gamma1
plot( eta1, type = "l")
eta1.1 = B[,1:dim(b1)[2]] %*% gamma1[1:dim(b1)[2]]
plot(eta1.1, type = "l")
eta1.2 = B[,(dim(b1)[2] + 1):(dim(b1)[2] + dim(b2)[2]) ] %*% gamma1[(dim(b1)[2] + 1):(dim(b1)[2] + dim(b2)[2])]
plot(eta1.2, type = "l")
eta1.3 = B[,(dim(b1)[2] + dim(b2)[2] + 1):(dim(b1)[2] + dim(b2)[2] + dim(b3)[2]) ] %*% gamma1[(dim(b1)[2] + dim(b2)[2] + 1):(dim(b1)[2] + dim(b2)[2] + dim(b3)[2])]
plot(eta1.3, type = "l")



deta1 = d1B %*% gamma1
plot( deta1, type = "l")
Deta1 = D1B %*% gamma1
plot( Deta1, type = "l")

deta2 = d2B %*% gamma1
plot( deta2, type = "l")
# Deta2 = D2B %*% gamma1
# plot( Deta2, type = "l")

deta3 = d3B %*% gamma1
plot( deta3, type = "l")
# Deta3 = D3B %*% gamma1
# plot( Deta3, type = "l")

beta1 = data.frame(beta12 = basis %*% gamma1, set)
mbeta0 = melt(beta1, id.vars = "set")
mbeta1 = melt(beta1, id.vars = "set")
mbeta2 = melt(beta1, id.vars = "set")
mbeta3 = melt(beta1, id.vars = "set")

# make the plot -----------------------------------------------------------

## up to order 3
# w1 = 50, lambda = 1
g0 = ggplot(data = mbeta0, aes(x = set, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_line(data = c0.df, aes(x = grid, y = -y0),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = seq(from = -0.6, to = 0.6, by = 0.2),
                     labels = round(seq(from = -0.6, to = 0.6, by = 0.2),1)) +
  labs(x = "time", title = "Estimation of Beta12 by basis expansion 
       B-spline Basis: up to order 3")

g0


# order 1 basis
g1 = ggplot(data = mbeta1, aes(x = set, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_line(data = c0.df, aes(x = grid, y = -y0),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = seq(from = -0.6, to = 0.6, by = 0.2),
                     labels = round(seq(from = -0.6, to = 0.6, by = 0.2),1)) +
  labs(x = "time", title = "Estimation of Beta12 by basis expansion 
       B-spline Basis: order 1")

g1

# order 2 basis
g2 = ggplot(data = mbeta2, aes(x = set, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_line(data = c0.df, aes(x = grid, y = -y0),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = seq(from = -0.6, to = 0.6, by = 0.2),
                     labels = round(seq(from = -0.6, to = 0.6, by = 0.2),1)) +
  labs(x = "time", title = "Estimation of Beta12 by basis expansion 
       B-spline Basis: order 2")

g2

# order 3 basis
g3 = ggplot(data = mbeta3, aes(x = set, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_line(data = c0.df, aes(x = grid, y = -y0),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = seq(from = -0.6, to = 0.6, by = 0.2),
                     labels = round(seq(from = -0.6, to = 0.6, by = 0.2),1)) +
  labs(x = "time", title = "Estimation of Beta12 by basis expansion 
       B-spline Basis: order 3")

g3

pdf(file = "R/Simple_Case/plots/new_TV_order1.pdf", height = 8, width = 8)
grid.arrange(g0 ,g1, g2, g3,  nrow = 2)
dev.off()
