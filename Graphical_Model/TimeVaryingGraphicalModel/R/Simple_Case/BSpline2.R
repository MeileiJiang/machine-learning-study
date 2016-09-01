##########################################################################
## BSpline2.R
## BSpline Across different time on the coefficient for toyexample2.
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

load("R/Simple_Case/2dexample2.RData")
load("rdata/2dtoyexample4.Rdata")

Data = Toydata[,-3]
time = unique(Toydata$Time)
n = length(Toydata$Time)/length(time)

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


w0 = 0.5
w1 = 20
w2 = 1
A = rbind(w0*B, w1*D1B)

Ap = solve(t(A)%*% A ) %*% t(A)

## Set the parameter for plots

set = seq(0, 1, by = 0.005)
# basis = cbind(bs(set, knots = grid, intercept = F, degree = 1), bs(set, knots = grid, intercept = F, degree = 2),
#               bs(set, knots = grid, intercept = F, degree = 3))

basis = cbind(splineDesign(knots = grid, x = set, ord = 1, outer.ok = TRUE), splineDesign(knots = grid, x = set, ord = 2, outer.ok = TRUE), 
              splineDesign(knots = grid, x = set, ord = 3, outer.ok = TRUE), splineDesign(knots = grid, x = set, ord = 4, outer.ok = TRUE))

basis = cbind(splineDesign(knots = grid, x = set, ord = 1, outer.ok = TRUE))

# node 1 ------------------------------------------------------------------

X1 = as.matrix(Data[,-1] )
Y1 = as.matrix(Data[,1])


checkEquals( length(time)*n, dim(X1)[1])
U1 = matrix(nrow = length(time)*n, ncol = dim(X1)[2]*dim(B)[2])
for(t in 1:length(time)){
  U1[ (1+(t-1)*n) : (t*n),] = X1[(1+(t-1)*n) : (t*n),]%*%t(B[t,])
}


# directly solve the model ------------------------------------------------

fit1 = glmnet(x = U1, y = Y1, family = "gaussian", intercept = F)
cv.glmnet(x = U1, y = Y1, family = "gaussian", intercept = F)


gamma1 = as.matrix(coef.glmnet(fit1, s =0.08)[-1,1])

beta12 =  basis %*% gamma1[1:dim(B)[2],1]


beta1 = data.frame(beta12, set)
mbeta1 = melt(beta1, id.vars = "set")

g0 = ggplot(data = mbeta1, aes(x = set, y = value)) +
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
       B-spline basis by lasso")
pdf(file = "figures/2dexample_lasso.pdf")
g0
dev.off()

# solve by fused lasso ---------------------------------------------------------

D0 = diag(dim(U1)[2])
D1 = getDtf(dim(U1)[2], ord =0)
D2 = getDtf(dim(U1)[2], ord =1)
D = rbind(D0, D2)

out1 = fusedlasso1d(y = Y1, X = U1, gamma =  1)
out1 = trendfilter(y = Y1, X = U1, ord = 0.7)

plot(out1)


gamma1 = coef.genlasso(out1, lambda = 9)$beta
eta1 = basis %*% gamma1
plot( eta1, type = "l")


beta12 =  basis %*% gamma1[1:dim(B)[2],1]


beta1 = data.frame(beta12, set)
mbeta1 = melt(beta1, id.vars = "set")

g2 = ggplot(data = mbeta1, aes(x = set, y = value)) +
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
       B-spline basis by fused lasso")

pdf(file = "figures/2dexample_fused_lasso.pdf")
g2
dev.off()



# solve by generalized lasso ---------------------------------------------------------



out1 = genlasso(y = Y1, X = U1, D = A)

plot(out1)


gamma1 = coef.genlasso(out1, lambda = 2.2)$beta
eta1 = basis %*% gamma1
plot( eta1, type = "l")

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

# solve by FLiRTI ---------------------------------------------------------

V1 = U1 %*% Ap
dim(V1)

fit1.0 = glmnet(x = V1, y = Y1, family = "gaussian", intercept = F)
s0 = cv.glmnet(x = V1, y = Y1, family = "gaussian", intercept = F)$lambda.1se
eta1 = as.matrix(coef.glmnet(fit1.0, s = s0)[-1,1])
plot(eta1, type = "l")
gamma1 = Ap %*% eta1

eta1.0 = pmax(1, eta1)

W1 = diag(sqrt(1/c(eta1.0)))

V2 = V1 %*% W1

fit2.0 = glmnet(x = V2, y = Y1, family = "gaussian", intercept = F)

s1 = cv.glmnet(x = V2, y = Y1, family = "gaussian", intercept = F)$lambda.1se

eta2 = W1 %*% as.matrix(coef.glmnet(fit2.0, s = 0.02)[-1,1])
plot(eta2, type = "l")
gamma1 = Ap %*% eta2




# get the coeffiecient function beta(t)

beta12 =  basis %*% gamma1[1:dim(B)[2],1]


beta1 = data.frame(beta12, set)
mbeta1 = melt(beta1, id.vars = "set")

g3 = ggplot(data = mbeta1, aes(x = set, y = value)) +
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
       B-spline basis by reparameterized lasso")

pdf(file = "figures/2dexample_reparameterized_lasso.pdf")
g3
dev.off()

# solve by  elasticnet---------------------------------------------------------

V1 = U1 %*% Ap
dim(V1)

fit1.0 = glmnet(x = V1, y = Y1, family = "gaussian", intercept = F, alpha = 0.5)
s0 = cv.glmnet(x = V1, y = Y1, family = "gaussian", intercept = F, alpha = 0.5)$lambda.1se
eta1 = as.matrix(coef.glmnet(fit1.0, s = 2*s0)[-1,1])
plot(eta1, type = "l")
gamma1 = Ap %*% eta1

beta12 =  basis %*% gamma1[1:dim(B)[2],1]


beta1 = data.frame(beta12, set)
mbeta1 = melt(beta1, id.vars = "set")

g4 = ggplot(data = mbeta1, aes(x = set, y = value)) +
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
       B-spline basis by elastic net")

pdf(file = "figures/2dexample_elastic_net.pdf")
g4
dev.off()
