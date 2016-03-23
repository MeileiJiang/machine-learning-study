#####################################################################################
## B-spline.R
## Understand how to utilize B-spline.
## Author: Meilei
#####################################################################################
library(splines)
library(dplyr)
library(reshape2)
library(ggplot2)
library(genlasso)
library(glmnet)
library(Matrix)
library(RUnit)
require(stats); require(graphics)

load("rdata/toyexample4.Rdata")

Data = Toydata[,-11]
time = unique(Toydata$Time)
n = length(Toydata$Time)/length(time)

grid = seq(0, 1, by = 0.1)
S = 0.06
L = 65
# grid = time

# bspline basis matrix
B = cbind(splineDesign(knots = grid, x = time, ord = 1, outer.ok = TRUE), splineDesign(knots = grid, x = time, ord = 2,outer.ok = TRUE), 
      splineDesign(knots = grid, x = time, ord = 3, outer.ok = TRUE), splineDesign(knots = grid, x = time, ord = 4,outer.ok = TRUE))

B = cbind(bs(time, knots = grid, degree = 1, intercept = F), bs(time, knots = grid, degree = 2, intercept = F),
          bs(time, knots = grid, degree = 3, intercept = F))
dim(B)

# First order difference operator 
DifferenceOperator = function(B, time){
  DB = matrix(nrow = dim(B)[1]-1, ncol = dim(B)[2])
  for(t in 1:dim(DB)[1]){
    DB[t,] = (B[t+1,] - B[t])/(time[t+1] - time[t])
  }
  return(DB)
} 
D1B = DifferenceOperator(B, time)
D2B = DifferenceOperator(D1B, time)
D3B = DifferenceOperator(D2B, time)

Rowbdiag = function(M){
  newM = NULL
  for(t in 1:dim(M)[1]){
    tempM = as.matrix(t(bdiag(M[t,], M[t,], M[t,], 
                              M[t,], M[t,], M[t,], 
                              M[t,], M[t,], M[t,])))
    newM = rbind(newM, tempM)
  }
  return(newM)
}

A = Rowbdiag(B)
# A = Rowbdiag(rbind(B, D3B))
set = seq(0, 1, by = 0.005)
basis = cbind(bs(set, knots = grid, intercept = F, degree = 1), bs(set, knots = grid, intercept = F, degree = 2),
              bs(set, knots = grid, intercept = F, degree = 3))

basis = cbind(splineDesign(knots = grid, x = set, ord = 1), splineDesign(knots = set, x = time, ord = 2,outer.ok = TRUE), 
              splineDesign(knots = grid, x = set, ord = 3, outer.ok = TRUE), splineDesign(knots = set, x = time, ord = 4,outer.ok = TRUE))

# node 1 ------------------------------------------------------------------

X1 = as.matrix(Data[,-1] )
Y1 = as.matrix(Data[,1])


checkEquals( length(time)*n, dim(X1)[1])
U1 = matrix(nrow = length(time)*n, ncol = dim(X1)[2]*dim(B)[2])
for(t in 1:length(time)){
  diagB = as.matrix(bdiag(B[t,], B[t,], B[t,], B[t,], B[t,],
                B[t,], B[t,], B[t,], B[t,]))
  U1[ (1+(t-1)*n) : (t*n),] = X1[(1+(t-1)*n) : (t*n),]%*% t(diagB)
}


# directly solve the model ------------------------------------------------

fit1 = glmnet(x = U1, y = Y1, family = "gaussian", intercept = F)

gamma1 = as.matrix(coef.glmnet(fit1, s = S)[-1,1])


out1 = genlasso(y = Y1, X = U1, D = A)
plot(out1)

# get the estimation of gamma -- coefficients of basis functions for beta(t) -------------------------------------------

gamma1 = coef.genlasso(out, lambda = L)$beta

# get the coeffiecient function beta(t)

beta12 =  basis %*% gamma1[1:dim(B)[2],1]
beta13 =  basis %*% gamma1[(dim(B)[2]+1):(2*dim(B)[2]),1]
beta14 =  basis %*% gamma1[(2*dim(B)[2]+1):(3*dim(B)[2]),1]
beta15 =  basis %*% gamma1[(3*dim(B)[2]+1):(4*dim(B)[2]),1]
beta16 =  basis %*% gamma1[(4*dim(B)[2]+1):(5*dim(B)[2]),1]
beta17 =  basis %*% gamma1[(5*dim(B)[2]+1):(6*dim(B)[2]),1]
beta18 =  basis %*% gamma1[(6*dim(B)[2]+1):(7*dim(B)[2]),1]
beta19 =  basis %*% gamma1[(7*dim(B)[2]+1):(8*dim(B)[2]),1]
beta110 =  basis %*% gamma1[(8*dim(B)[2]+1):(9*dim(B)[2]),1]

beta1 = data.frame(beta12, beta13, beta14, beta15, beta16, 
                   beta17, beta18, beta19, beta110, set)
mbeta1 = melt(beta1, id.vars = "set")

g1 = ggplot(data = mbeta1, aes(x = set, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.15),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0),
                     labels = c(-0.5,-0.4,-0.3,-0.2,-0.1,0)) +
  labs(x = "time", title = "Coefficient of Node 1 by basis expansion (B-spline)")


# node 2 ------------------------------------------------------------------


X2 = as.matrix(Data[,-2] )
Y2 = as.matrix(Data[,2])




checkEquals( length(time)*n, dim(X2)[1])
U2 = matrix(nrow = length(time)*n, ncol = dim(X2)[2]*dim(B)[2])
for(t in 1:length(time)){
  diagB = as.matrix(bdiag(B[t,], B[t,], B[t,], B[t,], B[t,],
                          B[t,], B[t,], B[t,], B[t,]))
  U2[ (1+(t-1)*n) : (t*n),] = X2[(1+(t-1)*n) : (t*n),]%*% t(diagB)
}


# directly solve the model ------------------------------------------------

fit2 = glmnet(x = U2, y = Y2, family = "gaussian", intercept = F)

gamma2 = as.matrix(coef.glmnet(fit2, s = S)[-1,1])

# First order difference operator 

out2 = genlasso(y = Y2, X = U2, D = A)
plot(out2)

# get the estimation of gamma -- coefficients of basis functions for beta(t) -------------------------------------------

gamma2 = coef.genlasso(out2, lambda = 70)$beta

# get the coeffiecient function beta(t)


beta21 =  basis %*% gamma2[1:dim(B)[2],1]
beta23 =  basis %*% gamma2[(dim(B)[2]+1):(2*dim(B)[2]),1]
beta24 =  basis %*% gamma2[(2*dim(B)[2]+1):(3*dim(B)[2]),1]
beta25 =  basis %*% gamma2[(3*dim(B)[2]+1):(4*dim(B)[2]),1]
beta26 =  basis %*% gamma2[(4*dim(B)[2]+1):(5*dim(B)[2]),1]
beta27 =  basis %*% gamma2[(5*dim(B)[2]+1):(6*dim(B)[2]),1]
beta28 =  basis %*% gamma2[(6*dim(B)[2]+1):(7*dim(B)[2]),1]
beta29 =  basis %*% gamma2[(7*dim(B)[2]+1):(8*dim(B)[2]),1]
beta210 =  basis %*% gamma2[(8*dim(B)[2]+1):(9*dim(B)[2]),1]

beta2 = data.frame(beta21, beta23, beta24, beta25, beta26, 
                   beta27, beta28, beta29, beta210, set)
mbeta2 = melt(beta2, id.vars = "set")

g2 = ggplot(data = mbeta2, aes(x = set, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.15),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0),
                     labels = c(-0.5,-0.4,-0.3,-0.2,-0.1,0)) +
  labs(x = "time", title = "Coefficient of Node 2 by basis expansion (B-spline)")




# node 3 ------------------------------------------------------------------

X3 = as.matrix(Data[,-3] )
Y3 = as.matrix(Data[,3])


checkEquals( length(time)*n, dim(X3)[1])
U3 = matrix(nrow = length(time)*n, ncol = dim(X3)[2]*dim(B)[2])
for(t in 1:length(time)){
  diagB = as.matrix(bdiag(B[t,], B[t,], B[t,], B[t,], B[t,],
                          B[t,], B[t,], B[t,], B[t,]))
  U3[ (1+(t-1)*n) : (t*n),] = X3[(1+(t-1)*n) : (t*n),]%*% t(diagB)
}


# directly solve the model ------------------------------------------------

fit3 = glmnet(x = U3, y = Y3, family = "gaussian", intercept = F)

gamma3 = as.matrix(coef.glmnet(fit3, s = S)[-1,1])

# First order difference operator 

out3 = genlasso(y = Y3, X = U3, D = A)
plot(out3)

# get the estimation of gamma -- coefficients of basis functions for beta(t) -------------------------------------------

gamma3 = coef.genlasso(out3, lambda = 75)$beta

# get the coeffiecient function beta(t)


beta31 =  basis %*% gamma3[1:dim(B)[2],1]
beta32 =  basis %*% gamma3[(dim(B)[2]+1):(2*dim(B)[2]),1]
beta34 =  basis %*% gamma3[(2*dim(B)[2]+1):(3*dim(B)[2]),1]
beta35 =  basis %*% gamma3[(3*dim(B)[2]+1):(4*dim(B)[2]),1]
beta36 =  basis %*% gamma3[(4*dim(B)[2]+1):(5*dim(B)[2]),1]
beta37 =  basis %*% gamma3[(5*dim(B)[2]+1):(6*dim(B)[2]),1]
beta38 =  basis %*% gamma3[(6*dim(B)[2]+1):(7*dim(B)[2]),1]
beta39 =  basis %*% gamma3[(7*dim(B)[2]+1):(8*dim(B)[2]),1]
beta310 =  basis %*% gamma3[(8*dim(B)[2]+1):(9*dim(B)[2]),1]

beta3 = data.frame(beta31, beta32, beta34, beta35, beta36, 
                   beta37, beta38, beta39, beta310, set)
mbeta3 = melt(beta3, id.vars = "set")

g3 = ggplot(data = mbeta3, aes(x = set, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.15),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0),
                     labels = c(-0.5,-0.4,-0.3,-0.2,-0.1,0)) +
  labs(x = "time", title = "Coefficient of Node 3 by basis expansion (B-spline)")

# node 6 ------------------------------------------------------------------

X6 = as.matrix(Data[,-6])
Y6 = as.matrix(Data[,6])


checkEquals( length(time)*n, dim(X6)[1])
U6 = matrix(nrow = length(time)*n, ncol = dim(X6)[2]*dim(B)[2])
for(t in 1:length(time)){
  diagB = as.matrix(bdiag(B[t,], B[t,], B[t,], B[t,], B[t,],
                          B[t,], B[t,], B[t,], B[t,]))
  U6[ (1+(t-1)*n) : (t*n),] = X6[(1+(t-1)*n) : (t*n),]%*% t(diagB)
}


# directly solve the model ------------------------------------------------

fit6 = glmnet(x = U6, y = Y6, family = "gaussian", intercept = F)

gamma6 = as.matrix(coef.glmnet(fit6, s = S)[-1,1])

# First order difference operator 

out6 = genlasso(y = Y6, X = U6, D = A)
plot(out6)

# get the estimation of gamma -- coefficients of basis functions for beta(t) -------------------------------------------

gamma6 = coef.genlasso(out6, lambda = 210)$beta

# get the coeffiecient function beta(t)


beta61 =  basis %*% gamma6[1:dim(B)[2],1]
beta62 =  basis %*% gamma6[(dim(B)[2]+1):(2*dim(B)[2]),1]
beta63 =  basis %*% gamma6[(2*dim(B)[2]+1):(3*dim(B)[2]),1]
beta64 =  basis %*% gamma6[(3*dim(B)[2]+1):(4*dim(B)[2]),1]
beta65 =  basis %*% gamma6[(4*dim(B)[2]+1):(5*dim(B)[2]),1]
beta67 =  basis %*% gamma6[(5*dim(B)[2]+1):(6*dim(B)[2]),1]
beta68 =  basis %*% gamma6[(6*dim(B)[2]+1):(7*dim(B)[2]),1]
beta69 =  basis %*% gamma6[(7*dim(B)[2]+1):(8*dim(B)[2]),1]
beta610 =  basis %*% gamma6[(8*dim(B)[2]+1):(9*dim(B)[2]),1]

beta6 = data.frame(beta61, beta62, beta63, beta64, beta65, 
                   beta67, beta68, beta69, beta610, set)
mbeta6 = melt(beta6, id.vars = "set")

g6 = ggplot(data = mbeta6, aes(x = set, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.15),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0),
                     labels = c(-0.5,-0.4,-0.3,-0.2,-0.1,0)) +
  labs(x = "time", title = "Coefficient of Node 6 by basis expansion (B-spline)")


# save the plot -----------------------------------------------------------


pdf(file = "figures/b_spline.pdf")
print(g1)
print(g2)
print(g3)
print(g6)
dev.off()

pdf(file = "figures/b_spline2.pdf")
print(g1)
print(g2)
print(g3)
print(g6)
dev.off()
