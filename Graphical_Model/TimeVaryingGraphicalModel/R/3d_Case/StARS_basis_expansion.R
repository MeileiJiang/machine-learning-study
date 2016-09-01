#############################################################################
## StARS_basis_expansion.R
## Using Stability selection for tuning parameters.
## Author: Meilei
#############################################################################
library(dplyr)
library(ggplot2)
library(reshape2)
library(glmnet)
library(mvtnorm)
library(huge)
# set up paranell computing
library(doMC)
library(foreach)

registerDoMC(4)
# loading supporting function
source("R/3d_Case/subsample.R")
source("R/3d_Case/basis_expansion_tv.R")
source("R/3d_Case/basis_expansion_penalty.R")

load("rdata/3dtoyexample1.Rdata")

M = 100


time = unique(Toydata$Time)
n = length(Toydata$Time)/length(time)

delta = 0.05
grid = seq(0, 1, by = delta)

pen.mat  = basis_expansion_penalty(grid, time, order = 2, p = 2, w0 = 1, w = 0.05)
A = pen.mat$mat
A0 = pen.mat$mat0

set = seq(0, 1, by = 0.001)

L = c(seq(from = 1, to = 15, by = 2), 20, 50)

stars_result = foreach(i= 1:length(L), .combine = cbind) %dopar%{
  labmda = L[i]
  stars_basis_expansion(lambda = lambda, M = 10)
}

stars_result2 = foreach(i= 1:length(L), .combine = cbind) %dopar%{
  labmda = L[i]
  stars_basis_expansion(lambda = lambda, M = 50)
}

stars_result3 = foreach(i= 1:length(L), .combine = cbind) %dopar%{
  labmda = L[i]
  stars_basis_expansion(lambda = lambda, M = 60)
}
stars_basis_expansion = function(lambda, M = 20, epsilon = 10^(-6)){
  Data = Toydata[,-4]
  X1 = as.matrix(Data[,-1])
  Y1 = as.matrix(Data[,1])
  gamma0 = NULL
  u0 = NULL
  Beta1 = basis_expansion_tv(X = X1, Y = Y1, time = time, grid, order = 2, w = 0.05, set = set, lambda = lambda,
                             rho = 1, epsilon = epsilon, gamma0 = gamma0, u0 = u0, 
                             maxstep = 20000,minstep = 20, Print = F)
  gamma0 = as.vector(Beta1$gamma)
  u0 = Beta1$u
  count = NULL
  support = NULL
  for (k in 1:M) {
    Subtoydata = subsample(Toydata, frac = 0.5)
    Data = Subtoydata[,-4]
    
    # node 1
    X1 = as.matrix(Data[,-1])
    Y1 = as.matrix(Data[,1])
    
    Beta1 = basis_expansion_tv(X = X1, Y = Y1, time = time, grid, order = 2, w = 0.05, set = set, lambda = lambda,
                               rho = 1 , epsilon = epsilon, gamma0 = gamma0, u0 = u0, 
                               maxstep = 10000,minstep = 20, Print = F)
    gamma1 = Beta1$gamma
    if(is.null(count)){
      count = (gamma1 != 0)
    }else{
      count = count + (gamma1 != 0)
    }
    
    beta1 = Beta1$beta
    if(is.null(support)){
      support = (beta1 != 0)
    }else{
      support = support + (beta1 != 0)
    }
  #  if(k %%10 == 0) print(paste0("Finish ", sprintf("%.1f", round(100* k/M,1)), "%"))
  }
  
  count = count/M
  support = support/M
  
  stab_count = 2*count*(1-count)
  stab_count_val = mean(stab_count)
 
  stab_support = 2*support*(1-support)
  stab_support_val = mean(stab_support)
  return(c(stab_count = stab_count_val, stab_support = stab_support_val))
}

colnames(count) = c("beta12", "beta13")
rownames(count) = paste0("f", c(1:dim(count)[1]))
count.df = melt(count)
names(count.df) = c("basis","coef_func","value")

stab_count.df = melt(stab_count)
names(stab_count.df) = c("basis","coef_func","value")
stab_support.df = data.frame(beta12 = stab_support[,1], beta13 = stab_support[,2], time = set)

support.df = data.frame(beta12 = support[,1], beta13 = support[,2], time = set)
ggplot(data = count.df, aes(x = basis, y = coef_func, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "green", high = "blue", midpoint = 0.7)
ggplot(data = stab_count.df, aes(x = basis, y = coef_func, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "white",  high = "blue")

ggplot(data = support.df, aes(x = time, y = beta12)) + geom_point() + ylim(c(0, 1))
ggplot(data = support.df, aes(x = time, y = beta13)) + geom_point() + ylim(c(0, 1))

ggplot(data = stab_support.df, aes(x = time, y = beta12)) + geom_point() + ylim(c(0, 1/2))
ggplot(data = stab_support.df, aes(x = time, y = beta13)) + geom_point() + ylim(c(0, 1/2))

save(count, count.df, stab_count.df, support, support.df, stab_support.df, file = "stars_3d.rdata")


s1 = stars_basis_expansion(1,  M = 50)
s2 = stars_basis_expansion(5,  M = 50)
s3 = stars_basis_expansion(10, M = 50)
s4 = stars_basis_expansion(15, M = 50)
s5 = stars_basis_expansion(20, M = 50)
s6 = stars_basis_expansion(25, M = 50)
s7 = stars_basis_expansion(30, M = 50)
