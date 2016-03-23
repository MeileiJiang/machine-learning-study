#############################################################################################
## SimluationData.R
## Generating data from time varying graph. Namely, X(t) ~ N(0 , Sigma(t)). Each time point 
## we have one relization.
## Author: Meilei
#############################################################################################
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(glmnet) # Lasso and Elastic Net
library(splines) # B-spline

source("R/Visualize_Chain_Graph.R")

graphdata = function(t, n = 1){
  Omega = OmegaFun(t)
  Sigma = solve(Omega)
  time = rep(t, n)
  data = cbind(rmvnorm(n = n, mean = rep(0, nrow(Sigma)), sigma = Sigma), time)
  return(data)
}

# X = rbind( graphdata(0.1),  graphdata(0.3), graphdata(0.5), graphdata(0.7), graphdata(0.9))
# colnames(X) = paste0("X", c(1:10))
# 
# Data = data.frame(X)
# save(Data, file = "rdata/toyexample1.Rdata")
# 
# lasso = glmnet(x = X[,-1], y = X[,1], family = "gaussian")
# plot.glmnet(lasso)
# coef(lasso)
# cv.lasso = cv.glmnet(x = X[,-1], y = X[,1], family = "gaussian")


simulation = function(time, n = 1){
  X = matrix(nrow = length(time)*n, ncol = 10 + 1)
  for(t in 1:length(time)){
    X[(1+(t-1)*n):(t*n),] = graphdata(time[t], n)
  }
  return(X)
}


