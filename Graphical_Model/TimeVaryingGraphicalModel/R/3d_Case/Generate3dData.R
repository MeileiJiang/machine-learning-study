###############################################################################
## Generate3dData.R
## Generate 3-d guassian data.
## Author: Meilei
###############################################################################
library(dplyr)
library(Matrix)
library(mvtnorm)
source("R/3d_Case/3dexample.R")

graphdata = function(t, cfun, n = 1){
  Omega = Omegafun(t, cfun)
  Sigma = solve(Omega)
  time = rep(t, n)
  data = cbind(rmvnorm(n = n, mean = rep(0, nrow(Sigma)), sigma = Sigma), time)
  return(data)
}

simulation = function(time, cfun, n = 1){
  X = matrix(nrow = length(time)*n, ncol = 3 + 1)
  for(t in 1:length(time)){
    X[(1+(t-1)*n):(t*n),] = graphdata(time[t], cfun, n)
  }
  return(X)
}

time = seq(0, 1, by = 0.01)

Data = simulation(time = time, cfun = cfun, n = 100)

colnames(Data) = c(paste0("X",c(1:3)),"Time")
Toydata = data.frame(Data)
save(Toydata, file = "rdata/3dtoyexample1.Rdata")
