##############################################################################################
## GenerateToyExample.R
## Construct a function with piecewise polynomial.
## Author: Meilei
## Date: April 2016
##############################################################################################
library(dplyr)
library(ggplot2)
library(mvtnorm)

source("R/Simple_Case/2dexample2.R")

Omegafun = function(t, cfun){
  M = matrix(c(1,       cfun(t),   
               cfun(t), 1), 
             ncol = 2, nrow = 2)
  colnames(M) = paste0("X", c(1:2))
  rownames(M) = paste0("X", c(1:2))
  return(M)
}

 
graphdata = function(t, cfun, n = 1){
  Omega = Omegafun(t, cfun)
  Sigma = solve(Omega)
  time = rep(t, n)
  data = cbind(rmvnorm(n = n, mean = rep(0, nrow(Sigma)), sigma = Sigma), time)
  return(data)
}

simulation = function(time, cfun, n = 1){
  X = matrix(nrow = length(time)*n, ncol = 2 + 1)
  for(t in 1:length(time)){
    X[(1+(t-1)*n):(t*n),] = graphdata(time[t], cfun, n)
  }
  return(X)
}

time = seq(0, 1, by = 0.01)

Data0 = simulation(time = time, cfun = cfun0, n = 100)
Data1 = simulation(time = time, cfun = cfun1, n = 100)
Data2 = simulation(time = time, cfun = cfun2, n = 100)

colnames(Data0) = c(paste0("X",c(1:2)),"Time")
colnames(Data1) = c(paste0("X",c(1:2)),"Time")
colnames(Data2) = c(paste0("X",c(1:2)),"Time")

Toydata0 = data.frame(Data0)
Toydata1 = data.frame(Data1)
Toydata2 = data.frame(Data2)

save(Toydata0, Toydata1, Toydata2, file = "rdata/2dtoyexample4.Rdata")
