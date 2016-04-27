#############################################################################
## 3dexample.R
## Three dimensional graphical model. Only node 1 and node2 has interaction
## Author: Meilei
## Date: April 2016
#############################################################################
library(igraph)
library(ggplot2)
library(dplyr)
library(Matrix)
library(mvtnorm)

cfun = function(t){
  if(t >= 0 & t <= 0.2) return(min(1 - 10*abs(t-0.1), 0.5))
  if(t >= 0.4 & t <= 0.6) return(max(-1 + 10*abs(t-0.5), -0.5))
  if(t > 0.8 & t <= 1) return(min(1 - 10*abs(t-0.9), 0.5))
  return(0)
}

grid = seq(0, 1, by = 0.01)
y = rep(0, length(grid))

for(i in 1:length(grid)){
  y[i] = cfun(grid[i])
}

c.df = data.frame( y, grid)

ggplot(data = c.df, aes(x = grid, y = y)) +
  geom_line(col = "blue") +
  scale_y_continuous(limits = c(-1,1)) +
  labs(x = "time", y = "partical correlation", title = "Partial correlation over time")

save(c.df, file = "R/3d_Case/3dexample.RData")

Omegafun = function(t){
  M = matrix(c(1,       cfun(t),   0,
               cfun(t), 1,         0,
               0,       0,         1), 
             ncol = 3, nrow = 3)
  colnames(M) = paste0("X", c(1:3))
  rownames(M) = paste0("X", c(1:3))
  return(M)
}


graphdata = function(t, n = 1){
  Omega = Omegafun(t)
  Sigma = solve(Omega)
  time = rep(t, n)
  data = cbind(rmvnorm(n = n, mean = rep(0, nrow(Sigma)), sigma = Sigma), time)
  return(data)
}

simulation = function(time, n = 1){
  X = matrix(nrow = length(time)*n, ncol = 3 + 1)
  for(t in 1:length(time)){
    X[(1+(t-1)*n):(t*n),] = graphdata(time[t], n)
  }
  return(X)
}

time = seq(0, 1, by = 0.01)
Data = simulation(time = time, n = 100)

colnames(Data) = c(paste0("X",c(1:3)),"Time")
Toydata = data.frame(Data)
save(Toydata, file = "rdata/3dtoyexample1.Rdata")

