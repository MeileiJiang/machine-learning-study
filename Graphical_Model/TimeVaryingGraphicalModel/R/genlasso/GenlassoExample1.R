###############################################################################################
## GenlassoExample1.R
## Construct a signal approximation example for genlasso.
## Author: Meilei
## Date: April 2016
##############################################################################################
library(igraph)
library(ggplot2)
library(dplyr)
library(Matrix)
library(mvtnorm)

cfun = function(t){
  if(t >= 0 & t <= 0.2) return( 0.6*sin(10*pi*t))
  if(t >= 0.4 & t <= 0.6) return( - 0.6*sin(10*pi*t))
  if(t > 0.8 & t <= 1) return( 0.6*sin(10*pi*t))
  return(0)
}

grid = seq(0, 1, by = 0.01)
y = rep(0, length(grid))
beta = rep(0, length(grid))


for(i in 1:length(grid)){
  y[i] = cfun(grid[i]) + rnorm(1, 0 , 0.09)
  beta[i] = cfun(grid[i])
}

c.df = data.frame( y, beta, grid)

ggplot(data = c.df, aes(x = grid, y = y)) +
  geom_line(col = "blue") +
  geom_line(aes(x =grid, y = beta),col = "red", linetype = "dashed") +
  scale_y_continuous(limits = c(-1,1)) +
  labs(x = "time", y = "partical correlation", title = "Partial correlation over time")

save(c.df, file = "R/genlasso/genlassoExample.RData")
