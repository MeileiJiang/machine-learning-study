##############################################################################################
## 2dexample2.R
## Construct a function with piecewise polynomial.
## Author: Meilei
## Date: April 2016
##############################################################################################


library(igraph)
library(ggplot2)
library(dplyr)
library(Matrix)
source('R/Simple_Case/quad_fun.R')

grid = seq(0, 1, by = 0.001)

# different functions -----------------------------------------------------


# sparse piecewise constant
cfun0 = function(t){
  if(t > 0 & t <= 0.2) return( 0.6)
  if(t >= 0.4 & t <= 0.6) return( - 0.6)
  if(t > 0.8 & t < 1) return( 0.6)
  return(0)
}


# sparse piecewise linear

cfun1 = function(t){
  if(t > 0 & t <= 0.2) return( 0.6 - 6*abs(t-0.1))
  if(t >= 0.4 & t <= 0.6) return( - 0.6 + 6*abs(t-0.5))
  if(t > 0.8 & t < 1) return( 0.6-6*abs(t-0.9))
  return(0)
}




# sparse piecewise quadratic

cfun2 = function(t){
  if(t > 0 & t <= 0.2) return( 0.6*quad_fun(5*t) )
  if(t >= 0.4 & t <= 0.6) return( - 0.6*quad_fun(5*(t-0.4)))
  if(t > 0.8 & t < 1) return( 0.6*quad_fun(5*(t-0.8)))
  return(0)
}



# generate the function line ----------------------------------------------


# y0 = rep(0, length(grid))
# 
# for(i in 1:length(grid)){
#   y0[i] = cfun0(grid[i])
# }
# 
# c0.df = data.frame( y0, grid)
# 
# 
# y1 = rep(0, length(grid))
# for(i in 1:length(grid)){
#   y1[i] = cfun1(grid[i])
# }
# 
# c1.df = data.frame( y1, grid)
# 
# 
# y2 = rep(0, length(grid))
# 
# for(i in 1:length(grid)){
#   y2[i] = cfun2(grid[i])
# }
# 
# c2.df = data.frame( y2, grid)
# 
# 
# # visualize the function --------------------------------------------------
# 
# 
# ggplot(data = c0.df, aes(x = grid, y = y0)) +
#   geom_line(col = "blue") +
#   scale_y_continuous(limits = c(-1,1)) +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
#   labs(x = "time", y = "partical correlation", title = "Sparse Piecewise Constant")
# 
# ggplot(data = c1.df, aes(x = grid, y = y0)) +
#   geom_line(col = "blue") +
#   scale_y_continuous(limits = c(-1,1)) +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
#   labs(x = "time", y = "partical correlation", title = "Sparse Piecewise Linear")
# 
# ggplot(data = c2.df, aes(x = grid, y = y2)) +
#   geom_line(col = "blue") +
#   scale_y_continuous(limits = c(-1,1)) +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
#   labs(x = "time", y = "partical correlation", title = "Sparse Piecewise Quadratic")
# 
# save(c0.df, c1.df, c2.df, file = "R/Simple_Case/2dexample2.RData")

