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
  if(t >= 0.2 & t <= 0.3) return(0.6 - 12*abs(t-0.25))
  if(t >= 0.5 & t <= 0.6) return(-0.6 + 12*abs(t-0.55))
  if(t > 0.8 & t <= 0.9) return(0.6 - 12*abs(t-0.85))
  return(0)
}

grid = seq(0, 1, by = 0.01)
y = rep(0, length(grid))

for(i in 1:length(grid)){
  y[i] = cfun(grid[i])
}

# c.df = data.frame( y, grid)
# # The ture coefficient functions
# c1.df = data.frame(beta12 = c.df$y, beta13 = rep(0, length(c.df$grid)), grid = c.df$grid)
# mc1.df = melt(c1.df, id.vars = "grid")
# c2.df = data.frame(beta21 = c.df$y, beta23 = rep(0.6, length(c.df$grid)), grid = c.df$grid)
# mc2.df = melt(c2.df, id.vars = "grid")
# c3.df = data.frame(beta31 = rep(0, length(c.df$grid)), beta32 = rep(0.6, length(c.df$grid)), grid = c.df$grid)
# mc3.df = melt(c3.df, id.vars = "grid")
# 
# mc.df = rbind(mc1.df, mc2.df, mc3.df)
# 
# ggplot(data = mc.df, aes(x = grid, y = -value)) +
#   facet_wrap(~variable, ncol = 2) +
#   geom_line(col = "blue") +
#   scale_y_continuous(limits = c(-1,1)) +
#   labs(x = "time", y = "partical correlation", title = "Partial correlation over time for each node")
# 
# save(c.df, mc.df, file = "R/3d_Case/3dexample.RData")

Omegafun = function(t, cfun){
  M = matrix(c(1,       cfun(t),   0,   
               cfun(t), 1,         0.6,
               0,       0.6,         1), 
             ncol = 3, nrow = 3)
  colnames(M) = paste0("X", c(1:3))
  rownames(M) = paste0("X", c(1:3))
  return(M)
}




