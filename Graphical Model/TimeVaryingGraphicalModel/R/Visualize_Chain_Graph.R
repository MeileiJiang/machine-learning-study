#######################################################################################################################
## Visualize_Chain_Graph.R
## 
## Author: Meilei
## Date: March 2016
#######################################################################################################################
library(igraph)
library(ggplot2)
library(dplyr)
library(Matrix)



b12 = function(t){
  if(t >= 0 & t <= 0.2) return(0.5 - 5*abs(t-0.1))
  if(t > 0.8 & t <= 1) return(0.5 - 5*abs(t-0.9))
  return(0)
}

b23 = function(t){
  if(t >= 0 & t <= 0.2) return(0.5 - 5*abs(t-0.1))
  if(t > 0.2 & t <= 0.4) return(0.5 - 5*abs(t-0.3))
  return(0)
}

b34 = function(t){
  if(t >= 0.2 & t <= 0.4) return(0.5 - 5*abs(t-0.3))
  if(t > 0.4 & t <= 0.6) return(0.5 - 5*abs(t-0.5))
  return(0)
}

b45 = function(t){
  if(t >= 0.4 & t <= 0.6) return(0.5 - 5*abs(t-0.5))
  if(t > 0.6 & t <= 0.8) return(0.5 - 5*abs(t-0.7))
  return(0)
}

b51 = function(t){
  if(t >= 0.6 & t <= 0.8) return(0.5 - 5*abs(t-0.7))
  if(t>0.8 & t<=1) return(0.5 - 5*abs(t-0.9))
  return(0)
}







OmegaFun = function(t){
  ## Dynamic Part
  Omega1 = matrix(c(1,      b12(t), 0,      0,      b51(t),  
                    b12(t), 1,      b23(t), 0,      0, 
                    0,      b23(t), 1,      b34(t), 0,  
                    0,      0,      b34(t), 1,      b45(t), 
                    b51(t), 0,      0,      b45(t), 1), 
                  ncol = 5, nrow = 5)
  ## Static Part
  Omega2 = matrix(c(1, .5, 0, 0, 0.5,  
                    .5, 1, .5, 0, 0, 
                    0, .5, 1, .5, 0,  
                    0, 0, .5, 1, .5, 
                    0.5, 0, 0, .5, 1), ncol = 5, nrow = 5)
  ## Joint Graph
  M = as.matrix(bdiag(Omega1, Omega2))
  colnames(M) = paste0("X", c(1:10))
  rownames(M) = paste0("X", c(1:10))
  return(M)
}


# Omega = OmegaFun(0.9)
# 
# O = (Omega - diag(rep(1, 10)) != 0)
# 
# G = graph_from_adjacency_matrix(O, mode = "undirected")
# l <-layout.circle(G)
# plot.igraph(G,  vertex.size=10, vertex.label.dist = 0, vertex.color="red", layout = l)





