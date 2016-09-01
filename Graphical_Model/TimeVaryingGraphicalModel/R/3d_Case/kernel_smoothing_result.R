################################################################################
## kernel_smoothing_result.R
## Author: Meilei
## Applying kernel_smoothing on simulation data.
## Date: July 16th.
################################################################################
library(ggplot2)
library(dplyr)
library(glmnet)
library(gridExtra)
library(reshape2)

source("R/3d_Case/Kernel_Smoothing_glmnet.R")



# kernel smoothing result function ----------------------------------------


kernel_smoothing_result = function(Data, time, hh, mc.df){
  # node 1
  X1 = as.matrix(Data[,-1])
  Y1 = as.matrix(Data[,1])
  
  # node 2
  X2 = as.matrix(Data[,-2])
  Y2 = as.matrix(Data[,2])
  
  # node 3
  X3 = as.matrix(Data[,-3])
  Y3 = as.matrix(Data[,3])
  
  # estimating the coefficient function for each node
  beta1 = kern_smooth(h=hh, X = X1, Y = Y1)
  beta2 = kern_smooth(h=hh, X = X2, Y = Y2)
  beta3 = kern_smooth(h=hh, X = X3, Y = Y3)
  
  beta1df = data.frame(beta1, time)
  names(beta1df) = c("beta12", "beta13", "time")
  mbeta1df = melt(beta1df, id.vars = "time")
  
  beta2df = data.frame(beta2, time)
  names(beta2df) = c("beta21", "beta23", "time")
  mbeta2df = melt(beta2df, id.vars = "time")
  
  beta3df = data.frame(beta3, time)
  names(beta3df) = c("beta31", "beta32", "time")
  mbeta3df = melt(beta3df, id.vars = "time")
  
  mbetadf = rbind(mbeta1df, mbeta2df, mbeta3df)
  
  g = ggplot(data = mbetadf, aes(x = time, y = value),size = 1) +
    facet_wrap(~variable, ncol = 2) +
    geom_line() +
 #   geom_point(size = 1) +
    geom_line(data = mc.df, aes(x = grid, y = -value),linetype ="dashed", col = "blue") +
    scale_x_continuous(limits = c(0, 1), 
                       breaks = seq(0, 1, by = 0.2), 
                       labels =  seq(0, 1, by = 0.2)) +
    scale_y_continuous(limits = c(-1, 1),
                       seq(-1, 1, by = 0.2),
                       seq(-1, 1, by = 0.2)) + 
    scale_shape(solid = FALSE) +
    labs(x = "time", y ="value", title = paste0("Coefficient Estimation of Each Node \n
         Kernel smoothing, h = ", sprintf("%.3f", round(hh,3))))
  return(g)
}




