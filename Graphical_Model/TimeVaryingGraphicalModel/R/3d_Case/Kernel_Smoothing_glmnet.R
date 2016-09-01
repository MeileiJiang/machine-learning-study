##########################################################################
## Kernel_Smoothing_glmnet.R
## Author: Meilei
## Apply kernel smoothing on time varying 3d gaussian data.
## Date: April 2016
#########################################################################
library(dplyr)
library(glmnet)

# load the kernel function
source(file = "R/3d_Case/Epan_kernel.R")



# estimate the coefficient on each point of grid --------------------------


kern_smooth = function(h, X =X1, Y = Y1, Tick = time, Time = Toydata$Time){
  beta12 = NULL
  beta13 = NULL
  varnames = colnames(X)
  for(k in c(1:101)){
    tick = Tick[k]
    weights <- NULL
    for(t in c(1:length(Time))){
      weights[t] = epan(z = (tick - Time[t]), h)
    }
    weights = weights/sum(weights)
    wcv = cv.glmnet(x =X, y = Y, weights = weights)
    lamb = wcv$lambda.1se
    wfit = glmnet(x =X, y = Y, weights = weights, lambda = lamb)
    beta12[k] = coef(wfit)[varnames[1],]
    beta13[k] = coef(wfit)[varnames[2],]
  }
    beta1 = data.frame(beta12, beta13)
  return(beta1)
}


