############################################################################
## vanillaSimulation.R
## This is a simulation study on vanilla neural net for regression.
##
## Author: Meilei
############################################################################
library(dplyr)
library(ggplot2)
source("neural-network/R/sigmoid.R")

# generate simulation data ------------------------------------------------
# 2 random input varibles and 1 bias varibles, 100 training samples
train = c(rep(1, 100), rnorm(200, 0, 1)) 
trainX = matrix(train, ncol = 3, byrow = F)
trainY = matrix(rep(0, 100), ncol = 1)
a1 = matrix(c(0, 3, 3), ncol = 3)
a2 = matrix(c(0, 3,-3), ncol = 3)

for(i in 1:100){
  trainY[i] = sigma(a1 %*% trainX[i,]) + sigma(a2 %*% trainX[i,]) + rnorm(1,0,0.01)
}


# set up hidden variables -------------------------------------------------

## initialize the alpha (parameter of hidden variable functions): set them close to zeros
# 5 hidden variables and 1 bias variables
alpha = matrix(rnorm(15, 0, 0.01), ncol = 3, byrow = T)

## initialize the beta (parameter of output variable functions): set them close to zeros
# 1 output varibles: beta[6] is the bias
beta = matrix(rnorm(6, 0, 0.01), ncol = 6, byrow = T)


# back-propagation --------------------------------------------------------

for(n in 1:100000){
  r = 1/n^0.3 # learning rate
  
  ## current value of hidden variables and bias variable(always 1)
  tempZ = matrix(c(rep(0, 500), rep(1,100)), ncol = 6)
  for(i in 1:100){
    for(m in 1:5){
      tempZ[i,m] = sigma(alpha[m,] %*% trainX[i,] )    
    }
  }
  
  ## current output
  tempY =  tempZ %*% t(beta)
  temdY = sum((tempY - trainY)^2)
  
  ## current derivatives of beta w.r.t. squared error loss and update the beta
  deltaB  = matrix(rep(0, 600), ncol = 6)
  for(i in 1:100){
    for(m in 1:6){
      deltaB[i,m] = -2*(trainY[i,1] - tempY[i,1])*tempZ[i,m]       
    }
  }
  newbeta = beta - r * colMeans(deltaB)
  
  dB = sum((beta - newbeta)^2)
  ## current derivatives of alpha w.r.t. squared error loss and update the alpha
  newalpha = matrix(rep(0, 15), ncol = 3, byrow = T)
  for(m in 1:5){
    deltaA = matrix(rep(0, 300), ncol = 3) 
    for(i in 1:100){
      for(l in 1:3){
        deltaA[i,l] = deltaB[i,m] * beta[m] * (1 - tempZ[i,m]) * trainX[i,l]       
      }
    }
    newalpha[m,] = alpha[m,] - r * colMeans(deltaA)
  }
  dA = sum((alpha - newalpha)^2)
  if(dB  < 0.00000001 & dA < 0.00000001){
    print(n)
    break 
                                   }else{
    alpha = newalpha
    beta = newbeta
  } 
}



## update the value of hidden variables 

newZ = matrix(c(rep(0, 500), rep(1,100)), ncol = 6)
for(i in 1:100){
  for(m in 1:5){
    newZ[i,m] = sigma(newalpha[m,] %*% trainX[i,] )    
  }
}

## current output
newY = newZ %*% t(newbeta)
newdY = sum((newY - trainY)^2)/100
plot(trainY, newY)
abline(0,1)
