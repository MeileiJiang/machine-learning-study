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
# 2 random input varibles and 1 bias varibles
a1 = matrix(c(0, 3, 3), ncol = 3)
a2 = matrix(c(0, 3,-3), ncol = 3)

# 100 training samples
train = c(rep(1, 100), rnorm(200, 0, 1)) 
trainX = matrix(train, ncol = 3, byrow = F)
trainY = matrix(rep(0, 100), ncol = 1)
for(i in 1:100){
  trainY[i] = sigma(a1 %*% trainX[i,]) + sigma(a2 %*% trainX[i,]) + rnorm(1,0,1)
}
ggplot(data.frame(trainY), aes(x = trainY)) + geom_density()


test = c(rep(1, 10000), rnorm(20000, 0 ,1))
testX = matrix(test, ncol = 3, byrow =F)
testY = matrix(rep(0, 10000), ncol = 1)
for(i in 1:10000){
  testY[i] = sigma(a1 %*% testX[i,]) + sigma(a2 %*% testX[i,]) + rnorm(1,0,1)
}

ggplot(data.frame(testY), aes(x = testY)) + geom_density() +
  geom_density(data = data.frame(trainY), aes(x = trainY), col = "red")


# linear regression -------------------------------------------------------
lmtrain = data.frame(trainY, trainX)
lmtest = data.frame(testY, testX)
lm0 = lm(trainY~ X1 + X2 + X3 + 0, data = lmtrain)
summary(lm0)
lm.fit = predict(lm0)
mean((lm.fit - trainY)^2)
#[1] 0.05575285

lm.preY=predict(lm0, newdata = lmtest)
mean((lm.preY - testY)^2)
# [1] 0.05747757

# set up hidden variables -------------------------------------------------

## initialize the alpha (parameter of hidden variable functions): set them close to zeros
# 5 hidden variables and 1 bias variables
alpha = matrix(rnorm(15, 0, 0.01), ncol = 3, byrow = T)

## initialize the beta (parameter of output variable functions): set them close to zeros
# 1 output varibles: beta[6] is the bias
beta = matrix(rnorm(6, 0, 0.01), ncol = 6, byrow = T)


# back-propagation --------------------------------------------------------

for(n in 1:100000){
  r = 1/n^0.6 # learning rate
  
  ## current value of hidden variables and bias variable(always 1)
  tempZ = matrix(c(rep(0, 500), rep(1,100)), ncol = 6)
  for(i in 1:100){
    for(m in 1:5){
      tempZ[i,m] = sigma(alpha[m,] %*% trainX[i,] )    
    }
  }
  
  ## current output
  tempY =  tempZ %*% t(beta)
  
  ## current derivatives of beta w.r.t. squared error loss and update the beta
  deltaB  = matrix(rep(0, 600), ncol = 6)
  for(i in 1:100){
    for(m in 1:6){
      deltaB[i,m] = -2*(trainY[i,1] - tempY[i,1])*tempZ[i,m]       
    }
  }
  newbeta = beta - r * (colMeans(deltaB) - 0.04 * beta/(1 + beta^2)^2) # adding 2 beta to avoid overfitting
  
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
    newalpha[m,] = alpha[m,] - r * (colMeans(deltaA) - 0.04 * alpha[m,]/(1 + alpha[m,]^2)^2)# adding 2 alpha to avoid overfitting
  }
  dA = sum((alpha - newalpha)^2)
  if(dB  < 0.0000000001 & dA < 0.00000001){
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
# training error of NN
mean((newY - trainY)^2)
# [1] 0.001404938 
plot(trainY, newY)
abline(0,1)


# test the model ----------------------------------------------------------

testZ = matrix(c(rep(0, 50000), rep(1,10000)), ncol = 6)
for(i in 1:100){
  for(m in 1:5){
    testZ[i,m] = sigma(newalpha[m,] %*% testX[i,] )    
  }
}
preY = testZ %*% t(newbeta)
# test error 
mean((preY - testY)^2)
# [1] 0.9189864
plot(testY, preY)
abline(0,1)

ggplot(data.frame(testY), aes(x = testY)) + geom_density() +
  geom_density(data = data.frame(preY), aes(x = preY), col = "red")



