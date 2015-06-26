#######################################################
## TreeClassification.R
## Build up a two terminal-node classification tree.
## Author: Meilei Jiang
#######################################################
library(dplyr)
library(tree)
# build up a simulation data set. The features X1, ..., X10 are standard
# independent Gaussians, and the deterministic target Y is define by Y = 1
# if sum Xi^2 > 9.34; Y = -1 otherwise.
train = rnorm(20000, 0, 1) 
trainX = data.frame(matrix(train, ncol = 10, byrow = T))
colnames(trainX) = paste("X",c(1:10),sep ="")
test = rnorm(100000,0,1)
testX = data.frame(matrix(test, ncol = 10, byrow = T))
colnames(testX) = paste("X",c(1:10),sep = "")


classY= function(x){
  n = length(x)
  y = rep(0,n)
  for(i in 1:n){
    if(x[i] > qchisq(0.5,10)) {
      y[i] = 1 }else y[i] = -1
  } 
  y
}

trainY = classY(rowSums(trainX^2))
trainY = data.frame(Y = trainY)
table(trainY)
# trainY
# -1    1 
# 1018  982 
testY = classY(rowSums(testX^2))
testY = data.frame(Y = testY)
table(testY)

Train = cbind(trainY, trainX)
mytree = prune.tree(tree(as.factor(Y) ~ ., data=Train), best = 2)
text(mytree)
summary(mytree)
mytree.fit = predict(mytree, newdata= trainX)
fitY = as.numeric(mytree.fit[,1] < 0.5 )
fitY[fitY==0] = -1
sum((trainY- fitY)^2/4)

mytree.test = predict(mytree, newdata = testX)
predictY = as.numeric(mytree.test[,1] < 0.5 )
predictY[predictY==0] = -1
sum((testY- predictY)^2/4)/10000
# [1] 0.4614

