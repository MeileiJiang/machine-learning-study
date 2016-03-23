#########################################
## lasso.R
## Understand the drawback of lasso.
## Author: Meilei
#########################################
library(glmnet)
library(dplyr)
library(ggplot2)



# set the data ------------------------------------------------------------

n = 1000;
p = n + 1;
eps = 0.0001;

X = cbind(diag(rep(1, n)), c(rep(eps, n-1), 1))
beta = c(rep(0, n-1), -1, 1)/eps
Y = X %*% beta


# least absolute shrinkage and selection operator -------------------------------------------------

fit = glmnet(X, Y)
cvfit = cv.glmnet(X, Y)
plot(fit)
plot(cvfit)
print(fit)

# 

# high dimension case -----------------------------------------------------

x=cbind(rep(0, 10), matrix(rnorm(10*10),10,10), rep(0, 10), 
        matrix(rnorm(10*10),10,10), rep(0, 10), 
        matrix(rnorm(10*10),10,10), rep(0, 10),
        matrix(rnorm(10*10),10,10), rep(0, 10), 
        matrix(rnorm(10*10),10,10), rep(0, 10))
y=rnorm(10)

fit1=glmnet(x,y)
plot(fit1)
