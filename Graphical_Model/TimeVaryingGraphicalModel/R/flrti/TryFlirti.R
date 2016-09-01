###########################################################################################
## TryFlirti.R
## Try the code of FLiRTI
## Author: Meilei
## Date: April 2016
############################################################################################
library(dplyr)
library(ggplot2)
library(reshape2)
library(lpSolve)

source("R/flrti/flrti.R")
ls()
# [1] "flrti"         "flrti.boot"    "flrti.cv"      "flrti.perm"    "linearoptim"   "makeA"         "predict.flrti"
# [8] "testdata" 
attach(testdata)
testfit <- flrti(Y,X,sigma=.001,weight=.1,plot=T, deriv=3, wrap = T)
lines(seq(0,1,length=100),beta,col=2)
