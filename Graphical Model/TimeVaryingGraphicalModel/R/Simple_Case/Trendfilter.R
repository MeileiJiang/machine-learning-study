#####################################################################################
## Trendfilter.R
## Try to use trendfilter to solve the problem.
## Author: Meilei
## Date: April 2015
####################################################################################
library(ggplot2)
library(dplyr)
library(genlasso)
library(reshape2)

load("rdata/2dtoyexample4.Rdata")

Data = Toydata[,-3]
time = unique(Toydata$Time)
n = length(Toydata$Time)/length(time)

X1 = as.matrix(Data[,-1] )
Y1 = as.matrix(Data[,1])

a = trendfilter(y = Y1,  ord = 3)
plot(a, nlam = 2)
