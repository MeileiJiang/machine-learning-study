###################################################################################
## Toyexample.R
## Author: Meilei
## Date: March 2016
###################################################################################
library(glmnet)
library(ggplot2)
library(igraph)
library(dplyr)

source("R/SimluationData.R")

time = seq(0, 1, by = 0.01)
Data = simulation(time = time, n = 100)

colnames(Data) = c(paste0("X",c(1:10)),"Time")
Toydata = data.frame(Data)
save(Toydata, file = "rdata/toyexample4.Rdata")

time = seq(0, 1, by = 0.01)
Data = simulation(time = time, n = 1)

colnames(Data) = c(paste0("X",c(1:10)),"Time")
Toydata = data.frame(Data)
save(Toydata, file = "rdata/toyexample5.Rdata")

time = seq(0, 1, by = 0.05)
Data = simulation(time = time, n = 20)

colnames(Data) = c(paste0("X",c(1:10)),"Time")
Toydata = data.frame(Data)
save(Toydata, file = "rdata/toyexample6.Rdata")
