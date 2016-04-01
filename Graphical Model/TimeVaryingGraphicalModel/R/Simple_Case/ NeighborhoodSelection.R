###################################################################################
## NeighborhoodSelection.R
## Author: Meilei
## Apply nearest neighborhood selection at each time point and then catenate them together
## Date: April 2016
###################################################################################
library(glmnet)
library(ggplot2)
library(igraph)
library(dplyr)
library(reshape2)

source("R/Simple_Case/2dexample.R")

load("rdata/2dtoyexample1.Rdata")

Toydata %>% 
  group_by(Time) %>%
  summarise(n = n())
table(Toydata$Time)

time = unique(Toydata$Time)

S = 0.18

# node 1 ------------------------------------------------------------------

beta12 = rep(0, length(time))

for(t in 1: length(time)){
  subdata = Toydata%>%filter(Time == time[t]) 
  subdata = subdata[,-11]
  X = as.matrix(subdata[,-1])
  Y = as.matrix(subdata[, 1])
  fit.l1 = glmnet(x= X, y = Y, intercept = F)
  beta12[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X2",]
}

beta1.df = data.frame(beta12, time)
mbeta1 = melt(beta1.df, id.vars = "time")

g1 = ggplot(mbeta1, aes(x = time, y = value)) +
  facet_wrap(~variable, ncol = 3) +
  geom_point(col= "blue") +
  geom_line(col ="red") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0),
                     labels = c(-0.5,-0.4,-0.3,-0.2,-0.1,0)) +
  labs(x = "time", title = "Coefficients of Node 1")
g1

