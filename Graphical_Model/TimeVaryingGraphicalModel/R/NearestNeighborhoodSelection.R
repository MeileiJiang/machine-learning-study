###################################################################################
## NearestNeighborhoodSelection.R
## Author: Meilei
## Apply nearest neighborhood selection at each time point and then catenate them together
## Date: March 2016
###################################################################################
library(glmnet)
library(ggplot2)
library(igraph)
library(dplyr)
library(reshape2)

source("R/Visualize_Chain_Graph.R")
source("R/SimluationData.R")

load("rdata/toyexample4.Rdata")

Toydata %>% 
  group_by(Time) %>%
  summarise(n = n())
table(Toydata$Time)

time = unique(Toydata$Time)

S = 0.18

# node 1 ------------------------------------------------------------------

beta12 = rep(0, length(time))
beta13 = rep(0, length(time))
beta14 = rep(0, length(time))
beta15 = rep(0, length(time))
beta16 = rep(0, length(time))
beta17 = rep(0, length(time))
beta18 = rep(0, length(time))
beta19 = rep(0, length(time))
beta110 = rep(0, length(time))

for(t in 1: length(time)){
  subdata = Toydata%>%filter(Time == time[t]) 
  subdata = subdata[,-11]
  X = as.matrix(subdata[,-1])
  Y = as.matrix(subdata[, 1])
  fit.l1 = glmnet(x= X, y = Y, intercept = F)
  beta12[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X2",]
  beta13[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X3",]
  beta14[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X4",]
  beta15[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X5",]
  beta16[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X6",]
  beta17[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X7",]
  beta18[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X8",]
  beta19[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X9",]
  beta110[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X10",]
}

beta1.df = data.frame(beta12, beta13, beta14, beta15, beta16, 
                     beta17, beta18, beta19, beta110, time)
mbeta1 = melt(beta1.df, id.vars = "time")

g1 = ggplot(mbeta1, aes(x = time, y = value)) +
  facet_wrap(~variable, ncol = 3) +
  geom_point(col= "blue") +
  geom_line(col ="red") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.15),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0),
                     labels = c(-0.5,-0.4,-0.3,-0.2,-0.1,0)) +
  labs(x = "time", title = "Coefficients of Node 1")


# node 2 ------------------------------------------------------------------

beta21 = rep(0, length(time))
beta23 = rep(0, length(time))
beta24 = rep(0, length(time))
beta25 = rep(0, length(time))
beta26 = rep(0, length(time))
beta27 = rep(0, length(time))
beta28 = rep(0, length(time))
beta29 = rep(0, length(time))
beta210 = rep(0, length(time))

for(t in 1: length(time)){
  subdata = Toydata%>%filter(Time == time[t]) 
  subdata = subdata[,-11]
  X = as.matrix(subdata[,-2])
  Y = as.matrix(subdata[, 2])
  fit.l1 = glmnet(x= X, y = Y, intercept = F)
  beta21[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X1",]
  beta23[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X3",]
  beta24[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X4",]
  beta25[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X5",]
  beta26[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X6",]
  beta27[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X7",]
  beta28[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X8",]
  beta29[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X9",]
  beta210[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X10",]
}

beta2.df = data.frame(beta21, beta23, beta24, beta25, beta26, 
                     beta27, beta28, beta29, beta210, time)
mbeta2 = melt(beta2.df, id.vars = "time")

g2 = ggplot(mbeta2, aes(x = time, y = value)) +
  facet_wrap(~variable, ncol = 3) +
  geom_point(col= "blue") +
  geom_line(col ="red") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.15),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0),
                     labels = c(-0.5,-0.4,-0.3,-0.2,-0.1,0)) +
  labs(x = "time", title = "Coefficients of Node 2")

# node 3 ------------------------------------------------------------------

beta31 = rep(0, length(time))
beta32 = rep(0, length(time))
beta34 = rep(0, length(time))
beta35 = rep(0, length(time))
beta36 = rep(0, length(time))
beta37 = rep(0, length(time))
beta38 = rep(0, length(time))
beta39 = rep(0, length(time))
beta310 = rep(0, length(time))

for(t in 1: length(time)){
  subdata = Toydata%>%filter(Time == time[t]) 
  subdata = subdata[,-11]
  X = as.matrix(subdata[,-3])
  Y = as.matrix(subdata[, 3])
  fit.l1 = glmnet(x= X, y = Y, intercept = F)
  beta31[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X1",]
  beta32[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X2",]
  beta34[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X4",]
  beta35[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X5",]
  beta36[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X6",]
  beta37[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X7",]
  beta38[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X8",]
  beta39[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X9",]
  beta310[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X10",]
}

beta3.df = data.frame(beta31, beta32, beta34, beta35, beta36, 
                      beta37, beta38, beta39, beta310, time)
mbeta3 = melt(beta3.df, id.vars = "time")

g3 = ggplot(mbeta3, aes(x = time, y = value)) +
  facet_wrap(~variable, ncol = 3) +
  geom_point(col= "blue") +
  geom_line(col ="red") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.15),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0),
                     labels = c(-0.5,-0.4,-0.3,-0.2,-0.1,0)) +
  labs(x = "time", title = "Coefficients of Node 3")


# node 6 -----------------------------------------------------------------

beta61 = rep(0, length(time))
beta62 = rep(0, length(time))
beta63 = rep(0, length(time))
beta64 = rep(0, length(time))
beta65 = rep(0, length(time))
beta67 = rep(0, length(time))
beta68 = rep(0, length(time))
beta69 = rep(0, length(time))
beta610 = rep(0, length(time))

for(t in 1: length(time)){
  subdata = Toydata%>%filter(Time == time[t]) 
  subdata = subdata[,-11]
  X = as.matrix(subdata[,-6])
  Y = as.matrix(subdata[, 6])
  fit.l1 = glmnet(x= X, y = Y, intercept = F)
  beta61[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X1",]
  beta62[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X2",]
  beta63[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X3",]
  beta64[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X4",]
  beta65[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X5",]
  beta67[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X7",]
  beta68[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X8",]
  beta69[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X9",]
  beta610[t] = as.matrix(coef.glmnet(fit.l1, s = S ))["X10",]
}

beta6.df = data.frame(beta61, beta62, beta63, beta64, beta65, 
                      beta67, beta68, beta69, beta610, time)
mbeta6 = melt(beta6.df, id.vars = "time")

g6 = ggplot(mbeta6, aes(x = time, y = value)) +
  facet_wrap(~variable, ncol = 3) +
  geom_point(col= "blue") +
  geom_line(col ="red") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.15),
                     breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0),
                     labels = c(-0.5,-0.4,-0.3,-0.2,-0.1,0)) +
  labs(x = "time", title = "Coefficients of Node 6")


pdf(file = "figures/nearest_neighbor_select.pdf")
print(g1)
print(g2)
print(g3)
print(g6)
dev.off()

# t = matrix(0.1, ncol = 1, nrow = 1)
# Omega = OmegaFun(t)
# Sigma = solve(Omega)
# O = (Omega - diag(rep(1, 10)) != 0)
# 
# G = graph_from_adjacency_matrix(O, mode = "undirected")
# l <-layout.circle(G)
# plot.igraph(G,  vertex.size=10, vertex.label.dist = 0, vertex.color="red", layout = l)
# 
# Data = graphdata(t, n = 50)
# colnames(Data) = paste0("X",c(1:10))
# X = Data[,-2]
# Y = Data[,2]
# Data.df = data.frame(Data)
# cv.fit1 = cv.glmnet(x = X, y =Y, intercept = F)
# 
# fit1 = glmnet(x = X, y = Y,  intercept = F)
# plot(fit1)
# coef(fit1)
# 
# lm(X1~.+0, data = Data.df)
