###############################################################################################
## basis_expansion_result.R
## author: Meilei Jiang
## Date: July 19th.
###############################################################################################
library(ggplot2)
library(dplyr)
library(reshape2)
source("R/3d_Case/basis_expansion_tv.R")

# loading the data --------------------------------------------------------


# load the true coefficent function
load(file = "R/3d_Case/3dexample.RData") 

# ggplot(data = mc.df, aes(x = grid, y = -value)) +
#   facet_wrap(~variable, ncol = 2) +
#   geom_line(col = "blue") +
#   scale_y_continuous(limits = c(-1,1)) +
#   labs(x = "time", y = "partical correlation", title = "Partial correlation over time for each node")
# load the simulation data
load("rdata/3dtoyexample1.Rdata")

# data processing ---------------------------------------------------------

Data = Toydata[,-4]
time = unique(Toydata$Time)
n = length(Toydata$Time)/length(time)

delta = 0.05
grid = seq(0, 1, by = delta)

# node 1
X1 = as.matrix(Data[,-1])
Y1 = as.matrix(Data[,1])

# node 2
X2 = as.matrix(Data[,-2])
Y2 = as.matrix(Data[,2])

# node 3
X3 = as.matrix(Data[,-3])
Y3 = as.matrix(Data[,3])




## Set the parameter for plots

set = seq(0, 1, by = 0.001)


# basis expasion for each node -------------------------------------------

Beta1 = basis_expansion_tv(X = X1, Y = Y1, time = time, grid, order = 4, w = 0.0005, set = set)
Beta2 = basis_expansion_tv(X = X2, Y = Y2, time = time, grid, order = 4, w = 0.0005, set = set)
Beta3 = basis_expansion_tv(X = X3, Y = Y3, time = time, grid, order = 4, w = 0.0005, set = set)

beta1 = Beta1$beta
beta2 = Beta2$beta
beta3 = Beta3$beta


beta1df = data.frame(beta1, set)
names(beta1df) = c("beta12", "beta13", "time")
mbeta1df = melt(beta1df, id.vars = "time")

beta2df = data.frame(beta2, set)
names(beta2df) = c("beta21", "beta23", "time")
mbeta2df = melt(beta2df, id.vars = "time")

beta3df = data.frame(beta3, set)
names(beta3df) = c("beta31", "beta32", "time")
mbeta3df = melt(beta3df, id.vars = "time")

mbetadf = rbind(mbeta1df, mbeta2df, mbeta3df)

ggplot(data = mbetadf, aes(x = time, y = value),size = 1) +
  facet_wrap(~variable, ncol = 2) +
  geom_line() +
  #  geom_point(size = 1) +
#  geom_line(data = mc.df, aes(x = grid, y = -value),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-1, 1),
                     seq(-1, 1, by = 0.2),
                     seq(-1, 1, by = 0.2)) + 
  scale_shape(solid = FALSE) +
  labs(x = "time", y = "partial correlation", title = "") +
  ylab("partial correlation")
