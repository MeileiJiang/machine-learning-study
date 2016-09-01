##########################################################################
## Kernel_Smoothing.R
## Author: Meilei
## Apply kernel smoothing on time varying data set.
## Date: April 2016
#########################################################################
library(ggplot2)
library(gridExtra)
library(dplyr)
library(glmnet)
library(lars)

# load the kernel function
source(file = "R/Simple_Case/Epan_kernel.R")

# load the true coefficent function
load(file = "R/Simple_Case/2dexample.RData")

# load the simulation data
load("rdata/2dtoyexample3.Rdata")


# data processing ---------------------------------------------------------

Data = Toydata[,-3]
time = unique(Toydata$Time)
n = length(Toydata$Time)/length(time)

grid = seq(0, 1, by = 0.05)

# node 1
X1 = as.matrix(Data[,-1] )
# node 2
Y1 = as.matrix(Data[,1])


# estimate the coefficient on each point of grid --------------------------

h = 1/40
kern_smooth = function(h, X =X1, Y = Y1, Tick = time, Time = Toydata$Time){
  beta12 = NULL
  for(k in c(1:101)){
    tick = Tick[k]
    weights <- NULL
    for(t in c(1:length(Time))){
      weights[t] = epan(z = (tick - Time[t]), h)
    }
    weights = weights/sum(weights)
    wX = X * sqrt(weights)
    wY = Y * sqrt(weights)
    wfit = lars(x= wX, y = wY, normalize = F, intercept = F)
    beta12[k] = coef(wfit)[2]
  }
  return(beta12)
}


# plot(beta12, type = "l")


# visualize the result ----------------------------------------------------

# Case 1: h = 1/40 
beta12 = kern_smooth(h=1/40)

beta1 = data.frame(beta12, time)
mbeta1 = melt(beta1, id.vars = "time")

g1.1 = ggplot(data = mbeta1, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point() +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Time Varying Coefficient Estimation of Node 1 \n
       Kernel smoothing, h = 1/40")
g1.1

# Case 2: h = 1/30 
beta12 = kern_smooth(h=1/30)

beta1 = data.frame(beta12, time)
mbeta1 = melt(beta1, id.vars = "time")

g1.2 = ggplot(data = mbeta1, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point() +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Time Varying Coefficient Estimation of Node 1 \n
       Kernel smoothing, h = 1/30")
g1.2

# Case 3: h = 1/20 
beta12 = kern_smooth(h=1/20)

beta1 = data.frame(beta12, time)
mbeta1 = melt(beta1, id.vars = "time")

g1.3 = ggplot(data = mbeta1, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point() +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Time Varying Coefficient Estimation of Node 1 \n
       Kernel smoothing, h = 1/20")
g1.3

# Case 4: h = 1/15 
beta12 = kern_smooth(h=1/15)

beta1 = data.frame(beta12, time)
mbeta1 = melt(beta1, id.vars = "time")

g1.4 = ggplot(data = mbeta1, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point() +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Time Varying Coefficient Estimation of Node 1 \n
       Kernel smoothing, h = 1/15")
g1.4

# Case 5: h = 1/10 
beta12 = kern_smooth(h=1/10)

beta1 = data.frame(beta12, time)
mbeta1 = melt(beta1, id.vars = "time")

g1.5 = ggplot(data = mbeta1, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point() +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Time Varying Coefficient Estimation of Node 1 \n
       Kernel smoothing, h = 1/10")
g1.5

# Case 6: h = 1/6
beta12 = kern_smooth(h=1/6)

beta1 = data.frame(beta12, time)
mbeta1 = melt(beta1, id.vars = "time")

g1.6 = ggplot(data = mbeta1, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point() +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Time Varying Coefficient Estimation of Node 1 \n
       Kernel smoothing, h = 1/6")
g1.6

pdf(file = "figures/kernel_smooth.pdf", height = 12, width = 10)
grid.arrange(g1.1, g1.2, g1.3, g1.4, g1.5, g1.6, ncol = 2 )
dev.off()
