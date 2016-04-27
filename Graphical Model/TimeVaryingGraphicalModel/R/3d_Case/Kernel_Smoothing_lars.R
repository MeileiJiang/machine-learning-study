##########################################################################
## Kernel_Smoothing_lars.R
## Author: Meilei
## Apply kernel smoothing on time varying 3d gaussian data by lars regression.
## Date: April 2016
#########################################################################
library(ggplot2)
library(gridExtra)
library(dplyr)
library(lars)
library(reshape2)

# load the kernel function
source(file = "R/3d_Case/Epan_kernel.R")

# load the true coefficent function
load(file = "R/3d_Case/3dexample.RData")

# load the simulation data
load("rdata/3dtoyexample1.Rdata")


# data processing ---------------------------------------------------------

Data = Toydata[,-4]
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
  beta13 = NULL
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
    beta12[k] = coef(wfit)[2,1]
    beta13[k] = coef(wfit)[2,2]
  }
  beta1 = data.frame(beta12, beta13)
  return(beta1)
}



# visualize the result ----------------------------------------------------

# Case 1: h = 1/40 
beta1 = kern_smooth(h=1/40)

beta12 = data.frame(beta12 = beta1$beta12, time)
mbeta12 = melt(beta12, id.vars = "time")
beta13 = data.frame(beta12 = beta1$beta13, time)
mbeta13 = melt(beta13, id.vars = "time")

g1.1 = ggplot(data = mbeta12, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 2\n
       Kernel smoothing, h = 1/40")
g1.1

g2.1 = ggplot(data = mbeta13, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = 0,linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 3\n
       Kernel smoothing, h = 1/40")
g2.1

# Case 2: h = 1/30 
beta1 = kern_smooth(h=1/30)

beta12 = data.frame(beta12 = beta1$beta12, time)
mbeta12 = melt(beta12, id.vars = "time")
beta13 = data.frame(beta12 = beta1$beta13, time)
mbeta13 = melt(beta13, id.vars = "time")

g1.2 = ggplot(data = mbeta12, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 2\n
       Kernel smoothing, h = 1/30")
g1.2

g2.2 = ggplot(data = mbeta13, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = 0,linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 3\n
       Kernel smoothing, h = 1/30")
g2.2


# Case 3: h = 1/20 
beta1 = kern_smooth(h=1/20)

beta12 = data.frame(beta12 = beta1$beta12, time)
mbeta12 = melt(beta12, id.vars = "time")
beta13 = data.frame(beta12 = beta1$beta13, time)
mbeta13 = melt(beta13, id.vars = "time")

g1.3 = ggplot(data = mbeta12, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 2\n
       Kernel smoothing, h = 1/20")
g1.3

g2.3 = ggplot(data = mbeta13, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = 0,linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 3\n
       Kernel smoothing, h = 1/20")
g2.3


# Case 4: h = 1/15 
beta1 = kern_smooth(h=1/15)

beta12 = data.frame(beta12 = beta1$beta12, time)
mbeta12 = melt(beta12, id.vars = "time")
beta13 = data.frame(beta12 = beta1$beta13, time)
mbeta13 = melt(beta13, id.vars = "time")

g1.4 = ggplot(data = mbeta12, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 2\n
       Kernel smoothing, h = 1/15")
g1.4

g2.4 = ggplot(data = mbeta13, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = 0,linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 3\n
       Kernel smoothing, h = 1/15")
g2.4


# Case 5: h = 1/10 
beta1 = kern_smooth(h=1/10)

beta12 = data.frame(beta12 = beta1$beta12, time)
mbeta12 = melt(beta12, id.vars = "time")
beta13 = data.frame(beta12 = beta1$beta13, time)
mbeta13 = melt(beta13, id.vars = "time")

g1.5 = ggplot(data = mbeta12, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 2\n
       Kernel smoothing, h = 1/10")
g1.5

g2.5 = ggplot(data = mbeta13, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = 0,linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 3\n
       Kernel smoothing, h = 1/10")
g2.5


# Case 6: h = 1/6
beta1 = kern_smooth(h=1/6)

beta12 = data.frame(beta12 = beta1$beta12, time)
mbeta12 = melt(beta12, id.vars = "time")
beta13 = data.frame(beta12 = beta1$beta13, time)
mbeta13 = melt(beta13, id.vars = "time")

g1.6 = ggplot(data = mbeta12, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_line(data = c.df, aes(x = grid, y = -y),linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 2\n
       Kernel smoothing, h = 1/6")
g1.6

g2.6 = ggplot(data = mbeta13, aes(x = time, y = value)) +
  facet_wrap(~variable) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = 0,linetype ="dashed", col = "blue") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, by = 0.2), 
                     labels =  seq(0, 1, by = 0.2)) +
  scale_y_continuous(limits = c(-0.65, 0.65),
                     breaks = c(-5:5)/10,
                     labels = c(-5:5)/10) +
  labs(x = "time", title = "Coefficient Estimation of Node 1 over Node 3\n
       Kernel smoothing, h = 1/6")
g2.6



pdf(file = "figures/kernel_smooth_lars.pdf", height = 12, width = 10)
grid.arrange(g1.1, g1.2, g1.3, g1.4, g1.5, g1.6, ncol = 2 )
dev.off()

pdf(file = "figures/kernel_smooth_lars2.pdf", height = 12, width = 10)
grid.arrange(g2.1, g2.2, g2.3, g2.4, g2.5, g2.6, ncol = 2 )
dev.off()
