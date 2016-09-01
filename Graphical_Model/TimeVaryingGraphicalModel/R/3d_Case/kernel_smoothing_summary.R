################################################################################
## kernel_smoothing_summary.R
## Author: Meilei
## Looking at a range of bandwidth values in the kernel smoothing methods on 3d
## data.
## Date: July 16th.
################################################################################
library(ggplot2)
library(dplyr)
library(glmnet)
library(gridExtra)
library(reshape2)
source("R/3d_Case/kernel_smoothing_result.R")


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


# visualizing the kernel smoothing result ---------------------------------


gh1 = kernel_smoothing_result(Data = Data, time = time, h = 1/40, mc.df = mc.df)
gh2 = kernel_smoothing_result(Data = Data, time = time, h = 1/30, mc.df = mc.df)
gh3 = kernel_smoothing_result(Data = Data, time = time, h = 1/20, mc.df = mc.df)
gh4 = kernel_smoothing_result(Data = Data, time = time, h = 1/15, mc.df = mc.df)
gh5 = kernel_smoothing_result(Data = Data, time = time, h = 1/10, mc.df = mc.df)
gh6 = kernel_smoothing_result(Data = Data, time = time, h = 1/6, mc.df = mc.df)



pdf(file = "figures/kernel_smooth_glmnet.pdf", height = 12, width = 8)
gh1
gh2
gh3
gh4
gh5
gh6
dev.off()

