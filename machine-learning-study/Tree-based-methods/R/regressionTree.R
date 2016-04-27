#################################################################
## regressionTree.R
## Apply regression trees on Hitters dataset.
## Author: Meilei
## Date: April 2016
#################################################################
library(tree)
library(ggplot2)
library(ggdendro)
library(dplyr)
library(ISLR)

# use classification trees to analyze Carseats data set.

Carseats = Carseats

High = ifelse(Sales <= 8, "No", "Yes")

Carseats = data.frame(Carseats, High)

tree.carseats = tree(High~.-Sales, Carseats)
summary(tree.carseats)

tree_data = dendro_data(tree.carseats)

ggplot(segment(tree_data)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend, size=n), 
               colour="lightblue") +
  scale_size("n") +
  geom_text(data=label(tree_data), 
            aes(x=x, y=y, label=label), vjust=-0.5, size=4) +
  geom_text(data=leaf_label(tree_data), 
            aes(x=x, y=y, label=label), vjust=0.5, size=3) +
  theme_dendro()

plot(tree.carseats)
text(tree.carseats, pretty = 0)
