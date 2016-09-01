#########################################################################################
## Epan_kernel.R
## Epanechnikov kernel.
## Author: Meilei
## Date: April 2016
##########################################################################################
library(ggplot2)
library(dplyr)

epan = function(z, h = 1){
  if(abs(z/h) <= 1){
    return(3/4*(1-z^2/h^2))
  } else{
    return(0)
  }
}
