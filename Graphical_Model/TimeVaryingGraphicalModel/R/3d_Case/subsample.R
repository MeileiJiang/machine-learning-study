#########################################################
## subsample.R
## Generate subsample for time varying graphical model
## dataset.
## Author: Meilei
#########################################################
library(dplyr)

# load("rdata/3dtoyexample1.Rdata")
# table(Toydata$Time)

subsample = function(Data, frac = 0.5){
  time = unique(Data$Time)
  newdata = NULL
  for(t in time){
    datat = Data %>%
      filter(Time == t) 
    nt = dim(datat)[1]
    datat = datat %>%
      sample_n(ceiling(nt*frac))
    newdata = rbind(newdata, datat)
  }
  return(newdata)
}

# subData = subsample(Data = Toydata, frac = 0.1)


