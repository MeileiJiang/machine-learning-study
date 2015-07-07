#########################################################################################################################
## TransformData.R
## Transform the daily category behavior to a 0-1 matrix.
##
## Author: Meilei
#########################################################################################################################
library(dplyr)
load("neural-network/Data/categoryEachDay.RData")

daydf = day_time_category_df
behavior = unique(unlist(daydf))
