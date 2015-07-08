#########################################################################################################################
## TransformData.R
## Transform the daily category behavior to a 0-1 matrix.
##
## Author: Qunqun
#########################################################################################################################
library(dplyr)
load("neural-network/Data/categoryEachDay.RData")

allCategory <- as.character(unique(as.vector(as.matrix(day_time_category_df))))
allCategory <- allCategory[2:13]

day_time_category_df <- as.data.frame(t(as.matrix(day_time_category_df)))
category_list <- lapply(day_time_category_df,function(x){
 category <-  sapply(allCategory,grepl, x, ignore.case = T)
 category_df <- as.data.frame(category)
 category_df$na <- apply(category_df,1,function(x) 1-sum(x))
 category_vector <- as.vector(t(as.matrix(category_df)))
})

category_matrix <- do.call(rbind,category_list)

save("category_matrix", file = "neural-network/Data/categoryMatrix.RData")

