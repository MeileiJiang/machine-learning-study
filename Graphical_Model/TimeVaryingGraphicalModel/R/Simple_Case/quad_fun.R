#############################################################
## quad_fun.R
## Quadratic-bump.
## Author: Meilei
#############################################################


quad_fun = function(x){
  if(x>= 0   && x<= 1/3) return(6*x^2);
  if(x>= 1/3 && x<= 2/3) return(1 - 12*(x - 1/2)^2);
  if(x>= 2/3 && x<= 1)   return(6*(1 - x)^2);
  return(-1);
}


# generate the function line ----------------------------------------------

# require(ggplot2)
# grid = seq(0, 1, by = 0.001)
# 
# y = rep(0, length(grid))
# 
# for(i in 1:length(grid)){
#   y[i] = quad_fun(grid[i])
# }
# 
# f.df = data.frame( y, grid)
# 
# ggplot(data = f.df, aes(x = grid, y = y)) +
#   geom_line(col = "blue") +
#   scale_y_continuous(limits = c(0,1)) +
#   scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
#   labs(x = "time", y = "", title = "B-Spline Quadratic Block")
