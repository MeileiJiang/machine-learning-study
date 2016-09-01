##################################################################################
## Soft_Thresholding_Operator.R
## 
## Author: Meilei
##################################################################################

soft_threshold = function(x, u){
  n = length(x)
  softx = NULL
  for (i in 1:n) {
    if(abs(x[i]) <= u){
      softx[i] = 0
    } else if(x[i] > u){
      softx[i] = x[i] - u
    } else {
      softx[i] = x[i] + u
    }
  }
  return(softx)
}