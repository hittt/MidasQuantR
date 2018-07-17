#-----------------
# MIDAS Beta weighting function
#-----------------
#' @import dplyr
#' @import pracma
#' @export MidasBeta

#--------------------
# Function to compute MIDAS Beta Weights
#--------------------
MidasBeta <- function(nlag,param1,param2){
  eps = 2.2204e-16
  seq = linspace(eps,1-eps,nlag);

  if(param1 == 1){
    weight = (1-seq)^(param2-1)
  } else{
    weight = ((1-seq)^(param2-1)) * (seq^(param1-1))
  }

weight = weight/sum(weight,na.rm = TRUE)
return(weight)
}
