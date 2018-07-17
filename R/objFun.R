#---------------
# Quantile Loss Function
#---------------
#' @export objFun
#' @useDynLib MidasQuantR
#' @importFrom Rcpp sourceCpp evalCpp
#--------------
# Compute the value of loss function given parameters
#-------------
objFun <- function(params,y,X,q,beta2para = NULL,out = NULL){
  if(is.null(beta2para)) beta2para = FALSE
  if(is.null(out)) out = 1
  intercept = params[1];
  slope = params[2];
  if(beta2para){
    k1 = params[3]
    k2 = params[4]
  } else{
    k1 = 1
    k2 = params[3]
  }
  nlag = dim(X)[2]
  weights = MidasBeta(nlag,k1,k2)
  condQuantile = intercept + slope * (X %*% weights)
  loss = y - condQuantile
  fval = sum(loss * (q - as.numeric((loss <= 0))))
  if(out == 1){
    return(fval)
  } else{
    return(condQuantile)
  }
}

.onUnload <- function (libpath) { library.dynam.unload("MidasQuantR", libpath)}
