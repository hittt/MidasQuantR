#-----------------------
# Function to get the initial parameters to start the optimization
#-----------------------

#' @export
GetIniParams <- function(y, X, X_neg, X_pos, q, numInitialsRand = 10000, numInitials = 10,beta2para = FALSE, As = FALSE){
  # Randomly sample second parameter of Beta polynomial
  InitialParamsVec = GetIniParamsC(y,X,X_neg, X_pos, q,numInitialsRand,numInitials,beta2para,As)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  if(As){
    beta = SortedResults[1:numInitials,2:(5 + beta2para)]
  } else{
    beta = SortedResults[1:numInitials,2:(4 + beta2para)]
  }
  return(beta)
}