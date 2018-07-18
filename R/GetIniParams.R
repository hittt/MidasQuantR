#-----------------------
# Function to get the initial parameters to start the optimization
#-----------------------

GetIniParams <- function(y, X, q, numInitialsRand = 10000, numInitials = 10,beta2para = FALSE){
  # Randomly sample second parameter of Beta polynomial
  InitialParamsVec = GetIniParamsC(y,X,q,numInitialsRand,numInitials,beta2para)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(4 + beta2para)]
  return(beta)
}