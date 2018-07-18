#-----------------------
# Function to get the initial parameters to start the optimization
#-----------------------

GetIniParams <- function(y, X, q, nlag, numInitialsRand = 10000, numInitials = 10){
  # Randomly sample second parameter of Beta polynomial
  InitialParamsVec = GetIniParamsC(y,X,q,numInitialsRand,numInitials)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1],decreasing = TRUE),]
  beta = SortedResults[1:numInitials,2:5]
  return(beta)
}