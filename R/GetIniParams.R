#-----------------------
# Function to get the initial parameters to start the optimization
#-----------------------

#' @export
GetIniParams <- function(y, X, X_neg, X_pos, q, numInitialsRand = 10000, 
                         numInitials = 10,beta2para = FALSE, As = FALSE){
  # Randomly sample second parameter of Beta polynomial
  InitialParamsVec = GetIniParamsUni(y,X,X_neg, X_pos, q,numInitialsRand,numInitials,beta2para,As)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(3 + As + beta2para + 1)]
  return(beta)
}

#' @export
GetIniParams_AL <- function(y, condMean, QuantEst, X, X_neg, X_pos, q,
                            numInitialsRand = 10000, numInitials = 10,
                            beta2para = FALSE, As = FALSE){
  # Randomly sample second parameter of Beta polynomial
  InitialParamsVec = GetIniParamsAL(yr = y, condmeanR = condMean, QuantEst = QuantEst,
                                    Xr = X, Xr_neg = X_neg,Xr_pos = X_pos, q = q,
                                    numInitialsRand = numInitialsRand, numInitials = numInitials,
                                    beta2para = beta2para, As = As)
  SortedResults = InitialParamsVec#[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(4 + As + beta2para + 1)]
  return(beta)
}
