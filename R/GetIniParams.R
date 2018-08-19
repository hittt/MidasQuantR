#-----------------------
# Function to get the initial parameters to start the optimization
#-----------------------

#-------------------------------------------------------
# Get initial parameters for the univariate midas quantile
#-------------------------------------------------------

#' @export
GetIniParams_midas <- function(y, X, X_neg, X_pos, q, numInitialsRand = 10000, 
                         numInitials = 10,beta2para = FALSE, As = FALSE){
  # Randomly sample second parameter of Beta polynomial
  InitialParamsVec = C_GetIniParams_midas(y,X,X_neg, X_pos, q,numInitialsRand,beta2para,As)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(3 + As + beta2para + 1)]
  return(beta)
}

#------------------------------------------------------------
# Get initial parameters for the midas quantile with AL distribution
#------------------------------------------------------------

#' @export
GetIniParamsAL_midas <- function(y, condMean, QuantEst, X, X_neg, X_pos, q,
                            numInitialsRand = 10000, numInitials = 10,
                            beta2para = FALSE, As = FALSE){
  # Randomly sample second parameter of Beta polynomial
  InitialParamsVec = C_GetIniParamsAL_midas(yr = y,condmeanR = condMean,QuantEst = QuantEst,Xr = X,Xr_neg = X_neg,
                                            Xr_pos = X_pos,q = q,numInitialsRand = numInitialsRand,
                                            beta2para = beta2para,As = As)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(4 + As + beta2para + 1)]
  return(beta)
}

#-------------------------------------------------------------
# Get initial parameters for the CAViaR model
#-------------------------------------------------------------
#' @export

GetIniParams_cav <- function(y,x, q, As, empQuant, Uni = TRUE, numInitialsRand = 10000, numInitials = 10){
  InitialParamsVec = C_GetIniParams_cav(yr = y, Xr = x, q = q, As = As, empQuant = empQuant,
                                        Uni = Uni, numInitialsRand = numInitialsRand)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(2 + As + 2)]
  return(beta)
}

#------------------------------------------------------------
# Get initial parameters for the CAV quantile with AL distribution
#------------------------------------------------------------

#' @export
GetIniParamsAL_cav <- function(y, condMean, QuantEst, X, As, empQuant,q,
                                 numInitialsRand = 10000, numInitials = 10,
                                 Uni = TRUE){
  # Randomly sample second parameter of Beta polynomial
  InitialParamsVec = C_GetIniParamsAL_cav(yr = y,condmeanR = condMean,QuantEst = QuantEst,Xr = X,
                                          q = q,numInitialsRand = numInitialsRand,As = As, empQuant = empQuant, Uni = Uni)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(2 + As + 1 + Uni + 1)]
  return(beta)
}
