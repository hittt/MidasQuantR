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
  InitialParamsVec = C_GetIniParams_midas(y,X,X_neg, X_pos, q,numInitialsRand,numInitials,beta2para,As)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(3 + As + beta2para + 1)]
  return(beta)
}

#------------------------------------------------------------
# Get initial parameters for the midas quantile with AL distribution
#------------------------------------------------------------

#' @export
GetIniParams_midasAL <- function(y, condMean, QuantEst, X, X_neg, X_pos, q,
                            numInitialsRand = 10000, numInitials = 10,
                            beta2para = FALSE, As = FALSE){
  # Randomly sample second parameter of Beta polynomial
  InitialParamsVec = C_GetIniParams_midasAL(yr = y, condmeanR = condMean, QuantEst = QuantEst,
                                    Xr = X, Xr_neg = X_neg,Xr_pos = X_pos, q = q,
                                    numInitialsRand = numInitialsRand, numInitials = numInitials,
                                    beta2para = beta2para, As = As)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(4 + As + beta2para + 1)]
  return(beta)
}

#-------------------------------------------------------------
# Get initial parameters for the CAViaR model
#-------------------------------------------------------------
#' @export

GetIniParams_cav <- function(y, x, q, model, empQuant, Uni = TRUE, numInitialsRand = 10000, numInitials = 10){
  InitialParamsVec = C_GetIniParams_cav(yr = y, Xr = x, q = q, model = model, empQuant = empQuant,
                                        Uni = Uni, numInitialsRand = numInitialsRand, numInitials = numInitials)
  SortedResults = InitialParamsVec[order(InitialParamsVec[,1]),]
  beta = SortedResults[1:numInitials,2:(2 + model + 1)]
  return(beta)
}
