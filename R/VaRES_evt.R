#--------------------------------
# MAIN FUNCTIONS: Estimate MIDAS quantile regression
#--------------------------------
#' @importFrom fExtremes gpdFit
#' @export VarEs_evt
VarEs_evt <- function(y,yDate,x = NULL, xDate = NULL, q = 0.01,qThreshold = 0.075,
                      horizon = 10, nlag = 100, ovlap = FALSE, numInitialsRand = 10000,
                      numInitials = 10, GetSe = FALSE, GetSeSim = 200, Params = NULL, startPars = NULL,
                      MainSolver = "ucminf",SecondSolver = "nmkb",quantEst = "midas", As = FALSE,
                      fitcontrol = list(rep = 5),beta2para = FALSE,warn = TRUE, simpleRet = FALSE,Uni = TRUE){
  #-- set up arguments ----
  if(length(yDate) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  y[is.na(y)] = mean(y,na.rm = TRUE)
  if(is.null(x)) x = abs(y) # If regressor is null, using absolute returns
  if(is.null(xDate)) xDate = yDate;
  if(length(x) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  if(is.na(match(quantEst,c("midas","cav")))){
    stop("\nError: the quantile regression should be either midas or cav..\n")
  }
  if(is.na(match(MainSolver,c("L-BFGS-B","bobyqa","nlminb","ucminf","Nelder-Mead","nmkb")))){
    stop("\nMidasQuantile-->error: available solvers are L-BFGS-B,bobyqa,nlminb,ucminf,Nelder-Mead,nmkb... \n")
  }
  if(!is.null(SecondSolver)){
    if(is.na(match(SecondSolver,c("L-BFGS-B","bobyqa","nlminb","Nelder-Mead","nmkb")))){
      stop("\nMidasQuantile-->error: Available solvers are L-BFGS-B,bobyqa,nlminb,Nelder-Mead,nmkb... \n")
    }
  }
  #------- In case of known parameters, just compute the condVaR and condEs and return
  if(!is.null(Params)){
    n = length(Params)
    quantPars = Params[1:(n-2)]
    gamma = Params[n-1]; beta = Params[n]
    VaRthreshold <- switch (quantEst,
                            cav = CAViaR(y = y, yDate = yDate, x = x, q = qThreshold, horizon = horizon, ovlap = ovlap, numInitialsRand = numInitialsRand,
                                         numInitials = numInitials, empQuant = empQuant, GetSe = FALSE, Params = quantPars, startPars = startPars,
                                         MainSolver = MainSolver,SecondSolver = SecondSolver, As = As, fitcontrol = list(rep = 5),
                                         warn = warn, simpleRet = simpleRet, Uni = Uni),
                            midas = MidasQuantile(y = y, yDate = yDate, x = x, xDate = xDate, q = qThreshold, horizon = horizon, nlag = nlag, ovlap = ovlap,
                                                  numInitialsRand = numInitialsRand, numInitials = numInitials, GetSe = FALSE, GetSeSim = GetSeSim, Params = quantPars,
                                                  startPars = startPars, MainSolver = MainSolver, SecondSolver = SecondSolver, As = As, fitcontrol = fitcontrol,
                                                  beta2para = beta2para, warn = warn, simpleRet = simpleRet)
    )
    condThres = VaRthreshold$condVaR
    y = VaRthreshold$yLowFreq
    StdExceed = y/condThres - 1
    Nexceed = length(StdExceed[StdExceed>0])
    nobs = length(y)
    StdQuantile = ((q*(nobs/Nexceed))^(-gamma)-1)*(beta/gamma)
    condVaR = condThres * (1 + StdQuantile)
    ExpectedMean_StdQuant = StdQuantile * (1/(1-gamma) + beta/((1-gamma)*StdQuantile))
    condES = condThres * (1 + ExpectedMean_StdQuant)
  } else{
  #------ Get the threshold for the extreme value theory-----
  VaRthreshold <- switch (quantEst,
    cav = CAViaR(y = y, yDate = yDate, x = x, q = qThreshold, horizon = horizon, ovlap = ovlap, numInitialsRand = numInitialsRand,
                   numInitials = numInitials, empQuant = empQuant, GetSe = GetSe, Params = NULL, startPars = startPars,
                   MainSolver = MainSolver,SecondSolver = SecondSolver, As = As, fitcontrol = list(rep = 5),
                   warn = warn, simpleRet = simpleRet, Uni = Uni),
    midas = MidasQuantile(y = y, yDate = yDate, x = x, xDate = xDate, q = qThreshold, horizon = horizon, nlag = nlag, ovlap = ovlap,
                          numInitialsRand = numInitialsRand, numInitials = numInitials, GetSe = GetSe, GetSeSim = GetSeSim, Params = NULL,
                          startPars = startPars, MainSolver = MainSolver, SecondSolver = SecondSolver, As = As, fitcontrol = fitcontrol,
                          beta2para = beta2para, warn = warn, simpleRet = simpleRet)
  )
  
  #------ Estimate the General Pareto Distribution parameters----
  if(VaRthreshold$conv == 1){
    stop("\nThe quantile regression for threshold level is not converged, refit the model with other solver...\n")
  } else{
    condThres = VaRthreshold$condVaR
    y = VaRthreshold$yLowFreq
    StdExceed = y/condThres - 1
    GPDfit = fExtremes::gpdFit(x = StdExceed,u = 0)
    GPDpars = GPDfit@fit$par.ests
    gamma = GPDpars[1]; beta = GPDpars[2]
    Nexceed = length(StdExceed[StdExceed>0])
    nobs = length(y)
    StdQuantile = ((q*(nobs/Nexceed))^(-gamma)-1)*(beta/gamma)
    condVaR = condThres * (1 + StdQuantile)
    ExpectedMean_StdQuant = StdQuantile * (1/(1-gamma) + beta/((1-gamma)*StdQuantile))
    condES = condThres * (1 + ExpectedMean_StdQuant)
  }
  }
  if(gamma >= 1 || any(condVaR >= 0)){
    stop("\n Unrealistic estimate of GPD paramters or positive condVaR, reestimate with another threshold level..\n")
  }
  out = list(thresholdEst = VaRthreshold, yLowFreq = y, yDate = yDate, condVaR = condVaR,condES = condES,
             GPDest = GPDfit, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),
             As = As, Uni = Uni)
  return(out)
}

