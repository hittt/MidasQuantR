#--------------------------------
# MAIN FUNCTIONS: Estimate MIDAS quantile regression
#--------------------------------
#' @importFrom forecast auto.arima Arima
#' @importFrom lmtest coeftest

#' @export VarEs_jointAL_cav
  VarEs_jointAL_cav <- function(y,yDate,x = NULL, xDate = NULL, q = 0.01,
                       armaOrder = NULL, horizon = 10, ovlap = FALSE, numInitialsRand = 10000,
                       numInitials = 10, GetSe = TRUE, GetSeSim = 200, Params = NULL, startPars = NULL,
                       MainSolver = "ucminf",SecondSolver = "Nelder-Mead",model = 1,Uni = TRUE,empQuant = NULL,
                       fitcontrol = list(rep = 5),warn = TRUE, simpleRet = FALSE,constrained = FALSE){
  nobs <- length(y)
  nobsShort = nobs-horizon+1
  yLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
  yDateLowFreq = matrix(NA,ncol = 1, nrow = nobsShort)
  for(t in 1:(nobs - horizon + 1)){
    yLowFreq[t,1] = sum(y[t:(t+horizon-1)]);
    yDateLowFreq[t,1] = yDate[t];
  }
  if(!ovlap){
    yLowFreq = yLowFreq[seq(1,nobsShort,horizon),1]
    yDateLowFreq = yDateLowFreq[seq(1,nobsShort,horizon),1]
  }
  if(simpleRet){
    y = exp(yLowFreq) - 1
  } else{
    y = yLowFreq
  }
  yDate = as.Date.numeric(yDateLowFreq,origin = "1970-01-01")
  if(is.null(x)){
    x = abs(y)
  }else {
    if(length(x) != length(y)){
      stop("\nCAviaR-->error: The explanatory variable should have the same length as the objective variable... \n")
    }
  }
  if(is.null(xDate)) xDate = yDate
  #-----Check the solvers input------
  # If we want to put constrains on the parameters, only the solvers that allow contrains is accepted.
  if(constrained){
    if(is.na(match(MainSolver,c("L-BFGS-B","bobyqa","nlminb","nmkb")))){
      stop("\nMidasQuantile-->error: only solver that allows for bounds can be used in case of constrained... \n")
    }
    if(is.na(match(SecondSolver,c("L-BFGS-B","bobyqa","nlminb","nmkb")))){
      stop("\nMidasQuantile-->error: only solver that allows for bounds can be used in case of constrained... \n")
    }
  }
  
  # The CAViaR model is sensitive to the choice of the empirical quantile to start the dynamics. Here, we use the empirical
  # quantile of the first 10% of the data sample to start the quantile dynamics.
  if(is.null(empQuant)) empQuant = unname(quantile(y[1:round(0.10*length(y))],q))
  
  # Set the bounds for the parameters. The autoregressive paramter should be between 0 and 1?
  tol = 1e-10
  if(constrained){
    if(model == 1){
      if(Uni){
        lb = c(-Inf,-Inf,tol,-Inf)
        ub = c(Inf,Inf,1-tol,Inf)
      } else{
        lb = c(-Inf,-Inf,-Inf,tol,-Inf)
        ub = c(Inf, Inf, Inf, 1-tol,Inf)
      }
    } else{
      if(Uni){
        lb = c(-Inf,-Inf,-Inf,tol,-Inf)
        ub = c(Inf, Inf, Inf, 1-tol,Inf)
      } else{
        lb = c(-Inf,-Inf,-Inf,-Inf,tol,-Inf)
        ub = c(Inf, Inf, Inf, Inf, 1-tol,Inf)
      }
    }
  } else{
    lb = -Inf
    ub = Inf
  }
  if(is.null(startPars)){
    UniQuantEst = CAViaR(y = y,yDate = yDate,x = x,xDate = xDate,q = q,horizon = 1,ovlap = FALSE,numInitialsRand = numInitialsRand,
                         numInitials = numInitials,empQuant = empQuant,GetSe = FALSE,Params = NULL,startPars = NULL,MainSolver = MainSolver,
                         SecondSolver = SecondSolver,model = model,fitcontrol = fitcontrol,warn = FALSE,simpleRet = simpleRet,
                        constrained = constrained,Uni = Uni)
    if(UniQuantEst$conv == 1){
      stop("\nCAViaR -->error: The univariate quantile does not converge, try other solvers.\n", call. = FALSE)
    }
  } 

  #------- Fit the conditional mean equation-------------
  if(is.null(armaOrder)){
    meanFit <- forecast::auto.arima(y, max.d = 0, max.D = 0)
    if(length(meanFit$coef) == 0){
      meanCoef = 0
      condMean = rep(0,length(y))
    } else{
      meanCoef <- lmtest::coeftest(meanFit)
      condMean = as.numeric(meanFit$fitted)
    }
  } else{
    meanFit = forecast::Arima(y,order = armaOrder)
    meanEst <- lmtest::coeftest(meanFit)
    meanCoef = unname(meanEst[,1])
    condMean = as.numeric(meanFit$fitted)
  }
  
  
  #----- Get the initial guess for the parameters-----
  betaIni = GetIniParamsAL_cav(y = y, condMean = condMean, QuantEst = UniQuantEst$estPars, X = x, 
                            q = q, numInitialsRand = numInitialsRand, model = model, empQuant = empQuant,
                            Uni = Uni)
  
  #----- Estimate the paramters -----------
  sol = .sol_cav(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFunAL_cav,
                 control = fitcontrol,lb = lb,ub = ub,y = y,x = x,model = model,q = q,empQuant = empQuant,
                 Uni = Uni,condMean = condMean)
  estPars = sol$par
  fval = sol$value
  convergeFlag = sol$convergence
  if(convergeFlag == 1){
    warnings("\nBoth Solvers failed to converge, try with other available solvers...\n")
    out = list(estPars = estPars, pval = NA, y = y, yDate = yDate, condES = NA,
               condVaR = NA, quantile = q, Uni = Uni, Solvers = c(MainSolver,SecondSolver),
               fval = fval, conv = convergeFlag, simpleRet = simpleRet,model = model)
  } else{
    VaRES = condVaRES_cav(params = estPars, yr = y, Xr = x, condmeanR = condMean, empQuant = empQuant, model = model, Uni = Uni)
    betaIniSim = matrix(estPars,nrow = 1)
    condVaR = VaRES$VaR
    condES = VaRES$ES
    if(Uni){
      if(model == 1){
        hypothesis = c(0,0,0,0) # intercept, absolute lag returns, autoregressive, phi
      } else{
        hypothesis = c(0,0,0,0,0)
      }
    } else{
      if(model == 1){
        hypothesis = c(0,0,0,0,0)
      } else{
        hypothesis = c(0,0,0,0,0,0)
      }
    }
    
    if(GetSe){
      resid = (y - condVaR)/abs(condVaR)
      parSim = matrix(0,nrow = length(estPars), ncol = GetSeSim)
      for(i in 1:GetSeSim){
        residSim = sample(resid,size = length(resid),replace = TRUE)
        xSim = sample(x, size = length(x),replace = TRUE)
        NegResidMean = mean(residSim[residSim < 0])
        ySim = cavSim(estPars, residSim, NegResidMean, y, condVaR, condES, model, Uni, xSim)
        ySim_filter = forecast::Arima(ySim,order = meanFit$arma[1:3])
        condMeanSim <- as.numeric(ySim_filter$fitted)
        estSim = .solverSwitch_cav(solver = MainSolver, pars = estPars, fun = objFunAL_cav, control = fitcontrol,
                                   model = model, lb = lb, ub = ub, y = ySim, x = xSim, empQuant = empQuant, Uni = Uni,
                                   condMean = condMeanSim, q = q)
        parSim[,i] = estSim$par
      }
      pval = rowMeans(abs(parSim - rowMeans(parSim,na.rm = TRUE) + hypothesis) > matrix(rep(abs(estPars),GetSeSim),nrow = length(estPars)),na.rm = TRUE)
    } else{
      pval = rep(NA,length(estPars))
    }
    out = list(estPars = estPars, pval = pval, y = y, yDate = yDate, condES = condES,
               condVaR = condVaR, quantile = q, Uni = Uni, Solvers = c(MainSolver,SecondSolver),
               fval = fval, conv = convergeFlag, simpleRet = simpleRet,model = model)
  }
  return(out)
}

