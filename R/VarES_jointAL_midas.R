#--------------------------------
# MAIN FUNCTIONS: Estimate MIDAS quantile regression
#--------------------------------
#' @importFrom forecast auto.arima Arima
#' @importFrom lmtest coeftest

#' @export VarEs_jointAL_midas
VarEs_jointAL_midas <- function(y,yDate,x = NULL, xDate = NULL, q = 0.01,
                          armaOrder = NULL, horizon = 10, nlag = 100, ovlap = FALSE, numInitialsRand = 10000,
                          numInitials = 10, GetSe = TRUE, GetSeSim = 200, Params = NULL, startPars = NULL,
                          MainSolver = "ucminf",SecondSolver = "nmkb",As = FALSE,
                          fitcontrol = list(rep = 5),beta2para = FALSE,warn = TRUE, simpleRet = FALSE){
  #-- set up arguments ----
  if(length(yDate) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  y[is.na(y)] = mean(y,na.rm = TRUE)
  if(is.null(x)) x = abs(y) # If regressor is null, using absolute returns
  if(is.null(xDate)) xDate = yDate;
  if(length(x) != length(y))  stop("\nMidasQuantile-->error: Length of y and X should be the same\n")
  if(is.na(match(MainSolver,c("L-BFGS-B","bobyqa","nlminb","ucminf","Nelder-Mead","nmkb")))){
    stop("\nMidasQuantile-->error: available solvers are L-BFGS-B,bobyqa,nlminb,ucminf,Nelder-Mead,nmkb... \n")
  }
  if(!is.null(SecondSolver)){
    if(is.na(match(SecondSolver,c("L-BFGS-B","bobyqa","nlminb","Nelder-Mead","nmkb")))){
      stop("\nMidasQuantile-->error: Available solvers are L-BFGS-B,bobyqa,nlminb,Nelder-Mead,nmkb... \n")
    }
  }
  if(!beta2para){
    if(As){
      lb = c(-Inf,-Inf,-Inf,0,-Inf)
      ub = c(Inf, Inf, Inf, 200, Inf)
    } else {
      lb = c(-Inf,-Inf,0,-Inf)
      ub = c(Inf,Inf,200, Inf)
    }
  }else{
    if(As){
      lb = c(-Inf,-Inf,-Inf,0, 0,-Inf)
      ub = c(Inf, Inf, Inf, 200, 200, Inf)
    } else {
      lb = c(-Inf,-Inf,0, 0,-Inf)
      ub = c(Inf,Inf,200, 200, Inf)
    }
  }
  if(is.null(startPars)){
    UniQuantEst = MidasQuantile(y = y, yDate = yDate, x = x, xDate = xDate, q = q,
                                horizon = horizon, nlag = nlag, ovlap = ovlap, numInitialsRand = numInitialsRand,
                                numInitials = numInitials, GetSe = FALSE, GetSeSim = NULL, Params = NULL, 
                                startPars = NULL, MainSolver = MainSolver, SecondSolver = SecondSolver, As = As,
                                fitcontrol = fitcontrol, beta2para = beta2para, warn = FALSE, simpleRet = simpleRet)
    if(UniQuantEst$conv == 1){
    stop("\nMidasQuantile -->error: The univariate quantile does not converge, try other solvers.\n", call. = FALSE)
    }
  } 
  #------ Get the mixed data to start estimation -----
  dataEst <- MixedFreqQuant(DataY = y,DataYdate = yDate,DataX = x,DataXdate = xDate,
                            xlag = nlag,period = horizon,ovlap = ovlap, simpleRet = simpleRet)
  dataHigh <- MixedFreqQuant(DataY = y,DataYdate = yDate,DataX = y,DataXdate = yDate,
                             xlag = nlag,period = horizon,ovlap = ovlap, simpleRet = simpleRet)
  y = dataEst$EstY;
  yDate <- dataEst$EstYdate
  x = dataEst$EstX
  xDate = dataEst$EstXdate
  yHigh = dataHigh$EstX
  x_neg = x
  x_pos = x
  x_neg[yHigh > 0] = 0
  x_pos[yHigh <= 0] = 0
  
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
  betaIni = GetIniParamsAL_midas(y = y, condMean = condMean, QuantEst = UniQuantEst$estPars, X = x, 
                         X_neg = x_neg,X_pos = x_pos, q = q, numInitialsRand = numInitialsRand,
                         numInitials = numInitials, beta2para = beta2para,As = As)
  
  #----- Estimate the paramters -----------
  sol = .sol(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFunAL_midas,
             condMean = condMean, y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, 
             lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As)
  estPars = sol$par
  fval = sol$value
  convergeFlag = sol$convergence
  if(convergeFlag == 1){
    warnings("\nBoth Solvers failed to converge, try with other available solvers...\n")
    out = list(estPars = estPars, pval = NA, yLowFreq = y, yDate = yDate, 
               condVaR = NA, condES = NA, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),
               fval = fval, conv = convergeFlag)
  } else{
    VaRES = condVaRES_midas(params = estPars, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos,beta2para = beta2para,As = As)
    betaIniSim = matrix(estPars,nrow = 1)
    condVaR = VaRES$VaR
    condES = VaRES$ES
    if(beta2para){
      if(As){
        hypothesis = c(0,0,0,1,1,0)
      } else {
        hypothesis = c(0,0,1,1,0)
      }
    } else{
      if(As){
        hypothesis = c(0,0,0,1,0)
      } else {
        hypothesis = c(0,0,1,0)
      }
    }
    if(GetSe){
      resid = (y - condVaR)/abs(condVaR)
      parSim = matrix(0,nrow = length(estPars), ncol = GetSeSim)
      for(i in 1:GetSeSim){
        residSim = sample(resid,size = length(resid),replace = TRUE)
        ySim = condVaR + residSim*abs(condVaR)
        ESsimMean = mean(residSim[residSim<0])
        Exceed = which(ySim < condVaR)
        ySim[Exceed] = condVaR[Exceed] + (residSim[Exceed]/ESsimMean)*(condES[Exceed]-condVaR[Exceed]);
        ySim_filter = forecast::Arima(ySim,order = meanFit$arma[1:3])
        condMeanSim <- as.numeric(ySim_filter$fitted)
        estSim = .solverSwitch(solver = MainSolver, pars = betaIniSim, y = ySim,condMean = condMeanSim, x = x,x_neg = x_neg,
                               x_pos = x_pos, fun = objFunAL_midas,q = q, beta2para = beta2para, lb = lb, ub = ub,
                               control = fitcontrol, As = As)
        parSim[,i] = estSim$par
      }
      pval = rowMeans(abs(parSim - rowMeans(parSim,na.rm = TRUE) + hypothesis) > matrix(rep(abs(estPars),GetSeSim),nrow = length(estPars)),na.rm = TRUE)
    } else{
      pval = rep(NA,length(estPars))
    }
    out = list(estPars = estPars, pval = pval, yLowFreq = y, yDate = yDate, condVaR = condVaR,condES = condES,
               quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),
               fval = fval, conv = convergeFlag, meanEst = meanCoef)
  }
  return(out)
}

