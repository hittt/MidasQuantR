#--------------------------------
# MAIN FUNCTIONS: Estimate MIDAS quantile regression
#--------------------------------

#' @export MidasQuantile
MidasQuantile <- function(y,yDate,x = NULL, xDate = NULL, q = 0.01,
                          horizon = 10, nlag = 100, ovlap = NULL, numInitialsRand = 10000,
                          numInitials = 10, GetSe = NULL, GetSeSim = NULL, Params = NULL, startPars = NULL,
                          MainSolver = "ucminf",SecondSolver = "nmkb",As = FALSE,
                          fitcontrol = list(rep = 5),beta2para = NULL,warn = TRUE){
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
  if(is.null(GetSe)) GetSe = TRUE
  if(is.null(GetSeSim)) GetSeSim = 200
  if(is.null(ovlap)) ovlap= TRUE
  if(is.null(beta2para)||!beta2para){
    beta2para = FALSE
    if(As){
      lb = c(-Inf,-Inf,-Inf,0)
      ub = c(Inf, Inf, Inf, 200)
    } else {
      lb = c(-Inf,-Inf,0)
      ub = c(Inf,Inf,200)
    }
  }else{
    if(As){
      lb = c(-Inf,-Inf,-Inf,0, 0)
      ub = c(Inf, Inf, Inf, 200, 200)
    } else {
      lb = c(-Inf,-Inf,0, 0)
      ub = c(Inf,Inf,200, 200)
    }
  }
  #------ Get the mixed data to start estimation -----
  dataEst <- MixedFreqQuant(DataY = y,DataYdate = yDate,DataX = x,DataXdate = xDate,
                            xlag = nlag,period = horizon,ovlap = ovlap)
  dataHigh <- MixedFreqQuant(DataY = y,DataYdate = yDate,DataX = y,DataXdate = yDate,
                             xlag = nlag,period = horizon,ovlap = ovlap)
  y = dataEst$EstY;
  yDate <- dataEst$EstYdate
  x = dataEst$EstX
  xDate = dataEst$EstXdate
  yHigh = dataHigh$EstX
  x_neg = x
  x_pos = x
  x_neg[yHigh > 0] = 0
  x_pos[yHigh <= 0] = 0
  
  #----- Get the initial guess for the parameters-----
  betaIni = GetIniParams(y = y, X = x,X_neg = x_neg,X_pos = x_pos, q = q, numInitialsRand = numInitialsRand,
                          numInitials = numInitials, beta2para = beta2para,As = As)
  
  #----- Estimate the paramters -----------
  sol = .sol(MainSolver = MainSolver,SecondSolver = SecondSolver,betaIni = betaIni,fun = objFun,
             y = y, x = x,x_neg = x_neg, x_pos = x_pos, q = q, beta2para = beta2para, 
             lb = lb, ub = ub, control = fitcontrol,warn = warn,As=As)
  estPars = sol$par
  fval = sol$value
  convergeFlag = sol$convergence
  if(convergeFlag == 1){
    warnings("\nBoth Solvers failed to converge, try with other available solvers...\n")
    out = list(estPars = estPars, pval = NA, yLowFreq = y, yDate = yDate, 
               condQuantile = NA, quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),
               fval = fval, conv = convergeFlag)
  } else{
    condQ = condQuantile(params = estPars, yr = y, Xr = x,Xr_neg = x_neg, Xr_pos = x_pos, beta2para = beta2para, As = As)
    betaIniSim = matrix(estPars,nrow = 1)
    if(beta2para){
      if(As){
        hypothesis = c(0,0,0,1,1)
      } else {
        hypothesis = c(0,0,1,1)
      }
    } else{
      if(As){
        hypothesis = c(0,0,0,1)
      } else {
        hypothesis = c(0,0,1)
      }
    }
    if(GetSe){
      resid = (y - condQ)/abs(condQ)
      parSim = matrix(0,nrow = length(estPars), ncol = GetSeSim)
      for(i in 1:GetSeSim){
        residSim = sample(resid,size = length(resid),replace = TRUE)
        ySim = condQ + residSim*abs(condQ)
        estSim = .solverSwitch(solver = MainSolver, pars = betaIniSim, y = ySim, x = x,x_neg = x_neg,
                               x_pos, fun = objFun,q = q, beta2para = beta2para, lb = lb, ub = ub,
                               control = fitcontrol, As = As)
        parSim[,i] = estSim$par
      }
      pval = rowMeans(abs(parSim - rowMeans(parSim,na.rm = TRUE) + hypothesis) > matrix(rep(abs(estPars),GetSeSim),nrow = length(estPars)),na.rm = TRUE)
    } else{
      pval = rep(NA,length(estPars))
    }
    out = list(estPars = estPars, pval = pval, yLowFreq = y, yDate = yDate, condQuantile = condQ,
               quantile = q, beta2para = beta2para, Solvers = c(MainSolver,SecondSolver),
               fval = fval, conv = convergeFlag)
  }
  return(out)
}

