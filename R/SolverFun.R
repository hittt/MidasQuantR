#--------------------
# Utility function to solve the objFun
#--------------------
#' @importFrom ucminf ucminf
#' @importFrom stats optim
#' @importFrom stats nlminb
#' @importFrom minqa bobyqa
# Decide the solver that is derivative-free. 
# CASE 1: The default case is holding the first parameters to be 1 and vary the second paramter of Beta Polynomial
# Good optimizer: Nelder-Mead (accurate but relatively slow) ; nlminb (good and quick); ucminf (good and quick); 
# BFGS (good and relatively quick)
# CASe 2: both parameters of the Beta Polynomials to vary. It seems that allowing for both two parameters of the 
# Beta polynomial to vary make the convergence very difficult.
# Relatively ok optimizer: L-BFGS-B; bobyqa; nlminb

# The optimization routine is to apply the MainSolver to solver for the problem without bounds. 
# If the solution of the MainSolver is not converged --> use the second method with bounds.
# If the SecondSolver still cannot solve the problem, save everything and then produce a warnings to reestimate with other solvers


.sol <- function(MainSolver,SecondSolver,betaIni,fun,y,X,q,beta2para,lb,ub,control){
  if(!is.na(match(MainSolver,c("L-BFGS-B","Nelder-Mead")))){
    control$method = MainSolver
    MainSolver = "optim"
  }
  if(!is.na(match(SecondSolver,c("L-BFGS-B","Nelder-Mead")))){
    control$method = SecondSolver
    SecondSolver = "optim"
  }
  numPars = dim(betaIni)[2]
  retval = switch(MainSolver,
                  optim = .msoptim(betaIni, fun, control, lb, ub, y, x, q, beta2para),
                  ucminf = .msucminf(betaIni, fun, control, y, x, q, beta2para),
                  nlminb = .msnlminb(betaIni, fun, control, lb, ub, y, x, q, beta2para),
                  bobyqa = .msbobyqa(betaIni, fun, control, lb, ub, y, x, q, beta2para))
  # Using the Second Optimizer in case the first optimizer does not converge or violate the bounds
  
  # In case of non-converge: Reoptimize with the SecondSolver completely
  if(retval$convergence == 1){
    warning("\nMidasQuantile : MainSolver fails, switching to the SecondSolver..\n")
    retval = switch(SecondSolver,
                     optim = .msoptim(betaIni, fun, control, lb, ub, y, x, q, beta2para),
                     ucminf = .msucminf(betaIni, fun, control, y, x, q, beta2para),
                     nlminb = .msnlminb(betaIni, fun, control, lb, ub, y, x, q, beta2para),
                     bobyqa = .msbobyqa(betaIni, fun, control, lb, ub, y, x, q, beta2para))
    if(retval$convergence == 1){
      warning("\nMidasQuantile : Both Solvers fail..\n")
    }
  }
  # In case of converge but violate bounds, reoptimize with the second optimizer with the 
  # stating parameters inherit from the first round
  if(sum(retval$par > lb) != numPars ||  sum(retval$par < ub) != numPars){
    warning("\nMidasQuantile : MainSolver convergence violate bounds, switching to the SecondSolver..\n")
    for(i in 1:dim(betaIni[1])) betaIni[i,] = retval$par
    retval = switch(SecondSolver,
                    optim = .msoptim(betaIni, fun, control, lb, ub, y, x, q, beta2para),
                    nlminb = .msnlminb(betaIni, fun, control, lb, ub, y, x, q, beta2para),
                    bobyqa = .msbobyqa(betaIni, fun, control, lb, ub, y, x, q, beta2para))
    if(retval$convergence == 1){
      warning("\nMidasQuantile : Second solver fails to converge..\n")
    }
  }
  return(retval)
}

#-----------------
# SOLVER MAIN FUNCTIONS
#-----------------

.nlminbsolver = function (pars, fun, control, lb, ub, y, x, q, beta2para){
  control$method = NULL
  control = .nlminb.ctrl(control)
  rep = 10
  ans = try(nlminb(start = pars, objective = fun, control = control, 
                   yr = y, Xr = x, q = q, beta2para = beta2para, 
                   lower = lb, upper = ub), silent = TRUE)
  pscale = rep(1, length(pars))
  smin = 0.1
  maxtries = 1
  while(ans$convergence!=0 && maxtries < rep) {
    control$step.min = smin*0.1
    smin = smin*0.1
    pscale = 0.25*pscale
    ans = try(nlminb(start = pars, objective = fun, control = control, 
                     yr = y, Xr = x, q = q, beta2para = beta2para, 
                     lower = lb, upper = ub), silent = TRUE)
    maxtries = maxtries+1
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
  }
  else{
    sol = ans
  }
  return(sol = sol)
}

.ucminfsolver = function(pars, fun, control,  y, x, q, beta2para){
  control$method = NULL
  control = .ucminf.ctrl(control)
  ans = try(ucminf(fn = fun, par = pars, control = control, yr = y, Xr = x, q = q, beta2para = beta2para), silent = TRUE)
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
  }
  else{
    sol = ans
    if(ans$convergence>0) sol$convergence = 0 else sol$convergence = 1
  }
return(sol)
}

.optimsolver = function(pars, fun, control, lb, ub, y, x, q, beta2para){
  method = control$method
  control$method = NULL
  if(method == "L-BFGS-B"){
    ans = optim(fn = fun, par = pars, control = control, method = "L-BFGS-B",lower = lb, upper = ub,
                yr = y, Xr = x, q = q, beta2para = beta2para)
  } else {
    ans = optim(fn = fun, par = pars, control = control, method = "Nelder-Mead",
                yr = y, Xr = x, q = q, beta2para = beta2para)
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
  }
  else{
    sol = ans
  }
  return(sol)
}

.bobyqasolver = function(pars, fun, control, lb, ub, y, x, q, beta2para){
  control$method = NULL
  control = .minqa.ctrl(control,pars)
  ans = bobyqa(fn = fun, par = pars, control = control,lower = lb, upper = ub,
               yr = y, Xr = x, q = q, beta2para = beta2para)
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
  }
  else{
    sol = ans
    sol$value = sol$fval
    sol$message = sol$msg
    sol$convergence = ans$ierr
    sol$fval = NULL
    sol$ierr= NULL
    sol$msg = NULL
  }
  return(sol)
}

.msnlminb = function(pars, fun, control, lb, ub, y, x, q, beta2para){
  N = NROW(pars)
  xsol = vector(mode="list", length = N)
  convCheck = 0;
  for(i in 1:N){
    sol = .nlminbsolver(pars[i,], fun, control,  lb, ub, y, x, q, beta2para)
    if(sol$convergence == 0){
      xsol[[i]] = sol
    } else{
      convCheck = convCheck + 1
    }
  }
  if(convCheck == N){
    out = list()
    out$convergence = 1
    out$pars = rep(NA, length(pars))
    out$value = NA
    warning("\nMidasQuantile-->warning: no convergence using nlminb...\n")
  } else{
    best = sapply(xsol, function(x) 
      if(is.null(x)){
        NA} else {
          x$objective
        })
    best = which(best == min(best, na.rm=TRUE))
    out = xsol[[best]]
  }
  return(out)
}

.msucminf = function(pars, fun, control, y, x, q, beta2para){
  N = NROW(pars)
  xsol = vector(mode="list", length = N)
  convCheck = 0;
  for(i in 1:N){
    sol = .ucminfsolver(pars[i,], fun, control, y, x, q, beta2para)
    if(sol$convergence == 0){
      xsol[[i]] = sol
    } else{
      convCheck = convCheck + 1
    }
  }
  if(convCheck == N){
    out = list()
    out$convergence = 1
    out$pars = rep(NA, length(pars))
    out$value = NA
    warning("\nMidasQuantile-->warning: no convergence using ucminf...\n")
  } else{
    best = sapply(xsol, function(x) 
      if(is.null(x)){
        NA} else {
          x$value
        })
    best = which(best == min(best, na.rm=TRUE))
    out = xsol[[best]]
  }
  return(out)
}

.msoptim = function(pars, fun, control, lb, ub, y, x, q, beta2para){
  N = NROW(pars)
  xsol = vector(mode="list", length = N)
  convCheck = 0;
  for(i in 1:N){
    sol = .optimsolver(pars[i,], fun, control,  lb, ub, y, x, q, beta2para)
    if(sol$convergence == 0){
      xsol[[i]] = sol
    } else{
      convCheck = convCheck + 1
    }
  }
  if(convCheck == N){
    out = list()
    out$convergence = 1
    out$pars = rep(NA, length(pars))
    out$value = NA
    warning(paste("\nMidasQuantile->warning: no convergence in ", method, "..\n", sep = " "))
  } else{
    best = sapply(xsol, function(x) 
      if(is.null(x)){
        NA} else {
          x$value
        })
    best = which(best == min(best, na.rm=TRUE))
    out = xsol[[best]]
  }
  return(out)
}

.msbobyqa = function(pars, fun, control, lb, ub, y, x, q, beta2para){
  N = NROW(pars)
  xsol = vector(mode="list", length = N)
  convCheck = 0;
  for(i in 1:N){
    sol = .bobyqasolver(pars[i,], fun, control,  lb, ub, y, x, q, beta2para)
    if(sol$convergence == 0){
      xsol[[i]] = sol
    } else{
      convCheck = convCheck + 1
    }
  }
  if(convCheck == N){
    out = list()
    out$convergence = 1
    out$pars = rep(NA, length(pars))
    out$value = NA
    warning("\nMidasQuantile-->warning: no convergence using bobyqa...\n")
  } else{
    best = sapply(xsol, function(x) 
      if(is.null(x)){
        NA} else {
          x$fval
        })
    best = which(best == min(best, na.rm=TRUE))
    out = xsol[[best]]
  }
  return(out)
}
#######################################
# SOLVER CONTROLS
#######################################
# Solver control parameters
.ucminf.ctrl = function(control){
  if(is.null(control$trace)) control$trace = 0
  if(is.null(control$xtol)) control$xtol = 1e-12
  if(is.null(control$stepmax)) control$stepmax = 1
  if(is.null(control$maxeval)) control$maxeval = 1500
  mm = match(names(control), c("trace", "xtol", "stepmax", "maxeval"))
  if(any(is.na(mm))){
    idx = which(is.na(mm))
    wrong_opts = NULL
    for(i in 1:length(idx)) wrong_opts = c(wrong_opts, names(control)[idx[i]])
    warning(paste(c("\nunidentified option(s) in solver.control:\n", wrong_opts), sep="", collapse=" "), call. = FALSE, domain = NULL)
  }
  return(control)
}

.nlminb.ctrl = function(control){
  if(is.null(control$trace)) control$trace = 0
  if(is.null(control$eval.max)) control$eval.max = 1500
  if(is.null(control$iter.max)) control$iter.max = 500
  if(is.null(control$abs.tol)) control$abs.tol = 0
  if(is.null(control$rel.tol)) control$rel.tol = 1e-10
  if(is.null(control$x.tol)) control$x.tol = 2.2e-8
  if(is.null(control$xf.tol)) control$xf.tol = 2.2e-14
  if(is.null(control$step.min)) control$step.min = 0.1
  if(is.null(control$step.max)) control$step.max = 1
  mm = match(names(control), c("trace", "eval.max", "iter.max", "abs.tol", "rel.tol", "x.tol", "xf.tol",
                               "step.min", "step.max"))
  if(any(is.na(mm))){
    idx = which(is.na(mm))
    wrong_opts = NULL
    for(i in 1:length(idx)) wrong_opts = c(wrong_opts, names(control)[idx[i]])
    warning(paste(c("\nunidentified option(s) in solver.control:\n", wrong_opts), sep="", collapse=" "), call. = FALSE, domain = NULL)
  }
  return(control)
}

.minqa.ctrl = function(control, pars){
  n = length(pars)
  if(is.null(control$npt)) control$npt = min(n*2, n+2)
  if(is.null(control$iprint)) control$iprint = 0
  if(is.null(control$maxfun)) control$maxfun = 10000
  mm = match(names(control), c("npt", "rhobeg", "rhoend", "iprint", "maxfun"))
  if(any(is.na(mm))){
    idx = which(is.na(mm))
    wrong_opts = NULL
    for(i in 1:length(idx)) wrong_opts = c(wrong_opts, names(control)[idx[i]])
    warning(paste(c("\nunidentified option(s) in solver.control:\n", wrong_opts), sep="", collapse=" "), call. = FALSE, domain = NULL)
  }
  return(control)
}
################################################################################
