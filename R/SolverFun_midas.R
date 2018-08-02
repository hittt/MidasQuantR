#--------------------
# Utility function to solve the objFun
#--------------------
#' @importFrom ucminf ucminf
#' @importFrom stats optim
#' @importFrom stats nlminb
#' @importFrom minqa bobyqa
#' @importFrom dfoptim nmkb
# Decide the solver that is derivative-free. 
# CASE 1: The default case is holding the first parameters to be 1 and vary the second paramter of Beta Polynomial
# Good optimizer: Nelder-Mead (accurate but relatively slow) ; nlminb (good and quick); ucminf (good and quick); 
# BFGS (good and relatively quick)
# CASe 2: both parameters of the Beta Polynomials to vary. It seems that allowing for both two parameters of the 
# Beta polynomial to vary make the convergence very difficult.
# Relatively ok optimizer: L-BFGS-B; bobyqa; nlminb

# The optimization procedure requires 2 optimizers. For each candidate of initial parameters, the function first attempt
# to solve the problem using the MainSolver, then using the solution (converge or not) to solve the problem. The solution of
# the second solver is then used as initial parameters to resolve the problem. The use of multiple solvers is an attempt to 
# get the global optimization. The process is repeated 10 times over 10 inital paramters guess (default)
.sol <- function(MainSolver,SecondSolver, betaIni, fun, y, condMean = NULL, x, x_neg, x_pos,
                 q, beta2para, lb, ub, control, warn = TRUE, As = FALSE){
  rep = control$rep
  control$rep = NULL
  N = NROW(betaIni)
  xsol = vector(mode="list", length = N)
  convCheck = 0;
  for(i in 1:N){
    for(ii in 1:rep){
      sol = .solverSwitch(MainSolver, betaIni[i,], fun, control,  lb, ub, y, condMean, x,  x_neg, x_pos, q, beta2para, As)
      iniParsTemp = sol$par
      sol = .solverSwitch(SecondSolver, iniParsTemp, fun, control,  lb, ub, y, condMean, x,  x_neg, x_pos, q, beta2para, As)
      iniParsTemp = sol$par
      sol = .solverSwitch(MainSolver, iniParsTemp, fun, control,  lb, ub, y, condMean, x,  x_neg, x_pos, q, beta2para, As)
      if(sol$convergence == 0) break
    }
    if(sol$convergence == 0){
      xsol[[i]] = sol
    } else{
      convCheck = convCheck + 1
    }
  }
  if(convCheck == N){
    out = list()
    out$convergence = 1
    out$pars = rep(NA, dim(betaIni)[2])
    out$value = NA
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

.solverSwitch <- function(solver, pars, fun, control, lb, ub, y, condMean, x, x_neg, x_pos,
                       q, beta2para, As){
  control$rep = NULL
  if(!is.na(match(solver,c("L-BFGS-B","Nelder-Mead")))){
    control$method = solver
    solver = "optim"
  }
  solution = switch(solver,
                 nmkb = .nmkbsolver(pars, fun, control, lb, ub, y, condMean, x,x_neg,x_pos, q, beta2para,As),
                 optim = .optimsolver(pars, fun, control, lb, ub, y, condMean, x, x_neg,x_pos, q, beta2para,As),
                 ucminf = .ucminfsolver(pars, fun, control, y, condMean, x, x_neg,x_pos,q, beta2para,As),
                 nlminb = .nlminbsolver(pars, fun, control, lb, ub, y, condMean, x, x_neg,x_pos ,q, beta2para,As),
                 bobyqa = .bobyqasolver(pars, fun, control, lb, ub, y, condMean, x, x_neg,x_pos, q, beta2para,As))
  return(solution)
}
#-----------------
# SOLVER MAIN FUNCTIONS
#-----------------

.nlminbsolver = function (pars, fun, control, lb, ub, y, condMean, x, x_neg, x_pos, q, beta2para, As){
  control$method = NULL
  control = .nlminb.ctrl(control)
  rep = 10
  if(is.null(condMean)){
    ans = try(nlminb(start = pars, objective = fun, control = control,As = As,
                   yr = y, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, 
                   lower = lb, upper = ub), silent = TRUE)
  } else{
    ans = try(nlminb(start = pars, objective = fun, control = control,As = As,
                     yr = y, condmeanR = condMean, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, 
                     lower = lb, upper = ub), silent = TRUE)
  }
  pscale = rep(1, length(pars))
  smin = 0.1
  maxtries = 1
  while(ans$convergence!=0 && maxtries < rep) {
    control$step.min = smin*0.1
    smin = smin*0.1
    pscale = 0.25*pscale
    if(is.null(condMean)){
      ans = try(nlminb(start = pars, objective = fun, control = control,As = As,
                       yr = y, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, 
                       lower = lb, upper = ub), silent = TRUE)
    } else{
      ans = try(nlminb(start = pars, objective = fun, control = control,As = As,
                       yr = y, condmeanR = condMean, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, 
                       lower = lb, upper = ub), silent = TRUE)
    }
    maxtries = maxtries+1
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
  }
  else{
    sol = ans
    sol$value = sol$objective
    sol$objective = NULL
  }
  return(sol = sol)
}

.ucminfsolver = function(pars, fun, control,  y, condMean, x, x_neg, x_pos, q, beta2para, As){
  control = .ucminf.ctrl(control)
  if(is.null(condMean)){
    ans = try(ucminf(fn = fun, par = pars, control = control, yr = y, Xr = x, 
                   Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, As = As), silent = TRUE)
  } else{
    ans = try(ucminf(fn = fun, par = pars, control = control, yr = y, Xr = x, condmeanR = condMean,
                     Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, As = As), silent = TRUE)
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
    sol$par = rep(NA,length(pars))
  }
  else{
    sol = ans
    if(ans$convergence>0) sol$convergence = 0 else sol$convergence = 1
  }
return(sol)
}

.optimsolver = function(pars, fun, control, lb, ub, y, condMean, x,x_neg,x_pos, q, beta2para,As){
  method = control$method
  control$method = NULL
  if(method == "L-BFGS-B"){
    if(is.null(condMean)){
      ans = optim(fn = fun, par = pars, control = control, method = "L-BFGS-B",lower = lb, upper = ub,
                yr = y, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, As = As)
    } else {
      ans = optim(fn = fun, par = pars, control = control, method = "L-BFGS-B",lower = lb, upper = ub,condmeanR = condMean,
                  yr = y, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, As = As)
    }
  } else {
    if(is.null(condMean)){
      ans = optim(fn = fun, par = pars, control = control, method = "Nelder-Mead",lower = lb, upper = ub,
                  yr = y, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, As = As)
    } else {
      ans = optim(fn = fun, par = pars, control = control, method = "Nelder-Mead",lower = lb, upper = ub,condmeanR = condMean,
                  yr = y, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, As = As)
    }
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

.bobyqasolver = function(pars, fun, control, lb, ub, y, condMean, x, x_neg, x_pos, q, beta2para, As){
  control$method = NULL
  control = .minqa.ctrl(control,pars)
  if(is.null(condMean)){
    ans = bobyqa(fn = fun, par = pars, control = control,lower = lb, upper = ub,
               yr = y, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos,q = q, beta2para = beta2para, As = As)
  } else{
    ans = bobyqa(fn = fun, par = pars, control = control,lower = lb, upper = ub, condmeanR = condMean,
                 yr = y, Xr = x, Xr_neg = x_neg, Xr_pos = x_pos,q = q, beta2para = beta2para, As = As)
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
  }
  else{
    sol = ans
    sol$value = sol$fval
    sol$message = sol$msg
    sol$convergence = sol$ierr
    sol$fval = NULL
    sol$ierr= NULL
    sol$msg = NULL
  }
  return(sol)
}

.nmkbsolver = function(pars, fun, control, lb, ub, y, condMean, x, x_neg, x_pos, q, beta2para,As){
  control$method = NULL
  control = .dfoptim.ctrl(control)
  if(is.null(condMean)){
    ans = try(nmkb(fn = fun, par = pars, control = control, lower = lb, upper = ub, yr = y, 
                 Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, As = As), silent = TRUE)
  } else {
    ans = try(nmkb(fn = fun, par = pars, control = control, lower = lb, upper = ub, yr = y, condmeanR = condMean,
                   Xr = x, Xr_neg = x_neg, Xr_pos = x_pos, q = q, beta2para = beta2para, As = As), silent = TRUE)
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
  if(is.null(control$eval.max)) control$eval.max = 500
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

.dfoptim.ctrl = function(control){
  if(is.null(control$tol)) control$tol = 1e-08
  if(is.null(control$restarts.max)) control$restarts.max = 5
  mm = match(names(control), c("tol", "restarts.max"))
  if(any(is.na(mm))){
    idx = which(is.na(mm))
    wrong_opts = NULL
    for(i in 1:length(idx)) wrong_opts = c(wrong_opts, names(control)[idx[i]])
    warning(paste(c("\nunidentified option(s) in solver.control:\n", wrong_opts), sep="", collapse=" "), call. = FALSE, domain = NULL)
  }
  return(control)
}
################################################################################
