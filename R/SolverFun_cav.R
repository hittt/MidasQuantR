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
.sol_cav <- function(MainSolver,SecondSolver, betaIni, fun, control, y, x, As, q, empQuant, Uni = FALSE,
                     warn = TRUE,condMean = NULL){
  rep = control$rep
  control$rep = NULL
  N = NROW(betaIni)
  xsol = vector(mode="list", length = N)
  convCheck = 0;
  for(i in 1:N){
    for(ii in 1:rep){
      sol = .solverSwitch_cav(MainSolver, betaIni[i,], fun, control, y, x, As, empQuant, Uni, q,condMean)
      iniParsTemp = sol$par
      sol = .solverSwitch_cav(SecondSolver, iniParsTemp, fun, control, y, x, As, empQuant, Uni, q,condMean)
      iniParsTemp = sol$par
      sol = .solverSwitch_cav(MainSolver, iniParsTemp, fun, control, y, x, As, empQuant, Uni, q,condMean)
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

.solverSwitch_cav <- function(solver, pars, fun, control, y, x, As, empQuant, Uni, q,condMean){
  if(!is.na(match(solver,c("L-BFGS-B","Nelder-Mead")))){
    control$method = solver
    solver = "optim"
  }
  control$rep = NULL
  solution = switch(solver,
                 nmkb = .nmkbsolver_cav(pars, fun, control,  y, x, As, empQuant, Uni, q,condMean),
                 optim = .optimsolver_cav(pars, fun, control,  y, x, As, empQuant, Uni, q,condMean),
                 ucminf = .ucminfsolver_cav(pars, fun, control, y, x, As, empQuant, Uni, q,condMean),
                 nlminb = .nlminbsolver_cav(pars, fun, control, y, x, As, empQuant, Uni, q,condMean),
                 bobyqa = .bobyqasolver_cav(pars, fun, control, y, x, As, empQuant, Uni, q,condMean))
  return(solution)
}
#-----------------
# SOLVER MAIN FUNCTIONS
#-----------------

.nlminbsolver_cav = function (pars, fun, control, y, x, As, empQuant, Uni, q,condMean){
  control$method = NULL
  control = .nlminb.ctrl(control)
  rep = 10
  if(is.null(condMean)){
    ans = try(nlminb(start = pars, objective = fun, control = control, yr = y, Xr = x, 
                     As = As, empQuant = empQuant, Uni = Uni, q = q), silent = TRUE)
  } else{
    ans = try(nlminb(start = pars, objective = fun, control = control, yr = y, Xr = x, condmeanR = condMean,
                     As = As, empQuant = empQuant, Uni = Uni, q = q), silent = TRUE)
  }
  
  pscale = rep(1, length(pars))
  smin = 0.1
  maxtries = 1
  while(ans$convergence!=0 && maxtries < rep) {
    control$step.min = smin*0.1
    smin = smin*0.1
    pscale = 0.25*pscale
    if(is.null(condMean)){
      ans = try(nlminb(start = pars, objective = fun, control = control, yr = y, Xr = x, 
                       As = As, empQuant = empQuant, Uni = Uni, q = q), silent = TRUE)
    } else{
      ans = try(nlminb(start = pars, objective = fun, control = control, yr = y, Xr = x, condmeanR = condMean,
                       As = As, empQuant = empQuant, Uni = Uni, q = q), silent = TRUE)
    }
    maxtries = maxtries+1
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
  } else{
    sol = ans
    sol$value = sol$objective
    sol$objective = NULL
  }
  return(sol = sol)
}

.ucminfsolver_cav = function(pars, fun, control, y, x, As, empQuant, Uni, q,condMean){
  control = .ucminf.ctrl(control)
  if(is.null(condMean)){
    ans = try(ucminf(fn = fun, par = pars, control = control, yr = y, Xr = x, 
                     As = As, empQuant = empQuant, Uni = Uni, q = q), silent = TRUE)
  } else{
    ans = try(ucminf(fn = fun, par = pars, control = control, yr = y, Xr = x,condmeanR = condMean,
                     As = As, empQuant = empQuant, Uni = Uni, q = q), silent = TRUE)
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
    sol$par = rep(NA,length(pars))
  } else{
    sol = ans
    if(ans$convergence>0) sol$convergence = 0 else sol$convergence = 1
  }
return(sol)
}

.optimsolver_cav = function(pars,  fun, control, y, x, As, empQuant, Uni, q,condMean){
  method = control$method
  control$method = NULL
  if(method == "L-BFGS-B"){
    if(is.null(condMean)){
      ans = try(optim(fn = fun, par = pars, control = control, method = "L-BFGS-B",As = As,yr = y, Xr = x,gr = NULL,
                   empQuant = empQuant, Uni = Uni, q = q),silent = TRUE)
    }else{
      ans = try(optim(fn = fun, par = pars, control = control, method = "L-BFGS-B",As = As,yr = y, Xr = x,gr = NULL,
                      empQuant = empQuant, Uni = Uni, q = q,condmeanR = condMean),silent = TRUE)
    }
   } else {
     if(is.null(condMean)){
       ans = try(optim(fn = fun, par = pars, control = control, method = "Nelder-Mead",As = As,yr = y, Xr = x,gr = NULL,
                       empQuant = empQuant, Uni = Uni, q = q),silent = TRUE)
     }else{
       ans = try(optim(fn = fun, par = pars, control = control, method = "Nelder-Mead",As = As,yr = y, Xr = x,gr = NULL,
                       empQuant = empQuant, Uni = Uni, q = q,condmeanR = condMean),silent = TRUE)
     }
  }
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

.bobyqasolver_cav = function(pars, fun, control,  y, x, As, empQuant, Uni, q, condMean){
  control = .minqa.ctrl(control,pars)
  if(is.null(condMean)){
    ans = try(bobyqa(fn = fun, par = pars, control = control,yr = y, Xr = x,
                     As = As, empQuant = empQuant, Uni = Uni, q = q),silent = TRUE)
  } else{
    ans = try(bobyqa(fn = fun, par = pars, control = control,yr = y, Xr = x,
                     As = As, empQuant = empQuant, Uni = Uni, q = q,condmeanR = condMean),silent = TRUE)
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
  }  else{
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

.nmkbsolver_cav = function(pars,  fun, control, y, x, As, empQuant, Uni, q,condMean){
  control = .dfoptim.ctrl(control)
  if(is.null(condMean)){
    ans = try(nmkb(fn = fun, par = pars, control = control, yr = y, 
                 Xr = x, As = As, empQuant = empQuant, Uni = Uni, q = q), silent = TRUE)
  } else{
    ans = try(nmkb(fn = fun, par = pars, control = control, yr = y, 
                   Xr = x, As = As, empQuant = empQuant, Uni = Uni, q = q,condmeanR = condMean), silent = TRUE)
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
  }  else{
    sol = ans
  }
  return(sol)
}
