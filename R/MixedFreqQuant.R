#-------------
# Mixed Frequency Function
#-------------
#' @export MixedFreqQuant
#' @useDynLib MidasQuantR
#' @import dplyr
#' @import pracma
#' @importFrom Rcpp sourceCpp evalCpp
#-------------
MixedFreqQuant <- function(DataY,DataYdate,DataX,DataXdate,xlag,period,ovlap = NULL){
  if(is.null(ovlap)) ovlap = TRUE
  nobs <- length(DataY)
  nobsShort = nobs-xlag-period+1
  DataYlowfreq = matrix(NA,ncol = 1, nrow = nobsShort)
  DataYDateLow = matrix(NA,ncol = 1, nrow = nobsShort)
  for(t in (xlag+1):(nobs - period + 1)){
  DataYlowfreq[t-xlag,1] = sum(DataY[t:(t+period-1)]);
  DataYDateLow[t-xlag,1] = DataYdate[t];
  }
  if(!ovlap){
    DataYlowfreq = DataYlowfreq[seq(1,nobsShort,period),1]
    DataYDateLow = DataYDateLow[seq(1,nobsShort,period),1]
  }
  # Set the start date and end date according to xlag, period and ylag
  minDateY = DataYDateLow[1]
  minDateX = DataXdate[xlag+1]
  if(minDateY > minDateX){
    estStart = as.numeric(minDateY)
  } else {
    estStart = as.numeric(minDateX)
  }
  maxDateY = tail(DataYDateLow,1);
  maxDateX = tail(DataXdate,1);
  if(maxDateY > maxDateX){
  estEnd = as.numeric(maxDateX)
  }else{
    estEnd = as.numeric(maxDateY)
  }
  # Construct Y data
  tol = 1e-10
  locStart = head(which(DataYDateLow >= estStart-tol),1)
  locEnd = head(which(DataYDateLow >= estEnd-tol),1)
  EstY = DataYlowfreq[locStart:locEnd]
  EstYdate = DataYDateLow[locStart:locEnd]

  # Construct X data
  nobsEst <- length(EstY)

  EstX = matrix(NA,nrow = nobsEst, ncol = xlag)
  EstXdate = matrix(NA,nrow = nobsEst, ncol = xlag)

  for(t in 1:nobsEst){
    loc = head(which(DataXdate >= EstYdate[t]-tol),1)
    if(length(loc) == 0){
      loc = length(DataXdate)
    } else if(loc > length(DataX)){
      nobsEst = t - 1
      EstY = EstY[1:nobsEst,]
      EstYdate = EstYdate[1:nobsEst,]
      EstX = EstX[1:nobsEst,]
      EstXdate = EstXdate[1:nobsEst,]
      maxDate = tail(EstYdate,1)
      warning(paste("Observations are further truncated to", as.Date.numeric(maxDate,origin = "1970-01-01"), sep = " "))
    } else {
    EstX[t,] = DataX[seq(loc-1,loc-xlag,by = -1)]
    EstXdate[t,] = DataXdate[seq(loc-1,loc-xlag,by = -1)]
    }
  }
  output <- list(EstY = EstY, EstYdate = as.Date.numeric(EstYdate,origin = "1970-01-01"),
                 EstX = EstX, EstXdate = as.Date.numeric(EstXdate,origin = "1970-01-01"),
                 EstStart = estStart, EstEnd = estEnd)
}

.onUnload <- function (libpath) { library.dynam.unload("MidasQuantR", libpath)}
