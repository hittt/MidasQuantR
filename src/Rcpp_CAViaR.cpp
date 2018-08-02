#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// These above headers are neccessary to use RcppArmadillo

//' @export
// [[Rcpp::export]]
double CAViaRobjFun(arma::colvec params, Rcpp::NumericVector yr, Rcpp::NumericVector Xr,
                    double q, int model, double empQuant, bool Uni = true){
  int T = yr.size();
  Rcpp::NumericVector condQuant(T);
  Rcpp::NumericVector absYr = Rcpp::abs(yr);
  condQuant[0] = empQuant;
  if(Uni){
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condQuant[i] = beta0 + beta1 * absYr[i-1] + beta2 * condQuant[i-1] ;
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condQuant[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0)  + beta2 * condQuant[i-1];
      }
    }
  } else{
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Parameters for additional explanatory variable
      double beta3 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condQuant[i] = beta0 + beta1 * absYr[i-1] + beta2 * Xr[i-1] + beta3 * condQuant[i-1];
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Parameters for the additional explanatory variable
      double beta3 = params[4]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condQuant[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0) + beta2 * Xr[i-1] + beta3 * condQuant[i-1];
      }
    }
  }
  Rcpp::NumericVector loss = yr - condQuant;
  Rcpp::NumericVector fval(T);
  double out = 0;
  for(int i = 0; i < T; ++i){
    out = loss[i] * (q - (loss[i] <=0));
    if(Rcpp::traits::is_infinite<REALSXP>(out)){ // Using line search technique may lead to trials with infinite values
      fval[i] = 1e100;
    }else{
      fval[i] = out;
    }
  }
  return sum(fval);
  }

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector CAViaRcondQuant(Rcpp::NumericVector params, Rcpp::NumericVector yr, Rcpp::NumericVector Xr,
                                    double q, int model, double empQuant, bool Uni = true){
  int T = yr.size();
  Rcpp::NumericVector condQuant(T);
  Rcpp::NumericVector absYr = Rcpp::abs(yr);
  condQuant[0] = empQuant;
  if(Uni){
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condQuant[i] = beta0 + beta1 * absYr[i-1] + beta2 * condQuant[i-1] ;
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condQuant[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0)  + beta2 * condQuant[i-1];
      }
    }
  } else{
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Parameters for additional explanatory variable
      double beta3 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condQuant[i] = beta0 + beta1 * absYr[i-1] + beta2 * Xr[i-1] + beta3 * condQuant[i-1];
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Parameters for the additional explanatory variable
      double beta3 = params[4]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condQuant[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0) + beta2 * Xr[i-1] + beta3 * condQuant[i-1];
      }
    }
  }
  return condQuant;
}

//' @export
// [[Rcpp::export]]
NumericMatrix C_GetIniParams_cav(Rcpp::NumericVector yr, Rcpp::NumericVector Xr, double q, int model, 
                                 double empQuant, bool Uni = true, int numInitialsRand  = 10000, 
                                 int numInitials = 10){
  arma::colvec y(yr.begin(),yr.size(),false);
  arma::colvec X(Xr.begin(),Xr.size(),false);
  int numPars = 2;
  if(Uni){
    numPars = numPars + model;
  } else{
    numPars = numPars + model + 1;
  }
  NumericMatrix InitialParamesVec(numInitialsRand,numPars + 1);
  for(int i = 0; i < numInitialsRand; ++i){
    arma::colvec candidatePars = arma::join_cols(2 * arma::randu(numPars-1,1) - 1, arma::randu(1,1)); 
    // Generate random numbers uniformly distributed between (-1,1) and between (0, 1) for the autoregressive parameter
    double fval = CAViaRobjFun(candidatePars,yr,Xr,q,model,empQuant,Uni);
    NumericVector xx(numPars + 1);
    xx[0] = fval;
    std::copy(candidatePars.begin(),candidatePars.end(),xx.begin()+1);
    InitialParamesVec(i,_) = xx;
  }
  return InitialParamesVec;
}