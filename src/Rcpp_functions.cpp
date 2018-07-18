#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// These above headers are neccessary to use RcppArmadillo

//' @export
// [[Rcpp::export]]
arma::mat MidasBetaC(int nlag, double param1, double param2){
  double eps = 2.2204e-16;
  arma::mat seq = arma::linspace(eps,1-eps,nlag);
  arma::mat weights = pow(1.0 - seq,param2 - 1.0) % pow(seq,param1 - 1.0);
  double sumWeights = arma::accu(weights);
  return  weights*(1.0/sumWeights);
}


//' @export
// [[Rcpp::export]]
double objFunC(Rcpp::NumericVector params,Rcpp::NumericVector yr, 
                Rcpp::NumericMatrix Xr, double q, bool beta2para = false){
  double intercept = params[0];
  double slope = params[1];
  double k1 = 1.0, k2 = 1.0;
  if(beta2para){
    k1 = params[2];
    k2 = params[3];  
  } else {
    k2 = params[2];
  }
  int nlag = Xr.ncol(),T = Xr.nrow();
  arma::mat X(Xr.begin(), T, nlag, false);
  arma::colvec y(yr.begin(), yr.size(), false);
  arma::mat weights = MidasBetaC(nlag,k1,k2);
  arma::colvec condQuantile = intercept + slope * (X*weights);
  arma::colvec loss = y - condQuantile;
  double fval = 0; 
  for(int i = 0; i < T; ++i){
    fval += loss[i] * (q - (loss[i] <= 0));
  }
  return fval;
}

//' @export
// [[Rcpp::export]]
arma::colvec condQuantileC(Rcpp::NumericVector params,Rcpp::NumericVector yr, 
               Rcpp::NumericMatrix Xr, bool beta2para = false){
  double intercept = params[0];
  double slope = params[1];
  double k1 = 1.0, k2 = 1.0;
  if(beta2para){
    k1 = params[2];
    k2 = params[3];  
  } else {
    k2 = params[2];
  }
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::mat weights = MidasBetaC(nlag,k1,k2);
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  arma::colvec condQuantile = intercept + slope * (X*weights);
  return condQuantile;
}

// [[Rcpp::export]]
NumericMatrix GetIniParamsC(Rcpp::NumericVector yr, Rcpp::NumericMatrix Xr, double q, 
                     int numInitialsRand = 10000, int numInitials = 10, bool beta2para = false){
  double T = Xr.nrow(), nlag = Xr.ncol();
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  arma::colvec nInitalBeta = 99 * arma::randu(numInitialsRand,1) + 1;
  // Generate initial guess for the second parameters of the Beta polynomial based on Uniform distribution
  NumericMatrix  InitialParamsVec(numInitialsRand,5);
  for(int i = 0; i < numInitialsRand;++i){
    arma::vec weights = MidasBetaC(nlag, 1 , nInitalBeta[i]);
    arma::mat Xweighted = X * weights;
    arma::mat X0 = arma::join_rows(arma::ones(T,1),Xweighted);
    arma::mat olsIni = solve(X0,y);
    if(beta2para){
      NumericVector xx = NumericVector::create(olsIni[0],olsIni[1],1,nInitalBeta[i]);
      double fval = objFunC(xx,yr,Xr,q,beta2para);
      InitialParamsVec(i,_) = NumericVector::create(fval,olsIni[0],olsIni[1],1,nInitalBeta[i]);
    } else{
      NumericVector xx = NumericVector::create(olsIni[0],olsIni[1],nInitalBeta[i]);
      double fval = objFunC(xx,yr,Xr,q,beta2para);
      InitialParamsVec(i,_) = NumericVector::create(fval,olsIni[0],olsIni[1],nInitalBeta[i]);  
    }
  }
  return InitialParamsVec;
}