#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// These above headers are neccessary to use RcppArmadillo


// [[Rcpp::export]]
arma::vec MidasBetaC(int nlag, double param1, double param2){
  double eps = 2.2204e-16;
  arma::vec seq = arma::linspace(eps,1-eps,nlag);
  arma::vec weights = pow(1.0 - seq,param2 - 1.0) % pow(seq,param1 - 1.0);
  double sumWeights = arma::accu(weights);
  return  weights/sumWeights; // Return a matrix of nlag x 1
}


//' @export
// [[Rcpp::export]]
double objFunC(Rcpp::NumericVector params,Rcpp::NumericVector yr, 
                Rcpp::NumericMatrix Xr, double q){
  double intercept = params[0];
  double slope = params[1];
  double k1 = params[2];
  double k2 = params[3];
  int nlag = Xr.ncol(),T = Xr.nrow();
  arma::mat X(Xr.begin(), T, nlag, false);
  arma::colvec y(yr.begin(), yr.size(), false);
  arma::mat weights = MidasBetaC(nlag,k1,k2);
  arma::colvec condQuantile = intercept + slope * (X*weights);
  arma::colvec loss = y - condQuantile;
  arma::colvec tickloss = loss % (q - (loss <= 0));
  double fval = arma::accu(tickloss);
  return fval;
}

//' @export
// [[Rcpp::export]]
arma::colvec condQuantileC(Rcpp::NumericVector params,Rcpp::NumericVector yr, 
               Rcpp::NumericMatrix Xr){
  double intercept = params[0];
  double slope = params[1];
  double k1 = params[2];
  double k2 = params[3];
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::mat weights = MidasBetaC(nlag,k1,k2);
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  arma::colvec condQuantile = intercept + slope * (X*weights);
  return condQuantile;
}

// [[Rcpp::export]]
NumericMatrix GetIniParamsC(Rcpp::NumericVector yr, Rcpp::NumericMatrix Xr, double q, 
                     int numInitialsRand = 10000, int numInitials = 10){
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
    NumericVector xx = NumericVector::create(olsIni[0],olsIni[1],1,nInitalBeta[i]);
    double fval = objFunC(xx,yr,Xr,q);
    NumericVector fvalAll = NumericVector::create(fval,olsIni[0],olsIni[1],1,nInitalBeta[i]);
    InitialParamsVec(i,_) = fvalAll;
  }
  return InitialParamsVec;
}

