#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::vec MidasBetaC(int nlag, double param1, double param2){
  double eps = 2.2204*pow(10.0,-16.0);
  arma::vec seq = arma::linspace(eps,1-eps,nlag);
  arma::vec weights = pow(1.0 - seq,param2 - 1.0) % pow(seq,param1 - 1.0);
  double sumWeights = arma::accu(weights);
  for(int i = 0; i < nlag; ++i){
    weights[i] = weights[i]/sumWeights;
  }
  return weights; // Return a matrix of nlag x 1
}

//' @export
arma::vec condQuantileC(arma::vec params, arma::vec y,
                        arma::mat X, bool beta2para = false){
  double intercept = params[0];
  double slope = params[1];
  int nrow = X.n_rows, ncol = X.n_cols;
  if (beta2pars){
  double k1 = params[2];
  double k2 = params[3];
  arma::vec weights = MidasBetaC(ncol,k1,k2);
  
  
  arma::vec condQuantile;
  if(beta2para){
    double k1 = params[2];
    double k2 = params[3];
    NumericVector weights = MidasBetaC(ncol,k1,k2);
    for(int i = 0; i < nrow; ++i){
      double rowMul = 0;
      for(int j = 0; j < ncol; ++j){
        rowMul += X(i,j) * weights[j];
      }
      condQuantile[i] = intercept + slope * rowMul;
    }
  } else {
    double k1 = 1;
    double k2 = params[2];
    NumericVector weights = MidasBetaC(ncol,k1,k2);
    for(int i = 0; i < nrow; ++i){
      double rowMul = 0;
      for(int j = 0; j < ncol; ++j){
        rowMul += X(i,j) * weights[j];
      }
      condQuantile[i] = intercept + slope * rowMul;
    }
  }
  return condQuantile;
}
