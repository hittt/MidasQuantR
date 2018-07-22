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
double objFun(Rcpp::NumericVector params,Rcpp::NumericVector yr, 
                Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, 
                Rcpp::NumericMatrix Xr_pos, double q,
                bool beta2para = false, bool As = false){
  double intercept = params[0];
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::colvec condQuantile = as<arma::colvec>(yr);
  arma::mat Xneg(Xr_neg.begin(),T,nlag,false);
  arma::mat Xpos(Xr_pos.begin(),T,nlag,false);
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  if(As){
    double slope_neg = params[1];
    double slope_pos = params[2];
    double k1 = 1.0, k2 = 1.0;
    if(beta2para){
      k1 = params[3];
      k2 = params[4];  
    } else {
      k2 = params[3];
    }
    arma::mat weights = MidasBetaC(nlag,k1,k2);
    condQuantile = intercept + slope_neg * (Xneg*weights) + slope_pos * (Xpos*weights);  
  } else {
  double slope = params[1];
  double k1 = 1.0, k2 = 1.0;
  if(beta2para){
    k1 = params[2];
    k2 = params[3];  
  } else {
    k2 = params[2];
  }
  arma::mat weights = MidasBetaC(nlag,k1,k2);
  condQuantile = intercept + slope * (X*weights);
  }
  arma::colvec loss = y - condQuantile;
  double fval = 0; 
  for(int i = 0; i < T; ++i){
    fval += loss[i] * (q - (loss[i] <= 0));
  }
  return fval;
}

//' @export
// [[Rcpp::export]]
arma::colvec condQuantile(Rcpp::NumericVector params,Rcpp::NumericVector yr, 
                          Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, 
                          Rcpp::NumericMatrix Xr_pos, bool beta2para = false, bool As = false){
  double intercept = params[0];
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::colvec condQuantile = as<arma::colvec>(yr);
  arma::mat Xneg(Xr_neg.begin(),T,nlag,false);
  arma::mat Xpos(Xr_pos.begin(),T,nlag,false);
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  if(As){
    double slope_neg = params[1];
    double slope_pos = params[2];
    double k1 = 1.0, k2 = 1.0;
    if(beta2para){
      k1 = params[3];
      k2 = params[4];  
    } else {
      k2 = params[3];
    }
    arma::mat weights = MidasBetaC(nlag,k1,k2);
    condQuantile = intercept + slope_neg * (Xneg*weights) + slope_pos * (Xpos*weights);  
  } else {
    double slope = params[1];
    double k1 = 1.0, k2 = 1.0;
    if(beta2para){
      k1 = params[2];
      k2 = params[3];  
    } else {
      k2 = params[2];
    }
    arma::mat weights = MidasBetaC(nlag,k1,k2);
    condQuantile = intercept + slope * (X*weights);
  }
  return condQuantile;
}

//' @export
// [[Rcpp::export]]
NumericMatrix GetIniParamsC(Rcpp::NumericVector yr, Rcpp::NumericMatrix Xr, 
                            Rcpp::NumericMatrix Xr_neg, Rcpp::NumericMatrix Xr_pos, double q, 
                            int numInitialsRand = 10000, int numInitials = 10, 
                            bool beta2para = false, bool As = false){
  double T = Xr.nrow(), nlag = Xr.ncol();
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::mat Xneg(Xr_neg.begin(),T,nlag,false);
  arma::mat Xpos(Xr_pos.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  arma::colvec nInitalBeta = 99 * arma::randu(numInitialsRand,1) + 1;
  // Generate initial guess for the second parameters of the Beta polynomial based on Uniform distribution
  NumericMatrix  InitialParamsVec(numInitialsRand,5);
  if(As){
    for(int i = 0; i < numInitialsRand;++i){
    arma::vec weights = MidasBetaC(nlag, 1 , nInitalBeta[i]);
    arma::mat Xweighted = Xneg*weights + Xpos*weights;
    arma::mat X0 = arma::join_rows(arma::ones(T,1),Xweighted);
    arma::mat olsIni = solve(X0,y);
    if(beta2para){
      NumericVector xx = NumericVector::create(olsIni[0],olsIni[1],olsIni[1],1,nInitalBeta[i]);
      double fval = objFun(xx,yr,Xr,Xr_neg,Xr_pos,q,beta2para,As);
      InitialParamsVec(i,_) = NumericVector::create(fval,olsIni[0],olsIni[1],olsIni[1],1,nInitalBeta[i]);
    } else{
      NumericVector xx = NumericVector::create(olsIni[0],olsIni[1],olsIni[1],nInitalBeta[i]);
      double fval = objFun(xx,yr,Xr,Xr_neg,Xr_pos,q,beta2para,As);
      InitialParamsVec(i,_) = NumericVector::create(fval,olsIni[0],olsIni[1],olsIni[1],nInitalBeta[i]);  
    }
  }
  } else{
    for(int i = 0; i < numInitialsRand;++i){
      arma::vec weights = MidasBetaC(nlag, 1 , nInitalBeta[i]);
      arma::mat Xweighted = X*weights;
      arma::mat X0 = arma::join_rows(arma::ones(T,1),Xweighted);
      arma::mat olsIni = solve(X0,y);
      if(beta2para){
        NumericVector xx = NumericVector::create(olsIni[0],olsIni[1],1,nInitalBeta[i]);
        double fval = objFun(xx,yr,Xr,Xr_neg,Xr_pos,q,beta2para,As);
        InitialParamsVec(i,_) = NumericVector::create(fval,olsIni[0],olsIni[1],1,nInitalBeta[i]);
      } else{
        NumericVector xx = NumericVector::create(olsIni[0],olsIni[1],nInitalBeta[i]);
        double fval = objFun(xx,yr,Xr,Xr_neg,Xr_pos,q,beta2para,As);
        InitialParamsVec(i,_) = NumericVector::create(fval,olsIni[0],olsIni[1],nInitalBeta[i]);  
      }
    }
  }
  return InitialParamsVec;
}