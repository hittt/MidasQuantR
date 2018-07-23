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

/////////////////////////////////////////////////////
// Functions for the univariate conditional quantiles////
/////////////////////////////////////////////////////

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
NumericMatrix GetIniParamsUni(Rcpp::NumericVector yr, Rcpp::NumericMatrix Xr, 
                            Rcpp::NumericMatrix Xr_neg, Rcpp::NumericMatrix Xr_pos, double q, 
                            int numInitialsRand = 10000, int numInitials = 10, 
                            bool beta2para = false, bool As = false){
  double T = Xr.nrow(), nlag = Xr.ncol();
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::mat Xneg(Xr_neg.begin(),T,nlag,false);
  arma::mat Xpos(Xr_pos.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  arma::colvec nInitalBeta = 99 * arma::randu(numInitialsRand,1) + 1;
  int numPars = 3 + beta2para + As + 1; 
  //  Depends on the model specification, we will have different number of parameters; + 1 is for fval
  // Generate initial guess for the second parameters of the Beta polynomial based on Uniform distribution
  NumericMatrix  InitialParamsVec(numInitialsRand, numPars);
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

////////////////////////////////////////////////////////
// Functions for the joint model with AL distribution///
///////////////////////////////////////////////////////


//' @export
// [[Rcpp::export]]
double objFunAL(Rcpp::NumericVector params, Rcpp::NumericVector yr, Rcpp::NumericVector condmeanR,
                 Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, 
                 Rcpp::NumericMatrix Xr_pos, double q,
                 bool beta2para = false, bool As = false){
  double intercept = params[0];
  double phi = 0.0;
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::colvec condQuantile = as<arma::colvec>(yr);
  arma::mat Xneg(Xr_neg.begin(),T,nlag,false);
  arma::mat Xpos(Xr_pos.begin(),T,nlag,false);
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  arma::colvec mu(condmeanR.begin(),condmeanR.size(),false);
  if(As){
    double slope_neg = params[1];
    double slope_pos = params[2];
    double k1 = 1.0, k2 = 1.0;
    if(beta2para){
      k1 = params[3];
      k2 = params[4];
      phi = params[5];
    } else {
      k2 = params[3];
      phi = params[4];
    }
    arma::mat weights = MidasBetaC(nlag,k1,k2);
    condQuantile = intercept + slope_neg * (Xneg*weights) + slope_pos * (Xpos*weights);  
  } else {
    double slope = params[1];
    double k1 = 1.0, k2 = 1.0;
    if(beta2para){
      k1 = params[2];
      k2 = params[3];
      phi = params[4];
    } else {
      k2 = params[2];
      phi = params[3];
    }
    arma::mat weights = MidasBetaC(nlag,k1,k2);
    condQuantile = intercept + slope * (X*weights);
  }
  arma::colvec es = (1.0 + exp(phi)) * condQuantile;
  arma::colvec loss = y - condQuantile;
  double fval = 0; 
  for(int i = 0; i < T; ++i){
    double hit = q - (loss[i] <= 0);
    double muAdj = mu[i] - es[i];
    double dist = ((1 - q)/muAdj) * exp((-loss[i]*hit)/(q*muAdj));
    if(dist <= 0){
      dist = 1e-100;
    }
    fval += -log(dist);
  }
  return fval;
}

//' @export
// [[Rcpp::export]]
Rcpp::List condVaRES(Rcpp::NumericVector params, Rcpp::NumericVector yr, Rcpp::NumericVector condmeanR,
                     Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, 
                     Rcpp::NumericMatrix Xr_pos,
                     bool beta2para = false, bool As = false){
  double intercept = params[0];
  double phi = 0.0;
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::colvec condQuantile = as<arma::colvec>(yr);
  arma::mat Xneg(Xr_neg.begin(),T,nlag,false);
  arma::mat Xpos(Xr_pos.begin(),T,nlag,false);
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  arma::colvec mu(condmeanR.begin(),condmeanR.size(),false);
  if(As){
    double slope_neg = params[1];
    double slope_pos = params[2];
    double k1 = 1.0, k2 = 1.0;
    if(beta2para){
      k1 = params[3];
      k2 = params[4];
      phi = params[5];
    } else {
      k2 = params[3];
      phi = params[4];
    }
    arma::mat weights = MidasBetaC(nlag,k1,k2);
    condQuantile = intercept + slope_neg * (Xneg*weights) + slope_pos * (Xpos*weights);  
  } else {
    double slope = params[1];
    double k1 = 1.0, k2 = 1.0;
    if(beta2para){
      k1 = params[2];
      k2 = params[3];
      phi = params[4];
    } else {
      k2 = params[2];
      phi = params[3];
    }
    arma::mat weights = MidasBetaC(nlag,k1,k2);
    condQuantile = intercept + slope * (X*weights);
  }
  arma::colvec es = (1.0 + exp(phi)) * condQuantile;
  return Rcpp::List::create(Rcpp::Named("VaR")=condQuantile,
                            Rcpp::Named("ES") =es);
}

//' @export
// [[Rcpp::export]]
NumericMatrix GetIniParamsAL(Rcpp::NumericVector yr, Rcpp::NumericVector condmeanR, Rcpp::NumericVector QuantEst, 
                              Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, Rcpp::NumericMatrix Xr_pos, double q,
                              int numInitialsRand, int numInitials, bool beta2para, bool As){
  double T = Xr.nrow(), nlag = Xr.ncol();
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::mat Xneg(Xr_neg.begin(),T,nlag,false);
  arma::mat Xpos(Xr_pos.begin(),T,nlag,false);
  arma::colvec y(yr.begin(),yr.size(),false);
  arma::colvec mu(condmeanR.begin(),condmeanR.size(),false);
  arma::colvec nInitalPhi = 3 * arma::randu(numInitialsRand,1) - 3;
  int numPars = 4 + beta2para + As + 1; 
  // Depends on the model specification, we will have different number of parameters; + 1 is for fval
  // Generate initial guess for the ES formulation based on Uniform distribution
  NumericMatrix  InitialParamsVec(numInitialsRand,numPars);
  for(int i = 0; i < numInitialsRand;++i){
    NumericVector xx(numPars-1);
    std::copy(QuantEst.begin(),QuantEst.end(),xx.begin());
    xx[numPars - 2] = nInitalPhi[i];
    double fval = objFunAL(xx,yr,condmeanR,Xr,Xr_neg,Xr_pos,q,beta2para,As);
    NumericVector xxFull(numPars);
    xxFull[0] = fval;
    std::copy(xx.begin(),xx.end(),xxFull.begin()+1);
    InitialParamsVec(i,_) = xxFull;
  }
  return InitialParamsVec;
}

//' @export
// [[Rcpp::export]]
