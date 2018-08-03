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
double objFun_midas(Rcpp::NumericVector params,Rcpp::NumericVector yr, 
                Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, 
                Rcpp::NumericMatrix Xr_pos, double q,
                bool beta2para = false, bool As = false){
  double intercept = params[0];
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::colvec condVaR(T);
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
    condVaR = intercept + slope_neg * (Xneg*weights) + slope_pos * (Xpos*weights);  
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
  condVaR = intercept + slope * (X*weights);
  }
  arma::colvec loss = y - condVaR;
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
arma::colvec condVaR_midas(Rcpp::NumericVector params,
                          Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, 
                          Rcpp::NumericMatrix Xr_pos, bool beta2para = false, bool As = false){
  double intercept = params[0];
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::mat Xneg(Xr_neg.begin(),T,nlag,false);
  arma::mat Xpos(Xr_pos.begin(),T,nlag,false);
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::colvec condVaR(T);
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
    condVaR = intercept + slope_neg * (Xneg*weights) + slope_pos * (Xpos*weights);  
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
    condVaR = intercept + slope * (X*weights);
  }
  return condVaR;
}

// [[Rcpp::export]]
NumericMatrix C_GetIniParams_midas(Rcpp::NumericVector yr, Rcpp::NumericMatrix Xr, 
                            Rcpp::NumericMatrix Xr_neg, Rcpp::NumericMatrix Xr_pos, double q, 
                            int numInitialsRand = 10000,
                            bool beta2para = false, bool As = false){
  int T = Xr.nrow(), nlag = Xr.ncol();
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
        double fval = objFun_midas(xx,yr,Xr,Xr_neg,Xr_pos,q,beta2para,As);
        InitialParamsVec(i,_) = NumericVector::create(fval,olsIni[0],olsIni[1],olsIni[1],1,nInitalBeta[i]);
      } else{
        NumericVector xx = NumericVector::create(olsIni[0],olsIni[1],olsIni[1],nInitalBeta[i]);
        double fval = objFun_midas(xx,yr,Xr,Xr_neg,Xr_pos,q,beta2para,As);
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
        double fval = objFun_midas(xx,yr,Xr,Xr_neg,Xr_pos,q,beta2para,As);
        InitialParamsVec(i,_) = NumericVector::create(fval,olsIni[0],olsIni[1],1,nInitalBeta[i]);
      } else{
        NumericVector xx = NumericVector::create(olsIni[0],olsIni[1],nInitalBeta[i]);
        double fval = objFun_midas(xx,yr,Xr,Xr_neg,Xr_pos,q,beta2para,As);
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
double objFunAL_midas(Rcpp::NumericVector params, Rcpp::NumericVector yr, Rcpp::NumericVector condmeanR,
                 Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, 
                 Rcpp::NumericMatrix Xr_pos, double q,
                 bool beta2para = false, bool As = false){
  double intercept = params[0];
  double phi = 0.0;
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::colvec condVaR(T);
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
    condVaR = intercept + slope_neg * (Xneg*weights) + slope_pos * (Xpos*weights);  
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
    condVaR = intercept + slope * (X*weights);
  }
  arma::colvec condES = (1.0 + exp(phi)) * condVaR;
  arma::colvec loss = y - condVaR;
  double fval = 0; 
  for(int i = 0; i < T; ++i){
    double hit = q - (loss[i] <= 0);
    double muAdj = mu[i] - condES[i];
    double dist = ((1 - q)/muAdj) * exp((-loss[i]*hit)/(q*muAdj));
    if(dist < 0 || Rcpp::traits::is_infinite<REALSXP>(dist)){
      dist = 1e-100;
    }
    fval += -log(dist);
  }
  return fval;
}

//' @export
// [[Rcpp::export]]
Rcpp::List condVaRES_midas(Rcpp::NumericVector params,
                     Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, 
                     Rcpp::NumericMatrix Xr_pos,
                     bool beta2para = false, bool As = false){
  double intercept = params[0];
  double phi = 0.0;
  int T = Xr.nrow(), nlag = Xr.ncol();
  arma::mat Xneg(Xr_neg.begin(),T,nlag,false);
  arma::mat Xpos(Xr_pos.begin(),T,nlag,false);
  arma::mat X(Xr.begin(),T,nlag,false);
  arma::colvec condVaR(T);
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
    condVaR = intercept + slope_neg * (Xneg*weights) + slope_pos * (Xpos*weights);  
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
    condVaR = intercept + slope * (X*weights);
  }
  arma::colvec condES = (1.0 + exp(phi)) * condVaR;
  return Rcpp::List::create(Rcpp::Named("VaR")=condVaR,
                            Rcpp::Named("ES") =condES);
}

// [[Rcpp::export]]
NumericMatrix C_GetIniParamsAL_midas(Rcpp::NumericVector yr, Rcpp::NumericVector condmeanR, Rcpp::NumericVector QuantEst, 
                              Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Xr_neg, Rcpp::NumericMatrix Xr_pos, double q,
                              int numInitialsRand, bool beta2para, bool As){
  arma::colvec nInitalPhi = 3 * arma::randu(numInitialsRand) - 3;
  int numPars = 4 + beta2para + As + 1; 
  // Depends on the model specification, we will have different number of parameters; + 1 is for fval
  // Generate initial guess for the ES formulation based on Uniform distribution
  NumericMatrix  InitialParamsVec(numInitialsRand,numPars);
  for(int i = 0; i < numInitialsRand;++i){
    NumericVector xx(numPars-1);
    std::copy(QuantEst.begin(),QuantEst.end(),xx.begin());
    xx[numPars - 2] = nInitalPhi[i];
    double fval = objFunAL_midas(xx,yr,condmeanR,Xr,Xr_neg,Xr_pos,q,beta2para,As);
    NumericVector xxFull(numPars);
    xxFull[0] = fval;
    std::copy(xx.begin(),xx.end(),xxFull.begin()+1);
    InitialParamsVec(i,_) = xxFull;
  }
  return InitialParamsVec;
}

