#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// These above headers are neccessary to use RcppArmadillo

//' @export
// [[Rcpp::export]]
double objFun_cav(arma::colvec params, Rcpp::NumericVector yr, Rcpp::NumericVector Xr,
                    double q, int model, double empQuant, bool Uni = true){
  int T = yr.size();
  Rcpp::NumericVector condVaR(T);
  Rcpp::NumericVector absYr = Rcpp::abs(yr);
  condVaR[0] = empQuant;
  if(Uni){
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1 * absYr[i-1] + beta2 * condVaR[i-1] ;
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0)  + beta2 * condVaR[i-1];
      }
    }
  } else{
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Parameters for additional explanatory variable
      double beta3 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1 * absYr[i-1] + beta2 * Xr[i-1] + beta3 * condVaR[i-1];
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Parameters for the additional explanatory variable
      double beta3 = params[4]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0) + beta2 * Xr[i-1] + beta3 * condVaR[i-1];
      }
    }
  }
  Rcpp::NumericVector loss = yr - condVaR;
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
Rcpp::NumericVector condVaR_cav(Rcpp::NumericVector params, Rcpp::NumericVector yr, Rcpp::NumericVector Xr,
                                    int model, double empQuant, bool Uni = true){
  int T = yr.size();
  Rcpp::NumericVector condVaR(T);
  Rcpp::NumericVector absYr = Rcpp::abs(yr);
  condVaR[0] = empQuant;
  if(Uni){
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1 * absYr[i-1] + beta2 * condVaR[i-1] ;
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0)  + beta2 * condVaR[i-1];
      }
    }
  } else{
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Parameters for additional explanatory variable
      double beta3 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1 * absYr[i-1] + beta2 * Xr[i-1] + beta3 * condVaR[i-1];
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Parameters for the additional explanatory variable
      double beta3 = params[4]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0) + beta2 * Xr[i-1] + beta3 * condVaR[i-1];
      }
    }
  }
  return condVaR;
}

// [[Rcpp::export]]
NumericMatrix C_GetIniParams_cav(Rcpp::NumericVector yr, Rcpp::NumericVector Xr, double q, int model, 
                                 double empQuant, bool Uni = true, int numInitialsRand  = 10000){
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
    double fval = objFun_cav(candidatePars,yr,Xr,q,model,empQuant,Uni);
    NumericVector xx(numPars + 1);
    xx[0] = fval;
    std::copy(candidatePars.begin(),candidatePars.end(),xx.begin()+1);
    InitialParamesVec(i,_) = xx;
  }
  return InitialParamesVec;
}

////////////////////////////////////////////////////////
// Functions for the joint model with AL distribution///
///////////////////////////////////////////////////////


//' @export
// [[Rcpp::export]]
double objFunAL_cav(arma::colvec params, Rcpp::NumericVector yr, Rcpp::NumericVector Xr,
                    Rcpp::NumericVector condmeanR,double q, double empQuant,
                      int model = false, bool Uni = false){
  int T = yr.size();
  arma::colvec condVaR(T);
  Rcpp::NumericVector absYr = Rcpp::abs(yr);
  condVaR[0] = empQuant;
  arma::colvec mu(condmeanR.begin(),condmeanR.size(),false);
  arma::colvec y(yr.begin(),yr.size(),false);
  double phi = 0;
  if(Uni){
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Autoregressive parameter
      phi = params[3];
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1 * absYr[i-1] + beta2 * condVaR[i-1] ;
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0)  + beta2 * condVaR[i-1];
      }
    }
  } else{
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Parameters for additional explanatory variable
      double beta3 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1 * absYr[i-1] + beta2 * Xr[i-1] + beta3 * condVaR[i-1];
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Parameters for the additional explanatory variable
      double beta3 = params[4]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0) + beta2 * Xr[i-1] + beta3 * condVaR[i-1];
      }
    }
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
Rcpp::List condVaRES_cav(arma::colvec params, Rcpp::NumericVector yr, Rcpp::NumericVector Xr,
                           Rcpp::NumericVector condmeanR,double empQuant,
                           int model = false, bool Uni = false){
  int T = yr.size();
  arma::colvec condVaR(T);
  Rcpp::NumericVector absYr = Rcpp::abs(yr);
  condVaR[0] = empQuant;
  arma::colvec mu(condmeanR.begin(),condmeanR.size(),false);
  arma::colvec y(yr.begin(),yr.size(),false);
  double phi = 0;
  if(Uni){
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Autoregressive parameter
      phi = params[3];
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1 * absYr[i-1] + beta2 * condVaR[i-1] ;
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0)  + beta2 * condVaR[i-1];
      }
    }
  } else{
    if(model == 1){
      double beta0 = params[0]; // Intercept parameter
      double beta1 = params[1]; // Parameters with the lag absolute returns
      double beta2 = params[2]; // Parameters for additional explanatory variable
      double beta3 = params[3]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1 * absYr[i-1] + beta2 * Xr[i-1] + beta3 * condVaR[i-1];
      }
    }else{
      double beta0 = params[0]; // Intercept parameter
      double beta1_neg = params[1]; //Paramter with the lag negative returns
      double beta1_pos = params[2]; //Parameter with the lag positive returns
      double beta2 = params[3]; // Parameters for the additional explanatory variable
      double beta3 = params[4]; // Autoregressive parameter
      for(int i = 1; i < T; ++i){
        condVaR[i] = beta0 + beta1_neg * absYr[i-1] * (yr[i-1] <= 0) + beta1_pos * absYr[i-1] * (yr[i-1] > 0) + beta2 * Xr[i-1] + beta3 * condVaR[i-1];
      }
    }
  }
  arma::colvec condES = (1.0 + exp(phi)) * condVaR;
  return Rcpp::List::create(Rcpp::Named("VaR")=condVaR,
                            Rcpp::Named("ES") =condES);
}

// [[Rcpp::export]]
NumericMatrix C_GetIniParamsAL_cav(Rcpp::NumericVector yr, Rcpp::NumericVector Xr,
                                   Rcpp::NumericVector condmeanR, Rcpp::NumericVector QuantEst,
                                   double q, double empQuant,int numInitialsRand, int model = false, bool Uni = false){
  arma::colvec nInitalPhi = 3 * arma::randu(numInitialsRand) - 3;
  int numPars = 2 + model + Uni + 1; 
  // Depends on the model specification, we will have different number of parameters; + 1 is for fval
  // Generate initial guess for the ES formulation based on Uniform distribution
  NumericMatrix  InitialParamsVec(numInitialsRand,numPars);
  for(int i = 0; i < numInitialsRand;++i){
    NumericVector xx(numPars-1);
    std::copy(QuantEst.begin(),QuantEst.end(),xx.begin());
    xx[numPars - 2] = nInitalPhi[i];
    double fval = objFunAL_cav(xx,yr,Xr,condmeanR,q,empQuant,model,Uni);
    NumericVector xxFull(numPars);
    xxFull[0] = fval;
    std::copy(xx.begin(),xx.end(),xxFull.begin()+1);
    InitialParamsVec(i,_) = xxFull;
  }
  return InitialParamsVec;
}

// [[Rcpp::export]]
NumericVector cavSim(Rcpp::NumericVector params, Rcpp::NumericVector ResidSim, 
                     double NegResidMean, Rcpp::NumericVector y, Rcpp::NumericVector condVaR, 
                     Rcpp::NumericVector condES, int model, bool Uni, Rcpp::NumericVector Xsim){
  int T = ResidSim.size();
  NumericVector ySim(T);
  NumericVector VaRsim(T);
  NumericVector ESsim(T);
  NumericVector absY = Rcpp::abs(y);
  ySim[0] = y[0];
  VaRsim[0] = condVaR[0]; 
  ESsim[0] = condES[0];
  for(int i = 1; i < T; ++i){
    if(Uni){
      if(model == 1){
        VaRsim[i] = params[0] + params[1] * absY[i-1] + params[2] * VaRsim[i-1];
        ESsim[i] = (1 + exp(params[3])) * VaRsim[i];
        double ySimDay = 0.0;
        if(VaRsim[i] <= 0){
          ySimDay = VaRsim[i] - ResidSim[i] * VaRsim[i];
        } else {
          ySimDay = VaRsim[i] + ResidSim[i] * VaRsim[i];
        }
        if(ySimDay <= VaRsim[i]){
          ySimDay = VaRsim[i] + (ResidSim[i]/NegResidMean) * (ESsim[i] - VaRsim[i]);
        }
        ySim[i] = ySimDay;
      } else{
        VaRsim[i] = params[0] + params[1] * absY[i-1] * (y[i-1] <= 0) + params[2] * absY[i-1] * (y[i-1] > 0) + params[3] * VaRsim[i-1];
        ESsim[i] = (1 + exp(params[4])) * VaRsim[i];
        double ySimDay = 0.0;
        if(VaRsim[i] <= 0){
          ySimDay = VaRsim[i] - ResidSim[i] * VaRsim[i];
        } else {
          ySimDay = VaRsim[i] + ResidSim[i] * VaRsim[i];
        }
        if(ySimDay <= VaRsim[i]){
          ySimDay = VaRsim[i] + (ResidSim[i]/NegResidMean) * (ESsim[i] - VaRsim[i]);
        }
        ySim[i] = ySimDay;
      }
    } else{
      if(model == 1){
        VaRsim[i] = params[0] + params[1] * absY[i-1] + params[2] * Xsim[i-1] + params[3] * VaRsim[i-1];
        ESsim[i] = (1 + exp(params[4])) * VaRsim[i];
        double ySimDay = 0.0;
        if(VaRsim[i] <= 0){
          ySimDay = VaRsim[i] - ResidSim[i] * VaRsim[i];
        } else {
          ySimDay = VaRsim[i] + ResidSim[i] * VaRsim[i];
        }
        if(ySimDay <= VaRsim[i]){
          ySimDay = VaRsim[i] + (ResidSim[i]/NegResidMean) * (ESsim[i] - VaRsim[i]);
        }
        ySim[i] = ySimDay;
      } else{
        VaRsim[i] = params[0] + params[1] * absY[i-1] * (y[i-1] <= 0) + params[2] * absY[i-1] * (y[i-1] > 0) + params[3] * Xsim[i-1] + params[4] * VaRsim[i-1];
        ESsim[i] = (1 + exp(params[5])) * VaRsim[i];
        double ySimDay = 0.0;
        if(VaRsim[i] <= 0){
          ySimDay = VaRsim[i] - ResidSim[i] * VaRsim[i];
        } else {
          ySimDay = VaRsim[i] + ResidSim[i] * VaRsim[i];
        }
        if(ySimDay <= VaRsim[i]){
          ySimDay = VaRsim[i] + (ResidSim[i]/NegResidMean) * (ESsim[i] - VaRsim[i]);
        }
        ySim[i] = ySimDay;
      }
    }
  }
  return ySim;
}