// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// MidasBetaC
arma::vec MidasBetaC(int nlag, double param1, double param2);
RcppExport SEXP _MidasQuantR_MidasBetaC(SEXP nlagSEXP, SEXP param1SEXP, SEXP param2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nlag(nlagSEXP);
    Rcpp::traits::input_parameter< double >::type param1(param1SEXP);
    Rcpp::traits::input_parameter< double >::type param2(param2SEXP);
    rcpp_result_gen = Rcpp::wrap(MidasBetaC(nlag, param1, param2));
    return rcpp_result_gen;
END_RCPP
}
// objFunC
double objFunC(Rcpp::NumericVector params, Rcpp::NumericVector yr, Rcpp::NumericMatrix Xr, double q);
RcppExport SEXP _MidasQuantR_objFunC(SEXP paramsSEXP, SEXP yrSEXP, SEXP XrSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type yr(yrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Xr(XrSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(objFunC(params, yr, Xr, q));
    return rcpp_result_gen;
END_RCPP
}
// condQuantileC
arma::colvec condQuantileC(Rcpp::NumericVector params, Rcpp::NumericVector yr, Rcpp::NumericMatrix Xr);
RcppExport SEXP _MidasQuantR_condQuantileC(SEXP paramsSEXP, SEXP yrSEXP, SEXP XrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type yr(yrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Xr(XrSEXP);
    rcpp_result_gen = Rcpp::wrap(condQuantileC(params, yr, Xr));
    return rcpp_result_gen;
END_RCPP
}
// GetIniParamsC
NumericMatrix GetIniParamsC(Rcpp::NumericVector yr, Rcpp::NumericMatrix Xr, double q, int numInitialsRand, int numInitials);
RcppExport SEXP _MidasQuantR_GetIniParamsC(SEXP yrSEXP, SEXP XrSEXP, SEXP qSEXP, SEXP numInitialsRandSEXP, SEXP numInitialsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type yr(yrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Xr(XrSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type numInitialsRand(numInitialsRandSEXP);
    Rcpp::traits::input_parameter< int >::type numInitials(numInitialsSEXP);
    rcpp_result_gen = Rcpp::wrap(GetIniParamsC(yr, Xr, q, numInitialsRand, numInitials));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MidasQuantR_MidasBetaC", (DL_FUNC) &_MidasQuantR_MidasBetaC, 3},
    {"_MidasQuantR_objFunC", (DL_FUNC) &_MidasQuantR_objFunC, 4},
    {"_MidasQuantR_condQuantileC", (DL_FUNC) &_MidasQuantR_condQuantileC, 3},
    {"_MidasQuantR_GetIniParamsC", (DL_FUNC) &_MidasQuantR_GetIniParamsC, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_MidasQuantR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
