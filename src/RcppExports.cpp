// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cov_cpp
double cov_cpp(const arma::rowvec& x1, const arma::rowvec& x2, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_cov_cpp(SEXP x1SEXP, SEXP x2SEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_cpp(x1, x2, covparms));
    return rcpp_result_gen;
END_RCPP
}
// cov_expo_iso_cpp
double cov_expo_iso_cpp(const arma::rowvec& x1, const arma::rowvec& x2, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_cov_expo_iso_cpp(SEXP x1SEXP, SEXP x2SEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_expo_iso_cpp(x1, x2, covparms));
    return rcpp_result_gen;
END_RCPP
}
// cov_expo_aniso_cpp
double cov_expo_aniso_cpp(const arma::rowvec& x1, const arma::rowvec& x2, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_cov_expo_aniso_cpp(SEXP x1SEXP, SEXP x2SEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_expo_aniso_cpp(x1, x2, covparms));
    return rcpp_result_gen;
END_RCPP
}
// sortSparse_Rcpp
Rcpp::List sortSparse_Rcpp(const arma::mat& x, const double& rho, const int& initInd, std::string fstr, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_sortSparse_Rcpp(SEXP xSEXP, SEXP rhoSEXP, SEXP initIndSEXP, SEXP fstrSEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int& >::type initInd(initIndSEXP);
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(sortSparse_Rcpp(x, rho, initInd, fstr, covparms));
    return rcpp_result_gen;
END_RCPP
}
// NNcheck_Rcpp
arma::rowvec NNcheck_Rcpp(const arma::rowvec& I, const arma::rowvec& J, const arma::rowvec& P, const arma::rowvec& distances, const arma::mat& x, const double rho, std::string fstr, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_NNcheck_Rcpp(SEXP ISEXP, SEXP JSEXP, SEXP PSEXP, SEXP distancesSEXP, SEXP xSEXP, SEXP rhoSEXP, SEXP fstrSEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type distances(distancesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(NNcheck_Rcpp(I, J, P, distances, x, rho, fstr, covparms));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_correlationVecchia_cov_cpp", (DL_FUNC) &_correlationVecchia_cov_cpp, 3},
    {"_correlationVecchia_cov_expo_iso_cpp", (DL_FUNC) &_correlationVecchia_cov_expo_iso_cpp, 3},
    {"_correlationVecchia_cov_expo_aniso_cpp", (DL_FUNC) &_correlationVecchia_cov_expo_aniso_cpp, 3},
    {"_correlationVecchia_sortSparse_Rcpp", (DL_FUNC) &_correlationVecchia_sortSparse_Rcpp, 5},
    {"_correlationVecchia_NNcheck_Rcpp", (DL_FUNC) &_correlationVecchia_NNcheck_Rcpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_correlationVecchia(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
