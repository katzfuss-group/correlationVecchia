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
// cov_bivariate_expo_latDim_cpp
double cov_bivariate_expo_latDim_cpp(const arma::rowvec& newx1, const arma::rowvec& newx2, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_cov_bivariate_expo_latDim_cpp(SEXP newx1SEXP, SEXP newx2SEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type newx1(newx1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type newx2(newx2SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_bivariate_expo_latDim_cpp(newx1, newx2, covparms));
    return rcpp_result_gen;
END_RCPP
}
// cov_matern_iso_cpp
double cov_matern_iso_cpp(const arma::rowvec& x1, const arma::rowvec& x2, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_cov_matern_iso_cpp(SEXP x1SEXP, SEXP x2SEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_matern_iso_cpp(x1, x2, covparms));
    return rcpp_result_gen;
END_RCPP
}
// cov_matern_aniso_cpp
double cov_matern_aniso_cpp(const arma::rowvec& x1, const arma::rowvec& x2, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_cov_matern_aniso_cpp(SEXP x1SEXP, SEXP x2SEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_matern_aniso_cpp(x1, x2, covparms));
    return rcpp_result_gen;
END_RCPP
}
// cov_matern_spacetime_cpp
double cov_matern_spacetime_cpp(const arma::rowvec& x1, const arma::rowvec& x2, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_cov_matern_spacetime_cpp(SEXP x1SEXP, SEXP x2SEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_matern_spacetime_cpp(x1, x2, covparms));
    return rcpp_result_gen;
END_RCPP
}
// cov_matern_2p5_cpp
double cov_matern_2p5_cpp(const arma::rowvec& x1, const arma::rowvec& x2, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_cov_matern_2p5_cpp(SEXP x1SEXP, SEXP x2SEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_matern_2p5_cpp(x1, x2, covparms));
    return rcpp_result_gen;
END_RCPP
}
// sortSparse_Rcpp
Rcpp::List sortSparse_Rcpp(const arma::mat& x, const double& rho, const int& initInd, std::string distype, std::string fstr, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_sortSparse_Rcpp(SEXP xSEXP, SEXP rhoSEXP, SEXP initIndSEXP, SEXP distypeSEXP, SEXP fstrSEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int& >::type initInd(initIndSEXP);
    Rcpp::traits::input_parameter< std::string >::type distype(distypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(sortSparse_Rcpp(x, rho, initInd, distype, fstr, covparms));
    return rcpp_result_gen;
END_RCPP
}
// predSortSparse_Rcpp
Rcpp::List predSortSparse_Rcpp(const arma::mat& xTrain, const arma::mat& xTest, const double& rho, const int& initInd, std::string distype, std::string fstr, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_predSortSparse_Rcpp(SEXP xTrainSEXP, SEXP xTestSEXP, SEXP rhoSEXP, SEXP initIndSEXP, SEXP distypeSEXP, SEXP fstrSEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xTrain(xTrainSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xTest(xTestSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int& >::type initInd(initIndSEXP);
    Rcpp::traits::input_parameter< std::string >::type distype(distypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(predSortSparse_Rcpp(xTrain, xTest, rho, initInd, distype, fstr, covparms));
    return rcpp_result_gen;
END_RCPP
}
// NNcheck_Rcpp
arma::rowvec NNcheck_Rcpp(const arma::rowvec& I, const arma::rowvec& J, const arma::rowvec& P, const arma::rowvec& distances, const arma::mat& x, const double rho, std::string distype, std::string fstr, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_NNcheck_Rcpp(SEXP ISEXP, SEXP JSEXP, SEXP PSEXP, SEXP distancesSEXP, SEXP xSEXP, SEXP rhoSEXP, SEXP distypeSEXP, SEXP fstrSEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type distances(distancesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< std::string >::type distype(distypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(NNcheck_Rcpp(I, J, P, distances, x, rho, distype, fstr, covparms));
    return rcpp_result_gen;
END_RCPP
}
// conditioning_Rcpp
arma::mat conditioning_Rcpp(const arma::rowvec& indvec, const arma::rowvec& condvec, const arma::rowvec& P, const int& maxsize, const arma::mat& x, std::string distype, std::string fstr, const arma::rowvec& covparms);
RcppExport SEXP _correlationVecchia_conditioning_Rcpp(SEXP indvecSEXP, SEXP condvecSEXP, SEXP PSEXP, SEXP maxsizeSEXP, SEXP xSEXP, SEXP distypeSEXP, SEXP fstrSEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type indvec(indvecSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type condvec(condvecSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxsize(maxsizeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type distype(distypeSEXP);
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(conditioning_Rcpp(indvec, condvec, P, maxsize, x, distype, fstr, covparms));
    return rcpp_result_gen;
END_RCPP
}
// index_smallest
int index_smallest(const arma::rowvec target);
RcppExport SEXP _correlationVecchia_index_smallest(SEXP targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec >::type target(targetSEXP);
    rcpp_result_gen = Rcpp::wrap(index_smallest(target));
    return rcpp_result_gen;
END_RCPP
}
// index_smallestk
arma::irowvec index_smallestk(const arma::rowvec target, const int k);
RcppExport SEXP _correlationVecchia_index_smallestk(SEXP targetSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(index_smallestk(target, k));
    return rcpp_result_gen;
END_RCPP
}
// conditioning_nn_Rcpp
arma::imat conditioning_nn_Rcpp(const int m, const arma::mat d);
RcppExport SEXP _correlationVecchia_conditioning_nn_Rcpp(SEXP mSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(conditioning_nn_Rcpp(m, d));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_correlationVecchia_cov_cpp", (DL_FUNC) &_correlationVecchia_cov_cpp, 3},
    {"_correlationVecchia_cov_expo_iso_cpp", (DL_FUNC) &_correlationVecchia_cov_expo_iso_cpp, 3},
    {"_correlationVecchia_cov_expo_aniso_cpp", (DL_FUNC) &_correlationVecchia_cov_expo_aniso_cpp, 3},
    {"_correlationVecchia_cov_bivariate_expo_latDim_cpp", (DL_FUNC) &_correlationVecchia_cov_bivariate_expo_latDim_cpp, 3},
    {"_correlationVecchia_cov_matern_iso_cpp", (DL_FUNC) &_correlationVecchia_cov_matern_iso_cpp, 3},
    {"_correlationVecchia_cov_matern_aniso_cpp", (DL_FUNC) &_correlationVecchia_cov_matern_aniso_cpp, 3},
    {"_correlationVecchia_cov_matern_spacetime_cpp", (DL_FUNC) &_correlationVecchia_cov_matern_spacetime_cpp, 3},
    {"_correlationVecchia_cov_matern_2p5_cpp", (DL_FUNC) &_correlationVecchia_cov_matern_2p5_cpp, 3},
    {"_correlationVecchia_sortSparse_Rcpp", (DL_FUNC) &_correlationVecchia_sortSparse_Rcpp, 6},
    {"_correlationVecchia_predSortSparse_Rcpp", (DL_FUNC) &_correlationVecchia_predSortSparse_Rcpp, 7},
    {"_correlationVecchia_NNcheck_Rcpp", (DL_FUNC) &_correlationVecchia_NNcheck_Rcpp, 9},
    {"_correlationVecchia_conditioning_Rcpp", (DL_FUNC) &_correlationVecchia_conditioning_Rcpp, 8},
    {"_correlationVecchia_index_smallest", (DL_FUNC) &_correlationVecchia_index_smallest, 1},
    {"_correlationVecchia_index_smallestk", (DL_FUNC) &_correlationVecchia_index_smallestk, 2},
    {"_correlationVecchia_conditioning_nn_Rcpp", (DL_FUNC) &_correlationVecchia_conditioning_nn_Rcpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_correlationVecchia(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
