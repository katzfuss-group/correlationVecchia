// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' @name fun_Rcpp_example
//'
//' @title Example function
//'
//' @param a First argument: real number
//' @param x Second argument: A vector
//'
//' @return Product of a and x
// [[Rcpp::export]]
arma::vec fun_Rcpp_example(const double & a, const arma::vec & x) { return(a * x); }