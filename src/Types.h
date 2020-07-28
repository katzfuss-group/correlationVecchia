#pragma once

#include "Rcpp.h"
#include "RcppArmadillo.h"

// using funcPtr = std::function<arma::vec(const arma::vec & x)>;

typedef arma::vec (*funcPtr)(const double & a, const arma::vec & x);
typedef Rcpp::XPtr<funcPtr> funcXptr;

typedef double (*covPtr) (const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms);
typedef Rcpp::XPtr<covPtr> covXptr;
