// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "rcpplist.h"
// using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

//' @name mat2vals
//'
//' @title mat2vals
//'
//' @param ptrs Pointers for sparse matrix
//' @param inds Indices for sparse matrix
//' @param covmat A covariance matrix
//'
//' @return A vector of covariance values
// [[Rcpp::export]]
Rcpp::NumericVector mat2vals(const Rcpp::NumericVector ptrs, const Rcpp::NumericVector inds, const Rcpp::NumericMatrix covmat) {

  const int nvals = inds.size();
  const int N = ptrs.size()-1;
  Rcpp::NumericVector vals(nvals);

  for(int i = 0; i < N; ++i){
    for(int j = ptrs[i]; j < ptrs[i+1]; ++j){

      vals[j] = covmat(i, inds[j]);
    }
  }

  return vals;
}

