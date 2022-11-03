// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "rcpplist.h"
// using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

//' @name conditioning_b1_sptm_Rcpp
//'
//' @title Nearest-Neighbor (NN) conditioning for baseline_1_spacetime_specify()
//'
//' @param n Number of locations
//' @param m Number of nearby points to condition on
//'
//' @return A matrix of indices giving NN conditioning sets
// [[Rcpp::export]]
Rcpp::NumericMatrix conditioning_b1_sptm_Rcpp(const unsigned int n, const unsigned int m)
{
  Rcpp::NumericMatrix condsets(n, m + 1);
  Rcpp::NumericVector cond(m + 1);

  for(unsigned int i = 0; i < n; i++) {

    if(i > m) {

      cond = Rcpp::rev(Rcpp::seq(i-m, i));

    } else {

      cond.fill(NA_REAL);

      for(unsigned int j = 0; j <= i; j++) {

        cond[j] = i-j;
      }
    }

    condsets.row(i) = cond;
  }

  return condsets;
}
