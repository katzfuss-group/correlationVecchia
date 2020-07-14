// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "vector"
#include "math.h"
#include "Types.h"
#include "SortSparse.h"
// using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec fun_cpp(const double & a, const arma::vec & x) { return(a * x); }

// [[Rcpp::export]]
double cov_cpp(const arma::rowvec & x1, const arma::rowvec & x2) { return(pow(norm(x1 - x2), 2.0)); }

funcXptr putFunPtrInXPtr(std::string fstr) { 
    
     if (fstr == "fun_cpp") {
         
     return(funcXptr(new funcPtr(&fun_cpp)));
         
     } else {
         
     return(funcXptr(R_NilValue)); // runtime error as NULL no XPtr
         
     }
}

covXptr putCovPtrInXptr(std::string fstr) {
    
    if (fstr == "cov_cpp") {
        
        return(covXptr(new covPtr(&cov_cpp)));
        
    } else {
        
        return(covXptr(R_NilValue));
        
    }
}

// [[Rcpp::export]]
Rcpp::List sortSparse_Rcpp(const arma::mat & x, const double & rho, const int & initInd, std::string fstr) {
    
    int n = x.n_rows;
    
    covXptr ptr = putCovPtrInXptr(fstr);
    covPtr cov = *ptr;
    
    // Rcpp::XPtr<funcPtr> xpfun(xpsexp);
    // funcPtr dist2Func = *xpfun;
    
    function<double(int, int)> dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); }; // auto dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); };
    
    output result = sortSparse(n, rho, dist2Func, initInd);
    
    vector<signed int> rowvalOut(result.rowval.size(), -1);
    for (int i = 0; i < result.rowval.size(); i++) {
        rowvalOut[i] = result.rowval[i].id;
    }
    
    return Rcpp::List::create(Rcpp::Named("P") = result.P,
                              Rcpp::Named("revP") = result.revP,
                              Rcpp::Named("colptr") = result.colptr,
                              Rcpp::Named("rowval") = rowvalOut,
                              Rcpp::Named("distances") = result.distances);
}

// [[Rcpp::export]]
arma::rowvec NNcheck_Rcpp(const arma::rowvec & I, const arma::rowvec & J, const arma::rowvec & P, const arma::rowvec & distances, const arma::mat & x, const double rho, std::string fstr) {
    
    arma::rowvec chk = arma::ones<arma::rowvec>(arma::size(I));
    
    covXptr ptr = putCovPtrInXptr(fstr);
    covPtr cov = *ptr;
    
    function<double(int, int)> dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); }; 
    
    for (int k = 0; k < I.n_cols; k++) {
        if ( sqrt(dist2Func(P[I[k]], P[J[k]])) > rho * std::min(distances[I[k]], distances[J[k]]) ) {
            chk[k] = 0;
        }
    }
    
    return(chk);
}

