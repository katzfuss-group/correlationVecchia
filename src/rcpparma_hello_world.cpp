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

arma::vec fun_cpp(const double & a, const arma::vec & x) { return(a * x); }

// [[Rcpp::export]]
double cov_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) { return(covparms[0] - pow(norm(x1 - x2), 2.0)); }

// [[Rcpp::export]]
double cov_expo_iso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) { return(covparms[0] * exp(- norm(x1 - x2) / covparms[1])); }

// [[Rcpp::export]]
double cov_expo_aniso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) { 
    
    arma::rowvec diff = x1 - x2;
    diff[0] = diff[0] * covparms[2];
    
    return( covparms[0] * exp(- norm(diff) / covparms[1]) ); 
    
    // arma::rowvec diagvector(size(x1)); diagvector.ones(); diagvector[0] = pow(covparms[2], 2);
    // return(covparms[0] * exp( - sqrt(dot(x1-x2, (x1-x2) * diagmat(diagvector))) / covparms[1] )); 
}

// [[Rcpp::export]]
double cov_matern_2p5_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) {
    
    double dist = norm(x1 - x2) / covparms[1];
    
    return( covparms[0] * (1 + sqrt(5) * dist + (5.0/3.0) * pow(dist, 2)) * exp(- sqrt(5) * dist) );
}

funcXptr putFunPtrInXPtr(std::string fstr) { 
    
     if (fstr == "fun_cpp") {
         
     return(funcXptr(new funcPtr(&fun_cpp)));
         
     } else {
         
     return(funcXptr(R_NilValue)); // runtime error as NULL no XPtr
         
     }
}

covXptr putCovPtrInXptr(std::string fstr) {
    
    if (fstr == "cov" || fstr == "cov_cpp") {
        
        return(covXptr(new covPtr(&cov_cpp)));
        
    } else if (fstr == "cov_expo_iso" || fstr == "cov_expo_iso_cpp") {
        
        return(covXptr(new covPtr(&cov_expo_iso_cpp)));
        
    } else if (fstr == "cov_expo_aniso" || fstr == "cov_expo_aniso_cpp") {
        
        return(covXptr(new covPtr(&cov_expo_aniso_cpp)));
        
    } else if (fstr == "cov_matern_2.5" || fstr == "cov_matern_2.5_cpp") {
        
        return(covXptr(new covPtr(&cov_matern_2p5_cpp)));
        
    } else {
        
        return(covXptr(R_NilValue));
        
    }
}

// [[Rcpp::export]]
Rcpp::List sortSparse_Rcpp(const arma::mat & x, const double & rho, const int & initInd, std::string distype, std::string fstr, const arma::rowvec & covparms) {
    
    int n = x.n_rows;
    
    covXptr ptr = putCovPtrInXptr(fstr);
    covPtr cov = *ptr;
    
    // Rcpp::XPtr<funcPtr> xpfun(xpsexp);
    // funcPtr dist2Func = *xpfun;
    
    function<double(int, int)> dist2Func = [&, distype, covparms](int i, int j) { 
        if (distype == "euclidean") {
            return(pow(norm(x.row(i) - x.row(j)), 2.0));
        } else if (distype == "correlation") {
            return(covparms[0] - pow(cov(x.row(i), x.row(j), covparms), 2));
        } else {
            return(pow(norm(x.row(i) - x.row(j)), 2.0));
        }
    };
    // function<double(int, int)> dist2Func = [&, covparms](int i, int j) { return(covparms[0] - pow(cov(x.row(i), x.row(j), covparms), 2)); }; // auto dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); };
    // function<double(int, int)> dist2Func = [&, covparms](int i, int j) { return(covparms[0] - abs(cov(x.row(i), x.row(j), covparms))); }; // auto dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); };
    
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
arma::rowvec NNcheck_Rcpp(const arma::rowvec & I, const arma::rowvec & J, const arma::rowvec & P, const arma::rowvec & distances, const arma::mat & x, const double rho, std::string distype, std::string fstr, const arma::rowvec & covparms) {
    
    arma::rowvec chk = arma::ones<arma::rowvec>(arma::size(I));
    
    covXptr ptr = putCovPtrInXptr(fstr);
    covPtr cov = *ptr;
    
    function<double(int, int)> dist2Func = [&, distype, covparms](int i, int j) { 
        if (distype == "euclidean") {
            return(pow(norm(x.row(i) - x.row(j)), 2.0));
        } else if (distype == "correlation") {
            return(covparms[0] - pow(cov(x.row(i), x.row(j), covparms), 2));
        } else {
            return(pow(norm(x.row(i) - x.row(j)), 2.0));
        }
    };
    // function<double(int, int)> dist2Func = [&, covparms](int i, int j) { return(covparms[0] - pow(cov(x.row(i), x.row(j), covparms), 2)); }; // auto dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); }; 
    // function<double(int, int)> dist2Func = [&, covparms](int i, int j) { return(covparms[0] - abs(cov(x.row(i), x.row(j), covparms))); }; // auto dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); };
    
    for (int k = 0; k < I.n_cols; k++) {
        if ( sqrt(dist2Func(P[I[k]], P[J[k]])) > rho * std::min(distances[I[k]], distances[J[k]]) ) {
            chk[k] = 0;
        }
    }
    
    return(chk);
}

// [[Rcpp::export]]
arma::mat conditioning_Rcpp(const arma::rowvec & indvec, const arma::rowvec & condvec, const arma::rowvec & P, const int & maxsize, const arma::mat & x, std::string distype, std::string fstr, const arma::rowvec & covparms) {
    
    int n = x.n_rows;
    int len = indvec.n_elem;
    int k = n - 1;
    
    covXptr ptr = putCovPtrInXptr(fstr);
    covPtr cov = *ptr;
    
    function<double(int, int)> dist2Func = [&, distype, covparms](int i, int j) { 
        if (distype == "euclidean") {
            return(pow(norm(x.row(i) - x.row(j)), 2.0));
        } else if (distype == "correlation") {
            return(covparms[0] - pow(cov(x.row(i), x.row(j), covparms), 2));
        } else {
            return(pow(norm(x.row(i) - x.row(j)), 2.0));
        }
    };
    // function<double(int, int)> dist2Func = [&, covparms](int i, int j) { return(covparms[0] - pow(cov(x.row(i), x.row(j), covparms), 2)); }; // auto dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); };
    // function<double(int, int)> dist2Func = [&, covparms](int i, int j) { return(covparms[0] - abs(cov(x.row(i), x.row(j), covparms))); }; // auto dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); };
    
    arma::mat condsets(n, maxsize); condsets.fill(NA_REAL);
    arma::rowvec freq(n); freq.zeros();
    
    for (int i = 0; i < len; i++) {
        if(indvec[i] == k) {
            condsets.at(k, freq[k]) = condvec.at(i);
            freq[k]++;
        } else {
            k--;
            condsets.at(k, freq[k]) = condvec.at(i);
            freq[k]++;
        }
    }
    
    arma::rowvec distvec(n);
    arma::uvec rankvec(n);
    
    for (int i = 0; i < n; i++) {
        
        distvec.resize(freq[i]); distvec.fill(std::numeric_limits<double>::max());
        for (int j = 0; j < freq[i]; j++) {
            distvec[j] = dist2Func(P.at(i), P.at(condsets.at(i, j)));
        }
        
        rankvec.resize(freq[i]); rankvec = arma::sort_index(distvec);
        
        for (int j = 0; j < freq[i]; j++) {
            distvec.at(j) = condsets.at(i, j);
        }
        
        for (int j = 0; j < freq[i]; j++) {
            condsets.at(i, j) = distvec.at(rankvec.at(j));
        } 
        
    }
    
    return(condsets);
}







