// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "rcpplist.h"
// using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

using covPtr = std::function<double(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms)>;

typedef Rcpp::XPtr<covPtr> covXptr;

covXptr putCovPtrInXptr(std::string fstr) {

  if (fstr == "cov_expo_iso" || fstr == "cov_expo_iso_cpp") {

    return(covXptr(new covPtr(&cov_expo_iso_cpp)));

  } else if (fstr == "cov_expo_aniso" || fstr == "cov_expo_aniso_cpp") {

    return(covXptr(new covPtr(&cov_expo_aniso_cpp)));

  } else if (fstr == "cov_expo_spacetime" || fstr == "cov_expo_spacetime_cpp") {

    return(covXptr(new covPtr(&cov_expo_spacetime_cpp)));

  } else if (fstr == "cov_expo_squared" || fstr == "cov_expo_squared_cpp") {

    return(covXptr(new covPtr(&cov_expo_squared_cpp)));

  } else if (fstr == "cov_matern_iso" || fstr == "cov_matern_iso_cpp") {

    return(covXptr(new covPtr(&cov_matern_iso_cpp)));

  } else if (fstr == "cov_matern_aniso" || fstr == "cov_matern_aniso_cpp") {

    return(covXptr(new covPtr(&cov_matern_aniso_cpp)));

  } else if (fstr == "cov_matern_scaledim" || fstr == "cov_matern_scaledim_cpp") {

    return(covXptr(new covPtr(&cov_matern_scaledim_cpp)));

  } else if (fstr == "cov_matern_spacetime" || fstr == "cov_matern_spacetime_cpp") {

    return(covXptr(new covPtr(&cov_matern_spacetime_cpp)));

  } else if (fstr == "cov_latentDim_biv" || fstr == "cov_latentDim_biv_cpp") {

    return(covXptr(new covPtr(&cov_latentDim_biv_cpp)));

  } else if (fstr == "cov_latentDim_triv" || fstr == "cov_latentDim_triv_cpp") {

    return(covXptr(new covPtr(&cov_latentDim_triv_cpp)));

  } else if (fstr == "GpGp_matern_spacetime" || fstr == "GpGp_matern_spacetime_cpp") {

    return(covXptr(new covPtr(&GpGp_matern_spacetime_cpp)));

  } else {

    return(covXptr(R_NilValue));
  }
}

//' @name sortSparse_Rcpp
//'
//' @title Florian's Algorithm without unmeasured inputs (C++ version)
//'
//' @param x A numerical matrix of inputs (locations)
//' @param rho hyperparameter rho (instead of m)
//' @param initInd initial input number
//' @param distype euclidean or correlation
//' @param fstr covariance function name
//' @param covparms A numerical vector with covariance parameters
//'
//' @return List. Check it out!
// [[Rcpp::export]]
Rcpp::List sortSparse_Rcpp(const arma::mat & x, const double & rho, const int & initInd, std::string distype, std::string fstr, const arma::rowvec & covparms) {

  int n = x.n_rows;

  covXptr ptr = putCovPtrInXptr(fstr);
  covPtr cov = *ptr;

  function<double(int, int)> dist2Func = [&, distype, covparms](int i, int j) {

    if (distype == "euclidean") {

      return(pow(norm(x.row(i) - x.row(j)), 2.0));

    } else if (distype == "correlation") {

      return(1 - cov(x.row(i), x.row(j), covparms) / sqrt( cov(x.row(i), x.row(i), covparms) ) / sqrt( cov(x.row(j), x.row(j), covparms) ) );

    } else {

      return(pow(norm(x.row(i) - x.row(j)), 2.0));
    }
  };

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

//' @name predSortSparse_Rcpp
//'
//' @title Florian's Algorithm with unmeasured inputs (C++ version)
//'
//' @param xTrain A numerical matrix of inputs (measured)
//' @param xTest A numerical matrix of inputs (not measured)
//' @param rho hyperparameter rho (instead of m)
//' @param initInd initial input number
//' @param distype euclidean or correlation
//' @param fstr covariance function name
//' @param covparms A numerical vector with covariance parameters
//'
//' @return List. Check it out!
// [[Rcpp::export]]
Rcpp::List predSortSparse_Rcpp(const arma::mat & xTrain, const arma::mat & xTest, const double & rho, const int & initInd, std::string distype, std::string fstr, const arma::rowvec & covparms) {

  int nTrain = xTrain.n_rows;
  int nTest = xTest.n_rows;
  int n = nTrain + nTest;

  arma::mat x = join_vert(xTrain, xTest);

  covXptr ptr = putCovPtrInXptr(fstr);
  covPtr cov = *ptr;

  // Rcpp::XPtr<funcPtr> xpfun(xpsexp);
  // funcPtr dist2Func = *xpfun;

  function<double(int, int)> dist2Func = [&, distype, covparms](int i, int j) {

    if (distype == "euclidean") {

      return(pow(norm(x.row(i) - x.row(j)), 2.0));

    } else if (distype == "correlation") {

      return(1 - cov(x.row(i), x.row(j), covparms) / sqrt( cov(x.row(i), x.row(i), covparms) ) / sqrt( cov(x.row(j), x.row(j), covparms) ) );

    } else {

      return(pow(norm(x.row(i) - x.row(j)), 2.0));
    }
  };
  // function<double(int, int)> dist2Func = [&, covparms](int i, int j) { return(covparms[0] - pow(cov(x.row(i), x.row(j), covparms), 2)); }; // auto dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); };
  // function<double(int, int)> dist2Func = [&, covparms](int i, int j) { return(covparms[0] - abs(cov(x.row(i), x.row(j), covparms))); }; // auto dist2Func = [&](int i, int j) { return(cov(x.row(i), x.row(j))); };

  output result = predSortSparse(nTrain, nTest, rho, dist2Func, initInd);

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

//' @name NNcheck_Rcpp
//'
//' @title NNcheck_Rcpp (only for Florian's algorithm)
//'
//' @param I I
//' @param J J
//' @param P P
//' @param distances distances
//' @param x x
//' @param rho rho
//' @param distype distype
//' @param fstr fstr
//' @param covparms covparms
//'
//' @return Matrix. Check it out!
// [[Rcpp::export]]
arma::rowvec NNcheck_Rcpp(const arma::rowvec & I, const arma::rowvec & J, const arma::rowvec & P, const arma::rowvec & distances, const arma::mat & x, const double rho, std::string distype, std::string fstr, const arma::rowvec & covparms) {

  arma::rowvec chk = arma::ones<arma::rowvec>(arma::size(I));

  covXptr ptr = putCovPtrInXptr(fstr);
  covPtr cov = *ptr;

  function<double(int, int)> dist2Func = [&, distype, covparms](int i, int j) {

    if (distype == "euclidean") {

      return(pow(norm(x.row(i) - x.row(j)), 2.0));

    } else if (distype == "correlation") {

      return(1 - cov(x.row(i), x.row(j), covparms) / sqrt( cov(x.row(i), x.row(i), covparms) ) / sqrt( cov(x.row(j), x.row(j), covparms) ) );

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

//' @name NNcheck_both_Rcpp
//'
//' @title creating both smaller and larger sparsity patterns (only for Florian's algorithm)
//'
//' @param I I
//' @param J J
//' @param P P
//' @param distances distances
//' @param x x
//' @param rho rho
//' @param distype distype
//' @param fstr fstr
//' @param covparms covparms
//'
//' @return list of rowvecs
// [[Rcpp::export]]
Rcpp::List NNcheck_both_Rcpp(const arma::rowvec & I, const arma::rowvec & J, const arma::rowvec & P, const arma::rowvec & distances, const arma::mat & x, const double rho, std::string distype, std::string fstr, const arma::rowvec & covparms) {

  arma::rowvec chkSmaller = arma::ones<arma::rowvec>(arma::size(I));
  arma::rowvec chkLarger = arma::ones<arma::rowvec>(arma::size(I));

  covXptr ptr = putCovPtrInXptr(fstr);
  covPtr cov = *ptr;

  function<double(int, int)> dist2Func = [&, distype, covparms](int i, int j) {

    if (distype == "euclidean") {

      return(pow(norm(x.row(i) - x.row(j)), 2.0));

    } else if (distype == "correlation") {

      return(1 - cov(x.row(i), x.row(j), covparms) / sqrt( cov(x.row(i), x.row(i), covparms) ) / sqrt( cov(x.row(j), x.row(j), covparms) ) );

    } else {

      return(pow(norm(x.row(i) - x.row(j)), 2.0));
    }
  };

  for (int k = 0; k < I.n_cols; k++) {

    // smaller sparsity pattern
    if ( sqrt(dist2Func(P[I[k]], P[J[k]])) > rho * std::min(distances[I[k]], distances[J[k]]) ) {
      chkSmaller[k] = 0;
    }

    // larger sparsity pattern
    if ( sqrt(dist2Func(P[I[k]], P[J[k]])) > rho * std::max(distances[I[k]], distances[J[k]]) ) {
      chkLarger[k] = 0;
    }
  }

  return Rcpp::List::create(Rcpp::Named("smaller") = chkSmaller, Rcpp::Named("larger") = chkLarger);
}

//' @name conditioning_rho_Rcpp
//'
//' @title conditioning based on Florian's algorithm (C++ version)
//'
//' @param indvec indvec
//' @param condvec condvec
//' @param P P
//' @param maxsize maxsize
//' @param x x
//' @param distype euclidean or correlation
//' @param fstr covariance function name
//' @param covparms A numerical vector with covariance parameters
//'
//' @return Matrix. Check it out!
// [[Rcpp::export]]
arma::mat conditioning_rho_Rcpp(const arma::rowvec & indvec, const arma::rowvec & condvec, const arma::rowvec & P, const int & maxsize, const arma::mat & x, std::string distype, std::string fstr, const arma::rowvec & covparms) {

  int n = x.n_rows;
  int len = indvec.n_elem;
  int k = n - 1;

  covXptr ptr = putCovPtrInXptr(fstr);
  covPtr cov = *ptr;

  function<double(int, int)> dist2Func = [&, distype, covparms](int i, int j) {

    if (distype == "euclidean") {

      return(pow(norm(x.row(i) - x.row(j)), 2.0));

    } else if (distype == "correlation") {

      return(1 - cov(x.row(i), x.row(j), covparms) / sqrt( cov(x.row(i), x.row(i), covparms) ) / sqrt( cov(x.row(j), x.row(j), covparms) ) );

    } else {

      return(pow(norm(x.row(i) - x.row(j)), 2.0));
    }
  };

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

int index_smallest(const arma::rowvec target) {

  int ind = target.n_elem - 1;

  for(int i = ind; i >= 0; i--) {

    if(target(i) < target(ind)) {
      ind = i;
    }
  }

  return(ind);
}

arma::irowvec index_smallestk(const arma::rowvec target, const int k) {

  int len = target.n_elem;
  arma::rowvec clone = target;
  arma::irowvec ord(k); ord.fill(NA_INTEGER);

  if (len == 1) {

    ord(0) = 0;
    return(ord);

  } else if (len == 2) {

    if(k == 1) {

      if(target(0) < target(1)) {

        ord(0) = 0;

      } else {

        ord(0) = 1;
      }

      return(ord);

    } else {

      if(target(0) < target(1)) {

        ord(0) = 0; ord(1) = 1;

      } else {

        ord(0) = 1; ord(1) = 0;
      }

      return(ord);
    }

  } else {

    if(k == 1) {

      ord(0) = index_smallest(clone);

    } else {

      for(int i = 0; i < std::min(len, k); i++) {

        ord(i) = index_smallest(clone);
        clone(ord(i)) = std::numeric_limits<double>::max();
      }
    }

    return(ord);
  }
}

//' @name conditioning_m_Rcpp
//'
//' @title Nearest-Neighbor (NN) conditioning with respect to user-specified distance matrix
//'
//' @param m Number of nearby points to condition on
//' @param d A matrix of distances between locations (for instance, 1 - cor)
//'
//' @return A matrix of indices giving NN conditioning sets
// [[Rcpp::export]]
arma::imat conditioning_m_Rcpp(const int m, const arma::mat d) {

  int n = d.n_rows;
  arma::imat NN(n, m + 1);

  arma::rowvec target(n);

  for(int i = 0; i < n; i++) {

    target = d.row(i).head(i+1);
    NN.row(i) = index_smallestk(target, m + 1);
  }

  return( NN );
}
