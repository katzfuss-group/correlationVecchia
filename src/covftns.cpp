// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "rcpplist.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

//' @name cov_expo_iso_cpp
//'
//' @title Isotropic exponential covariance function (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters in the form (variance, range)
//'
//' @return The isotropic exponenital covariance between the locations provided
// [[Rcpp::export]]
double cov_expo_iso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) { 
  
  return( covparms(0) * exp(- norm(x1 - x2) / covparms(1)) ); 
}

//' @name cov_expo_aniso_cpp
//'
//' @title Anisotropic exponential covariance function (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters in the form (variance, range, anisotropy for the first dimension)
//'
//' @return The anisotropic exponenital covariance between the locations provided
// [[Rcpp::export]]
double cov_expo_aniso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) { 
  
  arma::rowvec diff = x1 - x2;
  diff(0) = diff(0) * covparms(2);
  
  return( covparms(0) * exp(- norm(diff) / covparms(1)) ); 
}

//' @name cov_expo_spacetime_cpp
//'
//' @title Space-time exponential covariance function (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters in the form (variance, spatial range, temporal range)
//'
//' @return The space-time exponenital covariance between the locations provided
// [[Rcpp::export]]
double cov_expo_spacetime_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) { 

  arma::rowvec diff = x1 - x2;
  
  for(unsigned int i = 0; i < diff.n_elem - 1; i++) {
    diff(i) = diff(i) / covparms(1);
  }

  diff(diff.n_elem - 1) = diff(diff.n_elem - 1) / covparms(2);
  
  return( covparms(0) * exp(- norm(diff) ) ); 
}

//' @name cov_expo_squared_cpp
//'
//' @title Squared exponential covariance function (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters in the form (variance, range)
//'
//' @return The squared exponenital covariance between the locations provided
// [[Rcpp::export]]
double cov_expo_squared_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) {

  return( covparms(0) * exp(- pow(norm(x1 - x2), 2.0) / 2 / pow(covparms(1), 2.0) ) );
}

//' @name cov_matern_iso_cpp
//'
//' @title Isotropic Matern covariance function (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters in the form (sigma2, range, smoothness)
//'
//' @return The isotropic Matern covariance between the locations provided
// [[Rcpp::export]]
double cov_matern_iso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) {
  
  arma::rowvec diff = x1 - x2;
  double d = norm(diff);
  double covconst = covparms(0) / ( pow(2.0, covparms(2) - 1.0) * boost::math::tgamma(covparms(2)) );
  double cov = -1.0;
  
  if(d == 0.0) {
    
    cov = covparms(0);
    return(cov);
    
  } else {
    
    d = pow(2 * covparms(2), 0.5) * d / covparms(1);
    cov = covconst * pow(d, covparms(2)) * boost::math::cyl_bessel_k(covparms(2), d);
    return(cov);
  }
}

//' @name cov_matern_aniso_cpp
//'
//' @title Anisotropic Matern covariance function (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters in the form (sigma2, range, smoothness, anisotropy for the first dimension)
//'
//' @return The Anisotropic Matern covariance between the locations provided
// [[Rcpp::export]]
double cov_matern_aniso_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) {
  
  arma::rowvec diff = x1 - x2;
  diff(0) = diff(0) * covparms(3);
  
  double d = norm(diff);
  double covconst = covparms(0) / ( pow(2.0, covparms(2) - 1.0) * boost::math::tgamma(covparms(2)) );
  double cov = -1.0;
  
  if(d == 0.0) {
    
    cov = covparms(0);
    return(cov);
    
  } else {
    
    d = pow(2 * covparms(2), 0.5) * d / covparms(1);
    cov = covconst * pow(d, covparms(2)) * boost::math::cyl_bessel_k(covparms(2), d);
    return(cov);
  }
}

//' @name cov_matern_scaledim_cpp
//'
//' @title Anisotropic Matern covariance function (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters = (sigma2, range_1, ... , range_d, smoothness)
//'
//' @return The Anisotropic Matern covariance between the locations provided
// [[Rcpp::export]]
double cov_matern_scaledim_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) {
  
  // CAUTION: different argument setting: please see ?GpGp::matern_scaledim()
  
  arma::rowvec diff = x1 - x2;
  double nu = covparms(covparms.n_elem-1);
  
  for(unsigned int i = 0; i < diff.n_elem; i++) {
    diff(i) = diff(i) / covparms(i+1);
  }
  
  double d = norm(diff);
  double covconst = covparms(0) / ( pow(2.0, nu - 1.0) * boost::math::tgamma(nu) );
  double cov = -1.0;
  
  if(d == 0.0) {
    
    cov = covparms(0);
    return(cov);
    
  } else {
    
    d = pow(2 * nu, 0.5) * d;
    cov = covconst * pow(d, nu) * boost::math::cyl_bessel_k(nu, d);
    return(cov);
  }
}

//' @name cov_matern_spacetime_cpp
//'
//' @title Space-time Matern covariance function (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters = (sigma2, spatial range, temporal range, smoothness)
//'
//' @return The space-time Matern covariance between the locations provided
// [[Rcpp::export]]
double cov_matern_spacetime_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec & covparms) {
  
  arma::rowvec diff = x1 - x2;
  unsigned int dim = x1.n_elem - 1;
  arma::rowvec origin(dim + 1); origin.zeros();
  
  for(unsigned int i = 0; i < dim; i++) {
    diff(i) = diff(i) / covparms(1);
  }
  
  diff(dim) = diff(dim) / covparms(2);
  
  arma::rowvec newparms(3);
  newparms(0) = covparms(0); // variance parameter
  newparms(1) = 1.0; // range parameter = 1
  newparms(2) = covparms(3); // smoothness parameter
  
  double cov = cov_matern_iso_cpp(diff, origin, newparms);
  return(cov);
}

double spectralMatern(double d, double nu) {
  
  if(d == 0.0) return( 1.0 );
  
  double k = pow(2.0, 1-nu) / boost::math::tgamma( nu );
  
  k = k * pow(d, nu) * boost::math::cyl_bessel_k(nu, d);
  
  return( k );
}

//' @name cov_matern_ns_forR
//'
//' @title Nonstationary Matern covariance function (similar to R version, but slower than R version)
//'
//' @param locs1 A numerical matrix with \code{n} rows and \code{p} columns. Each row of locs1 gives a point of the first set in R^d
//' @param locs2 A numerical matrix with \code{n} rows and \code{p} columns. Each row of locs2 gives a point of the second set in R^d
//' @param sigma A numeric function for spatially varying standard deviation \eqn{\sigma( loc )}
//' @param smoothness A numeric function for spatially varying smoothness \eqn{\nu( loc )}
//' @param kernel A matrix-valued function for spatially varying (local) geometric anisotropy \eqn{\Sigma( loc )}
//'
//' @return The nonstationary Matern covariance between the locations provided
// [[Rcpp::export]]
arma::mat cov_matern_ns_forR(const arma::mat & locs1, const arma::mat & locs2, Rcpp::Function sigma, Rcpp::Function smoothness, Rcpp::Function kernel) {
  
  // https://stackoverflow.com/questions/17958168/rcpp-function-to-be-slower-than-same-r-function
  // https://stackoverflow.com/questions/27391472/passing-r-function-as-parameter-to-rcpp-function
  
  unsigned int n1 = locs1.n_rows;
  unsigned int n2 = locs2.n_rows;
  unsigned int d = locs1.n_cols;
  
  arma::rowvec tarvec1( d );
  arma::rowvec tarvec2( d );
  
  arma::colvec tarvec3( d );
  
  double sigma12;
  double smooth12;
  arma::mat kernel12( d, d );
  double q12;
  
  arma::mat covmat( n1, n2 );
  for(unsigned int i = 0; i < n1; i++) {
    for(unsigned int j = 0; j < n2; j++) {
      
      tarvec1 = locs1.row(i);
      tarvec2 = locs2.row(j);
      tarvec3 = arma::conv_to< arma::colvec >::from(locs1.row(i) - locs2.row(j)); 
      
      sigma12 = Rcpp::as<double>(sigma(tarvec1)) * Rcpp::as<double>(sigma(tarvec2));
      
      smooth12 = ( Rcpp::as<double>(smoothness(tarvec1)) + Rcpp::as<double>(smoothness(tarvec2)) ) / 2.0;
      
      // Rprintf("sigma[%i, %i] : %f \n", i, j, sigma12);
      // Rprintf("nu[%i, %i] : %f \n", i, j, smooth12);
      
      // kernel12 = ( Rcpp::as<Rcpp::NumericMatrix>(kernel(tarvec1)) + Rcpp::as<Rcpp::NumericMatrix>(kernel(tarvec2)) ) / 2.0;
      kernel12 = ( Rcpp::as<arma::mat>(kernel(tarvec1)) + Rcpp::as<arma::mat>(kernel(tarvec2)) ) / 2.0;
      
      q12 = arma::dot( tarvec3, arma::solve(kernel12, tarvec3) );
      
      // Rprintf("q[%i, %i] : %f \n", i, j, q12);

      // Rprintf("cov_matern_iso_cpp[%i, %i]: %f \n", i, j, cov_matern_iso_cpp(x1 = x1, x2 = x2, covparms = covparms));
      // Rprintf("det of kernel[%i %i]: %f \n", i, j, arma::det(kernel12));
      
      covmat(i,j) = sigma12 * spectralMatern(sqrt(q12), smooth12) / sqrt( arma::det(kernel12) );
    }
  }
  
  return( covmat );
}

//' @name cal_matern_ns
//'
//' @title For fast calculation of nonstationary Matern covariance (See cov_matern_ns)
//'
//' @param locs A numeric matrix of locations
//' @param sigvec A vector of variance parameters
//' @param nuvec A vector of smoothness parameters
//' @param kerList A list of kernel matrices
//'
//' @return The nonstationary Matern covariance matrix of the locations provided
// [[Rcpp::export]]
arma::mat cal_matern_ns(const arma::mat & locs, const arma::rowvec & sigvec, const arma::rowvec & nuvec, const Rcpp::List& kerList) {

  unsigned int n = locs.n_rows;
  unsigned int d = locs.n_cols;
  
  double sigma;
  double smooth;
  double q;
  
  arma::colvec tarvec( d );
  arma::mat kernel( d, d );

  arma::mat covmat( n, n ); covmat.fill(0.0);
  for(unsigned int i = 0; i < n; i++) {
    for(unsigned int j = i; j < n; j++) {
      
      tarvec = arma::conv_to< arma::colvec >::from(locs.row(i) - locs.row(j));
      
      sigma = sigvec(i) * sigvec(j);
      smooth = ( nuvec(i) + nuvec(j) ) / 2.0;
      
      kernel = ( Rcpp::as<arma::mat>(kerList[i]) + Rcpp::as<arma::mat>(kerList[j]) ) / 2.0;
      q = arma::dot( tarvec, arma::solve(kernel, tarvec) );
      
      covmat(i,j) = sigma * spectralMatern(sqrt(q), smooth) / sqrt( arma::det(kernel) );
    }
  }
  
  return( covmat );
}

//' @name cov_matern_ns_cpp
//'
//' @title Nonstationary Matern covariance function (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param sigma A numeric function for spatially varying standard deviation \eqn{\sigma( loc )}
//' @param smoothness A numeric function for spatially varying smoothness \eqn{\nu( loc )}
//' @param kernel A matrix-valued function for spatially varying (local) geometric anisotropy \eqn{\Sigma( loc )}
//'
//' @return The nonstationary Matern covariance between the locations provided
// [[Rcpp::export]]
double cov_matern_ns_cpp(const arma::rowvec & x1, const arma::rowvec & x2, Rcpp::Function sigma, Rcpp::Function smoothness, Rcpp::Function kernel) {
  
  // https://stackoverflow.com/questions/17958168/rcpp-function-to-be-slower-than-same-r-function
  // https://stackoverflow.com/questions/27391472/passing-r-function-as-parameter-to-rcpp-function
  
  unsigned int d = x1.n_elem;
  
  double cov = -1.0;
  arma::colvec tarvec(d); tarvec = arma::conv_to< arma::colvec >::from(x1 - x2); 
  
  double sigma12 = Rcpp::as<double>(sigma(x1)) * Rcpp::as<double>(sigma(x2));
  double smooth12 = ( Rcpp::as<double>(smoothness(x1)) + Rcpp::as<double>(smoothness(x2)) ) / 2.0;
  arma::mat kernel12( d, d ); kernel12 = ( Rcpp::as<arma::mat>(kernel(x1)) + Rcpp::as<arma::mat>(kernel(x2)) ) / 2.0;
  
  double q12 = arma::dot( tarvec, arma::solve(kernel12, tarvec) );
  
  arma::rowvec y1(1); y1 = {0.0};
  arma::rowvec y2(1); y2 = {sqrt(q12)}; 
  arma::rowvec covparms(3); covparms(0) = 1.0; covparms(1) = 1.0; covparms(2) = smooth12;
  
  cov = sigma12 * cov_matern_iso_cpp(y1, y2, covparms) / sqrt( arma::det(kernel12) );
  
  return(cov);
}

//' @name cov_latentDim_biv_cpp
//'
//' @title Isotropic exponential covariance function for bivariate Gaussian Processes using latent dimensions (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters = (sigma2, range, latent dimension, x1's id, x2's id). ID can be either 1 or 2.
//'
//' @return The isotropic exponential covariance between the locations provided with respect to latent dimension and their ids
// [[Rcpp::export]]
double cov_latentDim_biv_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec covparms) {
  
  unsigned int d = x1.n_elem;
  
  arma::rowvec newx1 = x1; newx1.resize(d + 1); newx1(d) = 0.0;
  arma::rowvec newx2 = x2; newx2.resize(d + 1);

  if(covparms(3) == covparms(4)) {
    
    newx2(d) = 0.0;
    
  } else {
    
    newx2(d) = covparms(2);
  }
  
  return( covparms(0) * exp(- norm(newx1 - newx2) / covparms(1)) );
}

//' @name cov_latentDim_triv_cpp
//'
//' @title Isotropic exponential covariance function for trivariate Gaussian Processes using latent dimensions (C++ version)
//'
//' @param x1 First location (vector)
//' @param x2 Second location (vector)
//' @param covparms A numerical vector with covariance parameters = (sigma2, range, latent dimension, x1's id, x2's id). ID can be 1 or 2, or 3.
//'
//' @return The isotropic exponential covariance between the locations provided with respect to latent dimension and their ids
// [[Rcpp::export]]
double cov_latentDim_triv_cpp(const arma::rowvec & x1, const arma::rowvec & x2, const arma::rowvec covparms) {
  
  unsigned int d = x1.n_elem;
  
  arma::rowvec newx1 = x1; newx1.resize(d + 1); newx1(d) = 0.0;
  arma::rowvec newx2 = x2; newx2.resize(d + 1);
  
  if(covparms(3) == covparms(4)) {
    
    newx2(d) = 0.0;
    
  } else if(abs(covparms(3) - covparms(4)) <= 1.0) {
    
    newx2(d) = covparms(2);
    
  } else {
    
    newx2(d) = 2.0 * covparms(2);
  }
  
  return( covparms(0) * exp(- norm(newx1 - newx2) / covparms(1)) );
}

