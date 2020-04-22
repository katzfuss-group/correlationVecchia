####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to perform parallel computing for the simulations.
###
###   Contents:
###
####################################################################################



#' @title Simulation on anisotropic cases
#'
#' @param cand.m A numeric vector of candidates of the size of conditioning sets (m)
#' @param cand.a A numeric vector of candidates of the degree of anisotropy (a) 
#' @param nsim A number of repeated simualtions for each case
#' @param n A number of locations
#' @param d A dimension of domain
#' @param covmodel A covariance function
#' @param covparms A numeric vector of covariance parameters
#' @param method.locs random or grid 
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default 
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param ncores A number of cores for parallel computing 
#'
#' @return list
#' 
#' @import foreach
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' out <- parallel_simulate_anisotropic_knownCovparms(cand.m = c(10, 20), 
#'                                                    cand.a = c(10, 20), 
#'                                                    nsim = 2, n = 10^2, d = 2)
#' out$kldiv
#' }
parallel_simulate_anisotropic_knownCovparms <- function(cand.m, cand.a, nsim, n, d, covmodel = cov_expo_aniso, covparms = c(1, 0.1), method.locs = 'random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()
  
  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")
  
  ### cand.all
  cand.all              <- expand.grid(cand.m, cand.a)
  cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
  colnames(cand.all)    <- c("index", "m", "a")
  
  ### ncores
  if(is.null(ncores)) {
    no_cores  <- parallel::detectCores() - 2
  } else {
    no_cores  <- ncores
  }
  
  ### simulation
  sim                   <- list()
  cl                    <- parallel::makeCluster(no_cores)
  
  doParallel::registerDoParallel(cl)
  sim <- foreach::foreach(m = cand.all$m, a = cand.all$a, .packages = c("correlationVecchia", "GPvecchia")) %dopar% simulate_univariate_knownCovparms(nsim = nsim, n = n, d = d, m = m, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = covparms, a = c(a, rep(1, d - 1)) )
  parallel::stopCluster(cl)
  
  ### KL divergence
  n.approx        <- sim[[1]]$n.approx
  kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
  kls[, 1]        <- seq(nrow(cand.all))
  for(k in 1:length(sim)) {
    kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
  }
  kls             <- as.data.frame(kls)
  colnames(kls)   <- c("index", names(sim[[1]]$kls.average))
  
  time.tot <- proc.time() - time.tot
  
  ### return
  return(list(vars = cand.all, kldiv = kls, time.tot = time.tot, nsim = nsim, n = n, d = d))
}



#' @title Simulation on nonstationary cases
#'
#' @param cand.m A numeric vector of candidates of the size of conditioning sets (m)
#' @param nsim A number of repeated simualtions for each case
#' @param n A number of locations
#' @param d A dimension of domain 
#' @param covmodel A covariance function
#' @param sigma  numeric function for spatially varying standard deviation \eqn{\sigma( loc )}
#' @param smoothness A numeric function for spatially varying smoothness \eqn{\nu( loc )}
#' @param kernel A matrix-valued function for spatially varying (local) geometric anisotropy \eqn{\Sigma( loc )}  
#' @param method.locs random or grid 
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default 
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param ncores A number of cores for parallel computing 
#'
#' @return list
#' 
#' @import foreach
#' 
#' @export
#'
#' @examples
#' ### comparison ###
#' 
#' locs <- matrix(runif(100), 50, 2)
#' 
#' kernel <- function(loc) {
#'   
#'   d         <- length(loc)
#'   
#'   a         <- function(loc) 10 * 10 # please use your own a function for the first coordinate
#'   b         <- function(loc) 10 # please use your own a function for the second coordinate
#'   angle     <- function(loc) 0 # please use your own spatially varying rotation angle
#'   
#'   eta       <- angle(loc)
#'   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
#'   
#'   range     <- c(a(loc)^(-2), b(loc)^(-2))
#'   
#'   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#' }
#' 
#' sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25
#' 
#' smoothness <- function(loc) 0.5
#' 
#' mat1 <- cov_matern_ns_bruteforce(locs1 = locs, locs2 = NULL, 
#'                                  sigma = sigma, smoothness = smoothness, kernel = kernel)
#' mat2 <- cov_expo_aniso(locs = locs, covparms = c(1, 0.1), a = c(10, 1))
#' 
#' sqrt(sum((mat1 - mat2)^2))
#' 
#' ### sim 1: smoothness ###
#' 
#' kernel <- function(loc) {
#'   
#'   d         <- length(loc)
#'   
#'   a         <- function(loc) 10
#'   b         <- function(loc) 10
#'   angle     <- function(loc) 0
#'   
#'   eta       <- angle(loc)
#'   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
#'   
#'   range     <- c(a(loc)^(-2), b(loc)^(-2))
#'   
#'   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#' }
#' 
#' sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25
#' 
#' smoothness <- function(loc) 0.2 + 1.3 * loc[1]
#' 
#'\dontrun{
#'  out <- parallel_simulate_nonstationary_knownCovparms(cand.m = c(10, 20), 
#'                                                       nsim = 2, n = 10^2, d = 2, 
#'                                                       covmodel = cov_matern_ns_bruteforce, 
#'                                                       sigma = sigma, smoothness = smoothness, 
#'                                                       kernel = kernel)
#'  out$kldiv
#'}
#' 
#' 
#' ### sim 2: range ###
#' 
#' kernel <- function(loc) {
#'   
#'   d         <- length(loc)
#'   
#'   a         <- function(loc) 10 * (1 + 10 * loc[1])
#'   b         <- function(loc) 10 * (1 + 10 * loc[1]) 
#'   angle     <- function(loc) 0 
#'   
#'   eta       <- angle(loc)
#'   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
#'   
#'   range     <- c(a(loc)^(-2), b(loc)^(-2))
#'   
#'   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#' }
#' 
#' sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25
#' 
#' smoothness <- function(loc) 0.5
#' 
#'\dontrun{
#'  out <- parallel_simulate_nonstationary_knownCovparms(cand.m = c(10, 20), 
#'                                                       nsim = 2, n = 10^2, d = 2, 
#'                                                       covmodel = cov_matern_ns_bruteforce, 
#'                                                       sigma = sigma, smoothness = smoothness, 
#'                                                       kernel = kernel)
#'  out$kldiv
#'}
#' 
#' 
#' ### sim 2 - 1: range ###
#' 
#' kernel <- function(loc) {
#'   
#'   d         <- length(loc)
#'   
#'   a         <- function(loc) 10/(0.1 + 0.9 * loc[1]) 
#'   b         <- function(loc) 10/(0.1 + 0.9 * loc[1]) 
#'   angle     <- function(loc) 0 
#'   
#'   eta       <- angle(loc)
#'   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
#'   
#'   range     <- c(a(loc)^(-2), b(loc)^(-2))
#'   
#'   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#' }
#' 
#' sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25
#' 
#' smoothness <- function(loc) 0.5
#' 
#'\dontrun{
#'  out <- parallel_simulate_nonstationary_knownCovparms(cand.m = c(10, 20), 
#'                                                       nsim = 2, n = 10^2, d = 2, 
#'                                                       covmodel = cov_matern_ns_bruteforce, 
#'                                                       sigma = sigma, smoothness = smoothness, 
#'                                                       kernel = kernel)
#'  out$kldiv
#'}
#' 
#' 
#' ### sim 3: anisotropy ###
#' 
#' kernel <- function(loc) {
#'   
#'   d         <- length(loc)
#'   
#'   a         <- function(loc) 10 * (1 + 10 * loc[1]) 
#'   b         <- function(loc) 10 
#'   angle     <- function(loc) 0 
#'   
#'   eta       <- angle(loc)
#'   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
#'   
#'   range     <- c(a(loc)^(-2), b(loc)^(-2))
#'   
#'   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#' }
#' 
#' sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25
#' 
#' smoothness <- function(loc) 0.5
#' 
#'\dontrun{
#'  out <- parallel_simulate_nonstationary_knownCovparms(cand.m = c(10, 20), 
#'                                                       nsim = 2, n = 10^2, d = 2, 
#'                                                       covmodel = cov_matern_ns_bruteforce, 
#'                                                       sigma = sigma, smoothness = smoothness, 
#'                                                       kernel = kernel)
#'  out$kldiv
#'}
#' 
#' 
#' ### sim 4: rotation ###
#' 
#' kernel <- function(loc) {
#'   
#'   d         <- length(loc)
#'   
#'   a         <- function(loc) 10 * 10 
#'   b         <- function(loc) 10 
#'   angle     <- function(loc) pi * loc[1] / 2
#'   
#'   eta       <- angle(loc)
#'   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
#'   
#'   range     <- c(a(loc)^(-2), b(loc)^(-2))
#'   
#'   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#' }
#' 
#' sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25
#' 
#' smoothness <- function(loc) 0.5
#' 
#'\dontrun{
#'  out <- parallel_simulate_nonstationary_knownCovparms(cand.m = c(10, 20), 
#'                                                       nsim = 2, n = 10^2, d = 2, 
#'                                                       covmodel = cov_matern_ns_bruteforce, 
#'                                                       sigma = sigma, smoothness = smoothness, 
#'                                                       kernel = kernel)
#'  out$kldiv
#'}
parallel_simulate_nonstationary_knownCovparms <- function(cand.m, nsim, n, d, covmodel = cov_matern_ns_bruteforce, sigma, smoothness, kernel, method.locs = 'random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()
  
  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")
  
  ### cand.all
  cand.all              <- data.frame(cand.m)
  cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
  colnames(cand.all)    <- c("index", "m")
  
  ### ncores
  if(is.null(ncores)) {
    no_cores  <- parallel::detectCores() - 2
  } else {
    no_cores  <- ncores
  }
  
  ### simulation
  sim                   <- list()
  cl                    <- parallel::makeCluster(no_cores)
  
  doParallel::registerDoParallel(cl)
  sim <- foreach::foreach(m = cand.all$m, .packages = c("correlationVecchia", "GPvecchia")) %dopar% simulate_univariate_knownCovparms(nsim = nsim, n = n, d = d, m = m, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, sigma = sigma, smoothness = smoothness, kernel = kernel)
  parallel::stopCluster(cl)
  
  ### KL divergence
  n.approx        <- sim[[1]]$n.approx
  kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
  kls[, 1]        <- seq(nrow(cand.all))
  for(k in 1:length(sim)) {
    kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
  }
  kls             <- as.data.frame(kls)
  colnames(kls)   <- c("index", names(sim[[1]]$kls.average))
  
  time.tot <- proc.time() - time.tot
  
  ### return
  return(list(vars = cand.all, kldiv = kls, time.tot = time.tot, nsim = nsim, n = n, d = d))
}



#' @title Simulation on multivariate cases
#'
#' @param cand.m A numeric vector of candidates of the size of conditioning sets (m)
#' @param cand.d A numeric vector of candidates of the distance of latent coordiante (d)
#' @param nsim A number of repeated simualtions for each case
#' @param n A number of locations
#' @param d A dimension of domain
#' @param p A number of processes 
#' @param covmodel A covariance function
#' @param covparms A numeric vector of covariance parameters 
#' @param method.locs "random", "overlap", or "grid" 
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default 
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param ncores A number of cores for parallel computing 
#'
#' @return list
#' 
#' @import foreach
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' out <- parallel_simulate_multivariate_knownCovparms(cand.m = c(10, 20), 
#'                                                     cand.d = c(1, 2), 
#'                                                     nsim = 2, n = 10^2, 
#'                                                     d = 2, p = 2, 
#'                                                     covmodel = cov_bivariate_expo_latDim, 
#'                                                     covparms = c(1, 0.1))
#' out$kldiv
#' }
parallel_simulate_multivariate_knownCovparms <- function(cand.m, cand.d, nsim, n, d, p = 2, covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = 'random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()
  
  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")
  
  ### cand.all
  cand.all              <- expand.grid(cand.m, cand.d)
  cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
  colnames(cand.all)    <- c("index", "m", "d")
  
  ### ncores
  if(is.null(ncores)) {
    no_cores  <- parallel::detectCores() - 2
  } else {
    no_cores  <- ncores
  }
  
  ### simulation
  sim                   <- list()
  cl                    <- parallel::makeCluster(no_cores)
  
  doParallel::registerDoParallel(cl)
  sim <- foreach::foreach(m = cand.all$m, d.latent = cand.all$d, .packages = c("correlationVecchia", "GPvecchia")) %dopar% simulate_multivariate_knownCovparms(nsim = nsim, n = n, d = d, p = p, m = m, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = c(covparms, d.latent))
  parallel::stopCluster(cl)
  
  ### KL divergence
  n.approx        <- sim[[1]]$n.approx
  kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
  kls[, 1]        <- seq(nrow(cand.all))
  for(k in 1:length(sim)) {
    kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
  }
  kls             <- as.data.frame(kls)
  colnames(kls)   <- c("index", names(sim[[1]]$kls.average))
  
  time.tot <- proc.time() - time.tot
  
  ### return
  return(list(vars = cand.all, kldiv = kls, time.tot = time.tot, nsim = nsim, n = n, d = d))
}



#' @title Simulation on spacetime cases
#'
#' @param cand.m A numeric vector of candidates of the size of conditioning sets (m)
#' @param nsim A number of repeated simualtions for each case
#' @param n A number of locations
#' @param d A dimension of domain 
#' @param t.len A number of repeated measurement (= a number of different temporal locations)
#' @param covmodel A covariance function
#' @param covparms A numeric vector of covariance parameters
#' @param method.locs "all.random", "space.random.time.grid", "all.grid", or "satellite"
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param ncores A number of cores for parallel computing 
#'
#' @return list
#' 
#' @import foreach
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' out <- parallel_simulate_spacetime_knownCovparms(cand.m = c(10, 20), 
#'                                                  nsim = 2, n = 10^2, 
#'                                                  d = 2, t.len = 3, 
#'                                                  covmodel = cov_spacetime_expo, 
#'                                                  covparms = c(1, 0.75, 50, 25), 
#'                                                  method.locs = "all.random")
#' out$kldiv
#' }
parallel_simulate_spacetime_knownCovparms <- function(cand.m, nsim, n, d, t.len, covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25), method.locs = 'all.random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()
  
  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")
  
  ### cand.all
  cand.all              <- data.frame(cand.m)
  cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
  colnames(cand.all)    <- c("index", "m")
  
  ### ncores
  if(is.null(ncores)) {
    no_cores  <- parallel::detectCores() - 2
  } else {
    no_cores  <- ncores
  }
  
  ### simulation
  sim                   <- list()
  cl                    <- parallel::makeCluster(no_cores)
  
  doParallel::registerDoParallel(cl)
  sim <- foreach::foreach(m = cand.all$m, .packages = c("correlationVecchia", "GPvecchia")) %dopar% simulate_spacetime_knownCovparms(nsim = nsim, n = n, d = d, t.len = t.len, m = m, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = covparms)
  parallel::stopCluster(cl)
  
  ### KL divergence
  n.approx        <- sim[[1]]$n.approx
  kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
  kls[, 1]        <- seq(nrow(cand.all))
  for(k in 1:length(sim)) {
    kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
  }
  kls             <- as.data.frame(kls)
  colnames(kls)   <- c("index", names(sim[[1]]$kls.average))
  
  time.tot <- proc.time() - time.tot
  
  ### return
  return(list(vars = cand.all, kldiv = kls, time.tot = time.tot, nsim = nsim, n = n, d = d))
}



#' @title Simulation on derivative cases
#'
#' @param cand.m A numeric vector of candidates of the size of conditioning sets (m)
#' @param cand.r A numeric vector of candidates of range parameter (r)
#' @param nsim A number of repeated simualtions for each case
#' @param n A number of locations
#' @param d A dimension of domain. It must be 2 for now
#' @param covmodel A covariance function
#' @param covparms A numeric vector of covariance parameters = c(sigma2, range = NA)
#' @param method.locs random or grid 
#' @param abs.corr Logical at \code{FALSE} be default. If \code{TRUE} then distance = 1-|rho|. If \code{FALSE} then distane = 1-rho
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default 
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param ncores A number of cores for parallel computing 
#'
#' @return list
#' 
#' @import foreach
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' out <- parallel_simulate_derivative_knownCovparms(cand.m = c(10, 20), 
#'                                                   cand.r = c(0.1, 10), 
#'                                                   nsim = 2, n = 10^2, 
#'                                                   d = 2, pivot = FALSE, 
#'                                                   method.modify = "eigen-I", 
#'                                                   tol = 1e-6)
#' out$kldiv
#' }
parallel_simulate_derivative_knownCovparms <- function(cand.m, cand.r, nsim, n, d, abs.corr = FALSE, covmodel = cov_derivative_matern_2.5_2d, covparms = c(1, NA), method.locs = 'random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()
  
  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")
  
  ### cand.all
  cand.all              <- expand.grid(cand.m, cand.r)
  cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
  colnames(cand.all)    <- c("index", "m", "r")
  
  ### ncores
  if(is.null(ncores)) {
    no_cores  <- parallel::detectCores() - 2
  } else {
    no_cores  <- ncores
  }
  
  ### simulation
  sim                   <- list()
  cl                    <- parallel::makeCluster(no_cores)
  
  doParallel::registerDoParallel(cl)
  sim <- foreach::foreach(m = cand.all$m, r = cand.all$r, .packages = c("correlationVecchia", "GPvecchia")) %dopar% simulate_derivative_knownCovparms(nsim = nsim, n = n, d = d, m = m, method.locs = method.locs, abs.corr = abs.corr, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = c(covparms[1], r))
  parallel::stopCluster(cl)
  
  ### KL divergence
  n.approx        <- sim[[1]]$n.approx
  kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
  kls[, 1]        <- seq(nrow(cand.all))
  for(k in 1:length(sim)) {
    kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
  }
  kls             <- as.data.frame(kls)
  colnames(kls)   <- c("index", names(sim[[1]]$kls.average))
  
  time.tot <- proc.time() - time.tot
  
  ### return
  return(list(vars = cand.all, kldiv = kls, time.tot = time.tot, nsim = nsim, n = n, d = d))
}
