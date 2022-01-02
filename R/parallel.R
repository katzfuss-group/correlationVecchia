####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to use parallel computing for simulation.
###
###   Contents:
###     parSim_aniso_knownCovparms / parSim_nonst_knownCovparms / parSim_wave_knownCovparms
###     parSim_mulv_knownCovparms
###     parSim_sptm_knownCovparms
###     parSim_deriv_knownCovparms
###     parSim_sptm_Fisher
###     parSim_sptm_Bayes
###     parSim_sptm_predict
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
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|cor|. If \code{FALSE} then distane = 1-cor. At \code{FALSE} by default
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
#' @export
#'
#' @examples
#' \dontrun{
#'   out <- parSim_aniso_knownCovparms(cand.m = c(10, 20), cand.a = c(10, 20),
#'                                     nsim = 2, n = 10^2, d = 2)
#'   out$kldiv
#' }
parSim_aniso_knownCovparms <- function(cand.m, cand.a, nsim, n, d, covmodel = cov_expo_aniso, covparms = c(1, 0.1), abs.corr = FALSE, method.locs = 'random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()

  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")

  `%do%` <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`

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
  sim <- foreach::foreach(m = cand.all$m, a = cand.all$a, .packages = c("correlationVecchia")) %dopar% simulate_univ( nsim = nsim, n = n, d = d, m = m, abs.corr = abs.corr, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = covparms, a = c(a, rep(1, d-1)) )
  parallel::stopCluster(cl)

  ### KL divergence
  n.approx              <- sim[[1]]$n.approx
  kls                   <- matrix(NA, nrow(cand.all), 1 + n.approx)
  kls[, 1]              <- seq(nrow(cand.all))
  for(k in 1:length(sim)) {
    kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
  }
  kls                   <- as.data.frame(kls)
  colnames(kls)         <- c("index", names(sim[[1]]$kls.average))

  time.tot              <- proc.time() - time.tot

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
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|cor|. If \code{FALSE} then distane = 1-cor. At \code{FALSE} by default
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
#' @export
#'
#' @examples
#' \dontrun{
#'   kernel <- function(loc) {
#'
#'     d         <- length(loc)
#'
#'     a         <- function(loc) 10
#'     b         <- function(loc) 10
#'     angle     <- function(loc) 0
#'
#'     eta       <- angle(loc)
#'     rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)),
#'                         nrow = d, ncol = d, byrow = TRUE)
#'
#'     range     <- c(a(loc)^(-2), b(loc)^(-2))
#'
#'     return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#'   }
#'
#'   sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25
#'
#'   smoothness <- function(loc) 0.2 + 1.3 * loc[1]
#'
#'   out <- parSim_nonst_knownCovparms(cand.m = c(10, 20), nsim = 2, n = 10^2, d = 2,
#'                                     covmodel = cov_matern_ns,
#'                                     sigma = sigma, smoothness = smoothness, kernel = kernel)
#'   out$kldiv
#' }
parSim_nonst_knownCovparms <- function(cand.m, nsim, n, d, covmodel = cov_matern_ns, sigma, smoothness, kernel, abs.corr = FALSE, method.locs = 'random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()

  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")

  `%do%` <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`

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
  sim <- foreach::foreach(m = cand.all$m, .packages = c("correlationVecchia")) %dopar% simulate_univ(nsim = nsim, n = n, d = d, m = m, abs.corr = abs.corr, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, sigma = sigma, smoothness = smoothness, kernel = kernel)
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
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|cor|. If \code{FALSE} then distane = 1-cor. At \code{FALSE} by default
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
#' @export
#'
#' @examples
#' \dontrun{
#'   out <- parSim_mulv_knownCovparms(cand.m = c(10, 20),
#'                                    cand.d = c(1, 2),
#'                                    nsim = 2, n = 10^2,
#'                                    d = 2, p = 2,
#'                                    covmodel = cov_latentDim_biv,
#'                                    covparms = c(1, 0.1))
#'   out$kldiv
#' }
parSim_mulv_knownCovparms <- function(cand.m, cand.d, nsim, n, d, p = 2, covmodel = cov_latentDim_biv, covparms = c(1, 0.1), abs.corr = FALSE, method.locs = 'random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()

  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")

  `%do%` <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`

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
  sim <- foreach::foreach(m = cand.all$m, d.latent = cand.all$d, .packages = c("correlationVecchia")) %dopar% simulate_mulv( nsim = nsim, n = n, d = d, p = p, m = m, abs.corr = abs.corr, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = c(covparms, d.latent) )
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
#' @param t A number of repeated measurement (= a number of different temporal locations)
#' @param covmodel A covariance function
#' @param covparms A numeric vector of covariance parameters
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|cor|. If \code{FALSE} then distane = 1-cor. At \code{FALSE} by default
#' @param method.locs "random", "monitoring", "grid", or "satellite"
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
#' @export
#'
#' @examples
#' \dontrun{
#'   out <- parSim_sptm_knownCovparms(cand.m = c(10, 20),
#'                                    nsim = 2, n = 10^2,
#'                                    d = 2, t = 2,
#'                                    covmodel = cov_expo_spacetime,
#'                                    covparms = c(1, 0.1, 1.0))
#'   out$kldiv
#' }
parSim_sptm_knownCovparms <- function(cand.m, nsim, n, d, t, covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0), abs.corr = FALSE, method.locs = 'random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()

  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")

  `%do%` <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`

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
  sim <- foreach::foreach(m = cand.all$m, .packages = c("correlationVecchia")) %dopar% simulate_spacetime(nsim = nsim, n = n, d = d, t = t, m = m, abs.corr = abs.corr, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = covparms)
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
#' @param covtype sqexpo or matern
#' @param covparms A numeric vector of covariance parameters
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|cor|. If \code{FALSE} then distane = 1-cor. At \code{FALSE} by default
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
#' @export
#'
#' @examples
#' \dontrun{
#'   out <- parSim_deriv_knownCovparms(cand.m = c(10, 20), cand.r = c(0.1, 10),
#'                                     nsim = 2, n = 10^2, d = 2,
#'                                     covmodel = cov_sqexpo_deriv, covparms = c(1),
#'                                     abs.corr = FALSE, method.locs = 'random',
#'                                     method.modify = "eigen-I", pivot = FALSE,
#'                                     tol = 1e-6, ncores = NULL)
#'   out$kldiv
#' }
parSim_deriv_knownCovparms <- function(cand.m, cand.r, nsim, n, d, covmodel = cov_sqexpo_deriv, covtype = "sqexpo", covparms = c(1), abs.corr = FALSE, method.locs = 'random', method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()

  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")

  `%do%` <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`

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
  if(covtype == "sqexpo") {

    sim                   <- list()
    cl                    <- parallel::makeCluster(no_cores)

    doParallel::registerDoParallel(cl)
    sim <- foreach::foreach(m = cand.all$m, r = cand.all$r, .packages = c("correlationVecchia")) %dopar% simulate_deriv(nsim = nsim, n = n, d = d, m = m, abs.corr = abs.corr, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = c(covparms, r))
    parallel::stopCluster(cl)

  } else if(covtype == "matern") {

    sim                   <- list()
    cl                    <- parallel::makeCluster(no_cores)

    doParallel::registerDoParallel(cl)
    sim <- foreach::foreach(m = cand.all$m, r = cand.all$r, .packages = c("correlationVecchia")) %dopar% simulate_deriv(nsim = nsim, n = n, d = d, m = m, abs.corr = abs.corr, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = c(covparms[1], r, covparms[length(covparms)]))
    parallel::stopCluster(cl)

  } else {

    stop("Check the argument covtype!")
  }

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

#' @title Simulation on periodic cases
#'
#' @param cand.m A numeric vector of candidates of the size of conditioning sets (m)
#' @param cand.l A numeric vector of candidates of period or nu
#' @param nsim A number of repeated simualtions for each case
#' @param n A number of locations
#' @param d A dimension of domain
#' @param covmodel A covariance function
#' @param covtype A method (function) to construct a wave model. It can be "Dampedsine," "Dampedcosine," or "BesselJ." At "Dampedsine" by default
#' @param covparms A numeric vector of covariance parameters
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|cor|. If \code{FALSE} then distane = 1-cor. At \code{FALSE} by default
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
#' @export
#'
#' @examples
#' \dontrun{
#'   out <- parSim_wave_knownCovparms(cand.m = c(10, 20),
#'                                    cand.l = c(0.1, 0.2),
#'                                    nsim = 2, n = 10^2, d = 2,
#'                                    covmodel = cov_wave, covtype = "Dampedsine",
#'                                    covparms = c(1), abs.corr = TRUE,
#'                                    method.modify = "eigen-I", tol = 1e-5)
#'   out$kldiv
#' }
parSim_wave_knownCovparms <- function(cand.m, cand.l, nsim, n, d, covmodel = cov_wave, covtype = "Dampedsine", covparms = c(1), abs.corr = TRUE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  time.tot <- proc.time()

  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")

  `%do%` <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`

  ### cand.all
  cand.all              <- expand.grid(cand.m, cand.l)
  cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
  colnames(cand.all)    <- c("index", "m", "period")

  ### ncores
  if(is.null(ncores)) {
    no_cores  <- parallel::detectCores() - 2
  } else {
    no_cores  <- ncores
  }

  covftn <- function(locs, covparms) covmodel(locs = locs, covparms = covparms, method = covtype)

  ### simulation
  sim                   <- list()
  cl                    <- parallel::makeCluster(no_cores)

  doParallel::registerDoParallel(cl)
  sim <- foreach::foreach(m = cand.all$m, period = cand.all$period, .packages = c("correlationVecchia")) %dopar% simulate_univ_cpt(nsim = nsim, n = n, d = d, m = m, abs.corr = abs.corr, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covftn, covparms = c(covparms, period))
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

#' @title Simulation on Fisher-scoring approach
#'
#' @param cand.m A numeric vector of candidates of the size of conditioning sets (m)
#' @param nsim A number of repeated simualtions for each case
#' @param n A number of locations
#' @param d A dimension of domain
#' @param t A number of repeated measurement (= a number of different temporal locations)
#' @param covmodel A covariance function
#' @param srange.ini Initial value of spatial range parameter
#' @param trange.ini Initial value of temporal range parameter
#' @param covparms A numeric vector of covariance parameters
#' @param method.locs "random", "monitoring", "grid", or "satellite"
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param ncores A number of cores for parallel computing
#' @param tol.dec tol.dec
#'
#' @return List
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   out <- parSim_sptm_Fisher(cand.m = c(10, 30), nsim = 2, n = 20^2, d = 2, t = 2,
#'                             covmodel = GpGp::matern_spacetime,
#'                             srange.ini = 0.5, trange.ini = 0.5,
#'                             covparms = c(1, 0.1, 1.0, 0.5, 0),
#'                             method.locs = "random",
#'                             method.modify = NULL, pivot = FALSE,
#'                             tol = .Machine$double.eps, ncores = NULL)
#'   out$kls.average
#' }
parSim_sptm_Fisher <- function(cand.m, nsim, n, d, t, covmodel = GpGp::matern_spacetime, srange.ini = 0.5, trange.ini = 0.5, covparms = c(1, 0.1, 1.0, 0.5, 0), method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL, tol.dec = 4)
{
  verbose     <- 0
  time.tot    <- proc.time()

  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")

  `%do%`      <- foreach::`%do%`
  `%dopar%`   <- foreach::`%dopar%`

  ### ncores
  if(is.null(ncores)) {
    no_cores    <- parallel::detectCores() - 2
  } else {
    no_cores    <- ncores
  }

  ### generate
  covftn      <- function(locs, covparms) covmodel(covparms = covparms, locs = locs)
  gen         <- generate_gp(nsim = nsim, n = n, d = d, domain = "spacetime", p = NULL, t = t, method.locs = method.locs, covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0))

  locs        <- list()
  for(i in 1:nsim) locs[[i]] <- as.matrix(gen$sim[[i]]$locs)

  rm(gen)

  covmat      <- list()
  for(i in 1:nsim) covmat[[i]] <- covftn(locs = locs[[i]], covparms = covparms)

  y           <- list()
  for(i in 1:nsim) y[[i]] <- t( factorize(covmat = covmat[[i]], pivot = pivot, method = method.modify, tol = tol, return.err = FALSE, verbose = FALSE)$covfactor ) %*% stats::rnorm(n = nrow(covmat[[i]]))

  ### fitting
  cl        <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)

  fit_b1 <- list()
  for(i in 1:length(cand.m)) {
    fit_b1[[i]] <- list()
    fit_b1[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia")) %dopar% fit_scaled_spacetime_range(method = "b1", y = y[[j]], inputs = locs[[j]], ms = cand.m[i], sig2 = covparms[1], srange.ini = srange.ini, trange.ini = trange.ini, nu = covparms[4], nug = 0, print.level = verbose, max.it = 100, tol.dec = tol.dec)

    if(verbose != 0) print(i)
  }

  fit_b2 <- list()
  for(i in 1:length(cand.m)) {
    fit_b2[[i]] <- list()
    fit_b2[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia")) %dopar% fit_scaled_spacetime_range(method = "b2", y = y[[j]], inputs = locs[[j]], ms = cand.m[i], sig2 = covparms[1], srange.ini = srange.ini, trange.ini = trange.ini, nu = covparms[4], nug = 0, print.level = verbose, max.it = 100, tol.dec = tol.dec)

    if(verbose != 0) print(i)
  }

  fit_b3 <- list()
  for(i in 1:length(cand.m)) {
    fit_b3[[i]] <- list()
    fit_b3[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia")) %dopar% fit_scaled_spacetime_range(method = "b3", y = y[[j]], inputs = locs[[j]], ms = cand.m[i], sig2 = covparms[1], srange.ini = srange.ini, trange.ini = trange.ini, nu = covparms[4], nug = 0, print.level = verbose, max.it = 100, tol.dec = tol.dec)

    if(verbose != 0) print(i)
  }

  fit_cc <- list()
  for(i in 1:length(cand.m)) {
    fit_cc[[i]] <- list()
    fit_cc[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia")) %dopar% fit_scaled_spacetime_range(method = "cc", y = y[[j]], inputs = locs[[j]], ms = cand.m[i], sig2 = covparms[1], srange.ini = srange.ini, trange.ini = trange.ini, nu = covparms[4], nug = 0, print.level = verbose, max.it = 100, tol.dec = tol.dec)

    if(verbose != 0) print(i)
  }

  fit_gp <- list()
  fit_gp <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia")) %dopar% fit_scaled_spacetime_range(method = "gp", y = y[[j]], inputs = locs[[j]], ms = cand.m[i], sig2 = covparms[1], srange.ini = srange.ini, trange.ini = trange.ini, nu = covparms[4], nug = 0, print.level = verbose, max.it = 100, tol.dec = tol.dec)

  fit_temp <- list()
  for(i in 1:length(cand.m)) {
    fit_temp[[i]] <- list()
    for(j in 1:nsim) fit_temp[[i]][[j]] <- fit_gp[[j]]

    if(verbose != 0) print(i)
  }

  parallel::stopCluster(cl)
  fit_gp <- fit_temp ; rm(fit_temp)

  ### KL divergence
  approxs <- c("b1", "b2", "b3", "cc")
  fit_list <- list()
  fit_list[[1]] <- fit_b1 ; rm(fit_b1)
  fit_list[[2]] <- fit_b2 ; rm(fit_b2)
  fit_list[[3]] <- fit_b3 ; rm(fit_b3)
  fit_list[[4]] <- fit_cc ; rm(fit_cc)

  cl        <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)

  kls.full <- foreach::foreach(j = seq(length(approxs)), .packages = c("correlationVecchia")) %dopar% .loop_kldiv(nsim = nsim, ms = cand.m, approx = approxs[j], locs.ls = locs, covmat.ls = covmat, fit.ls = fit_list[[j]], verbose = FALSE)

  kls.gp <- list()
  for(i in 1:length(cand.m)) {

    kls.gp[[i]] <- list()
    kls.gp[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia")) %dopar% kldiv(covmat0 = covmat[[j]], covmat1 = covftn(locs = locs[[j]], covparms = fit_gp[[i]][[j]]$covparms), mu0 = rep(0, length(y[[j]])), mu1 = rep(fit_gp[[i]][[j]]$betahat, length(y[[j]])))
  }

  parallel::stopCluster(cl)

  fit_list[[5]] <- fit_gp ; rm(fit_gp)
  kls.full[[5]] <- kls.gp ; rm(kls.gp)

  meankld_b1 <- .loop_meankldiv(nsim = nsim, m = cand.m, kldiv.ls = kls.full[[1]])
  meankld_b2 <- .loop_meankldiv(nsim = nsim, m = cand.m, kldiv.ls = kls.full[[2]])
  meankld_b3 <- .loop_meankldiv(nsim = nsim, m = cand.m, kldiv.ls = kls.full[[3]])
  meankld_cc <- .loop_meankldiv(nsim = nsim, m = cand.m, kldiv.ls = kls.full[[4]])
  meankld_gp <- .loop_meankldiv(nsim = nsim, m = cand.m, kldiv.ls = kls.full[[5]])

  kls.average <- data.frame(index = seq(length(meankld_b1)), m = cand.m, approx_1 = meankld_b1, approx_2 = meankld_b2, approx_3 = meankld_b3, approx_4 = meankld_cc, approx_5 = meankld_gp)
  rm(meankld_b1) ; rm(meankld_b2) ; rm(meankld_b3) ; rm(meankld_cc) ; rm(meankld_gp)

  ### Likelihood
  lik_b1 <- list()
  for(i in 1:length(cand.m)) {

    lik_b1[[i]] <- list()
    for(j in 1:nsim) {
      Sigma <- covftn(locs = locs[[j]], covparms = fit_list[[1]][[i]][[j]]$covparms)
      vecchia.approx <- baseline_1_spacetime_specify(locs = locs[[j]], m = cand.m[i])

      lik_b1[[i]][[j]] <- .loop_loglik(y = y[[j]], locs = locs[[j]], fit = fit_list[[1]][[i]][[j]], approx = vecchia.approx, covmodel = covftn)
    }

    if(verbose != 0) print(i)
  }

  lik_b2 <- list()
  for(i in 1:length(cand.m)) {

    lik_b2[[i]] <- list()
    for(j in 1:nsim) {
      Sigma <- covftn(locs = locs[[j]], covparms = fit_list[[2]][[i]][[j]]$covparms)
      vecchia.approx <- baseline_2_spacetime_specify(locs = locs[[j]], m = cand.m[i])

      lik_b2[[i]][[j]] <- .loop_loglik(y = y[[j]], locs = locs[[j]], fit = fit_list[[2]][[i]][[j]], approx = vecchia.approx, covmodel = covftn)
    }

    if(verbose != 0) print(i)
  }

  lik_b3 <- list()
  for(i in 1:length(cand.m)) {

    lik_b3[[i]] <- list()
    for(j in 1:nsim) {
      Sigma <- covftn(locs = locs[[j]], covparms = fit_list[[3]][[i]][[j]]$covparms)
      vecchia.approx <- baseline_3_spacetime_specify(locs = locs[[j]], m = cand.m[i], covmodel = covftn, covparms = fit_list[[3]][[i]][[j]]$covparms)

      lik_b3[[i]][[j]] <- .loop_loglik(y = y[[j]], locs = locs[[j]], fit = fit_list[[3]][[i]][[j]], approx = vecchia.approx, covmodel = covftn)
    }

    if(verbose != 0) print(i)
  }

  lik_cc <- list()
  for(i in 1:length(cand.m)) {

    lik_cc[[i]] <- list()
    for(j in 1:nsim) {
      Sigma <- covftn(locs = locs[[j]], covparms = fit_list[[4]][[i]][[j]]$covparms)
      vecchia.approx <- cvecchia_m_specify(locs = locs[[j]], cand.m[i], covmodel = Sigma, covparms = fit_list[[4]][[i]][[j]]$covparms)

      lik_cc[[i]][[j]] <- .loop_loglik(y = y[[j]], locs = locs[[j]], fit = fit_list[[4]][[i]][[j]], approx = vecchia.approx, covmodel = covftn)
    }

    if(verbose != 0) print(i)
  }

  lik_gp <- list()
  for(i in 1:length(cand.m)) {

    lik_gp[[i]] <- list()
    for(j in 1:nsim) {
      Sigma <- covftn(locs = locs[[j]], covparms = fit_list[[5]][[i]][[j]]$covparms)
      vecchia.approx <- baseline_1_spacetime_specify(locs = locs[[j]], m = length(y[[j]])-1)

      lik_gp[[i]][[j]] <- .loop_loglik(y = y[[j]], locs = locs[[j]], fit = fit_list[[5]][[i]][[j]], approx = vecchia.approx, covmodel = covftn)
    }

    if(verbose != 0) print(i)
  }

  meanlik_b1 <- rep(0, length(cand.m))
  for(i in 1:length(cand.m)) {
    for(j in 1:nsim) {
      meanlik_b1[i] <- meanlik_b1[i] + lik_b1[[i]][[j]]
    }

    meanlik_b1[i] <- meanlik_b1[i] / nsim
  }

  meanlik_b2 <- rep(0, length(cand.m))
  for(i in 1:length(cand.m)) {
    for(j in 1:nsim) {
      meanlik_b2[i] <- meanlik_b2[i] + lik_b2[[i]][[j]]
    }

    meanlik_b2[i] <- meanlik_b2[i] / nsim
  }

  meanlik_b3 <- rep(0, length(cand.m))
  for(i in 1:length(cand.m)) {
    for(j in 1:nsim) {
      meanlik_b3[i] <- meanlik_b3[i] + lik_b3[[i]][[j]]
    }

    meanlik_b3[i] <- meanlik_b3[i] / nsim
  }

  meanlik_cc <- rep(0, length(cand.m))
  for(i in 1:length(cand.m)) {
    for(j in 1:nsim) {
      meanlik_cc[i] <- meanlik_cc[i] + lik_cc[[i]][[j]]
    }

    meanlik_cc[i] <- meanlik_cc[i] / nsim
  }

  meanlik_gp <- rep(0, length(cand.m))
  for(i in 1:length(cand.m)) {
    for(j in 1:nsim) {
      meanlik_gp[i] <- meanlik_gp[i] + lik_gp[[i]][[j]]
    }

    meanlik_gp[i] <- meanlik_gp[i] / nsim
  }

  lik.full <- list()
  lik.full[[1]] <- lik_b1 ; rm(lik_b1)
  lik.full[[2]] <- lik_b2 ; rm(lik_b2)
  lik.full[[3]] <- lik_b3 ; rm(lik_b3)
  lik.full[[4]] <- lik_cc ; rm(lik_cc)
  lik.full[[5]] <- lik_gp ; rm(lik_gp)

  lik.average <- data.frame(index = seq(length(meanlik_b1)), m = cand.m, approx_1 = meanlik_b1, approx_2 = meanlik_b2, approx_3 = meanlik_b3, approx_4 = meanlik_cc, approx_5 = meanlik_gp)
  rm(meanlik_b1) ; rm(meanlik_b2) ; rm(meanlik_b3) ; rm(meanlik_cc) ; rm(meanlik_gp)

  ### MSE
  mse_srange_b1 <- c()
  mse_trange_b1 <- c()
  for(i in 1:length(cand.m)) {

    mse_srange_b1[i] <- 0
    mse_trange_b1[i] <- 0
    for(j in 1:nsim) {
      mse_srange_b1[i] <- mse_srange_b1[i] + (fit_list[[1]][[i]][[j]]$covparms[2] - covparms[2])^2
      mse_trange_b1[i] <- mse_trange_b1[i] + (fit_list[[1]][[i]][[j]]$covparms[3] - covparms[3])^2
    }

    mse_srange_b1[i] <- mse_srange_b1[i] / nsim
    mse_trange_b1[i] <- mse_trange_b1[i] / nsim
  }

  mse_srange_b2 <- c()
  mse_trange_b2 <- c()
  for(i in 1:length(cand.m)) {

    mse_srange_b2[i] <- 0
    mse_trange_b2[i] <- 0
    for(j in 1:nsim) {
      mse_srange_b2[i] <- mse_srange_b2[i] + (fit_list[[2]][[i]][[j]]$covparms[2] - covparms[2])^2
      mse_trange_b2[i] <- mse_trange_b2[i] + (fit_list[[2]][[i]][[j]]$covparms[3] - covparms[3])^2
    }

    mse_srange_b2[i] <- mse_srange_b2[i] / nsim
    mse_trange_b2[i] <- mse_trange_b2[i] / nsim
  }

  mse_srange_b3 <- c()
  mse_trange_b3 <- c()
  for(i in 1:length(cand.m)) {

    mse_srange_b3[i] <- 0
    mse_trange_b3[i] <- 0
    for(j in 1:nsim) {
      mse_srange_b3[i] <- mse_srange_b3[i] + (fit_list[[3]][[i]][[j]]$covparms[2] - covparms[2])^2
      mse_trange_b3[i] <- mse_trange_b3[i] + (fit_list[[3]][[i]][[j]]$covparms[3] - covparms[3])^2
    }

    mse_srange_b3[i] <- mse_srange_b3[i] / nsim
    mse_trange_b3[i] <- mse_trange_b3[i] / nsim
  }

  mse_srange_cc <- c()
  mse_trange_cc <- c()
  for(i in 1:length(cand.m)) {

    mse_srange_cc[i] <- 0
    mse_trange_cc[i] <- 0
    for(j in 1:nsim) {
      mse_srange_cc[i] <- mse_srange_cc[i] + (fit_list[[4]][[i]][[j]]$covparms[2] - covparms[2])^2
      mse_trange_cc[i] <- mse_trange_cc[i] + (fit_list[[4]][[i]][[j]]$covparms[3] - covparms[3])^2
    }

    mse_srange_cc[i] <- mse_srange_cc[i] / nsim
    mse_trange_cc[i] <- mse_trange_cc[i] / nsim
  }

  mse_srange_gp <- c()
  mse_trange_gp <- c()
  for(i in 1:length(cand.m)) {

    mse_srange_gp[i] <- 0
    mse_trange_gp[i] <- 0
    for(j in 1:nsim) {
      mse_srange_gp[i] <- mse_srange_gp[i] + (fit_list[[5]][[i]][[j]]$covparms[2] - covparms[2])^2
      mse_trange_gp[i] <- mse_trange_gp[i] + (fit_list[[5]][[i]][[j]]$covparms[3] - covparms[3])^2
    }

    mse_srange_gp[i] <- mse_srange_gp[i] / nsim
    mse_trange_gp[i] <- mse_trange_gp[i] / nsim
  }

  mse_srange <- list()
  mse_srange[[1]] <- mse_srange_b1 ; rm(mse_srange_b1)
  mse_srange[[2]] <- mse_srange_b2 ; rm(mse_srange_b2)
  mse_srange[[3]] <- mse_srange_b3 ; rm(mse_srange_b3)
  mse_srange[[4]] <- mse_srange_cc ; rm(mse_srange_cc)
  mse_srange[[5]] <- mse_srange_gp ; rm(mse_srange_gp)

  mse_trange <- list()
  mse_trange[[1]] <- mse_trange_b1 ; rm(mse_trange_b1)
  mse_trange[[2]] <- mse_trange_b2 ; rm(mse_trange_b2)
  mse_trange[[3]] <- mse_trange_b3 ; rm(mse_trange_b3)
  mse_trange[[4]] <- mse_trange_cc ; rm(mse_trange_cc)
  mse_trange[[5]] <- mse_trange_gp ; rm(mse_trange_gp)

  ### MAE
  mae_srange_b1 <- c()
  mae_trange_b1 <- c()
  for(i in 1:length(cand.m)) {

    mae_srange_b1[i] <- 0
    mae_trange_b1[i] <- 0
    for(j in 1:nsim) {
      mae_srange_b1[i] <- mae_srange_b1[i] + abs(fit_list[[1]][[i]][[j]]$covparms[2] - covparms[2])
      mae_trange_b1[i] <- mae_trange_b1[i] + abs(fit_list[[1]][[i]][[j]]$covparms[3] - covparms[3])
    }

    mae_srange_b1[i] <- mae_srange_b1[i] / nsim
    mae_trange_b1[i] <- mae_trange_b1[i] / nsim
  }

  mae_srange_b2 <- c()
  mae_trange_b2 <- c()
  for(i in 1:length(cand.m)) {

    mae_srange_b2[i] <- 0
    mae_trange_b2[i] <- 0
    for(j in 1:nsim) {
      mae_srange_b2[i] <- mae_srange_b2[i] + abs(fit_list[[2]][[i]][[j]]$covparms[2] - covparms[2])
      mae_trange_b2[i] <- mae_trange_b2[i] + abs(fit_list[[2]][[i]][[j]]$covparms[3] - covparms[3])
    }

    mae_srange_b2[i] <- mae_srange_b2[i] / nsim
    mae_trange_b2[i] <- mae_trange_b2[i] / nsim
  }

  mae_srange_b3 <- c()
  mae_trange_b3 <- c()
  for(i in 1:length(cand.m)) {

    mae_srange_b3[i] <- 0
    mae_trange_b3[i] <- 0
    for(j in 1:nsim) {
      mae_srange_b3[i] <- mae_srange_b3[i] + abs(fit_list[[3]][[i]][[j]]$covparms[2] - covparms[2])
      mae_trange_b3[i] <- mae_trange_b3[i] + abs(fit_list[[3]][[i]][[j]]$covparms[3] - covparms[3])
    }

    mae_srange_b3[i] <- mae_srange_b3[i] / nsim
    mae_trange_b3[i] <- mae_trange_b3[i] / nsim
  }

  mae_srange_cc <- c()
  mae_trange_cc <- c()
  for(i in 1:length(cand.m)) {

    mae_srange_cc[i] <- 0
    mae_trange_cc[i] <- 0
    for(j in 1:nsim) {
      mae_srange_cc[i] <- mae_srange_cc[i] + abs(fit_list[[4]][[i]][[j]]$covparms[2] - covparms[2])
      mae_trange_cc[i] <- mae_trange_cc[i] + abs(fit_list[[4]][[i]][[j]]$covparms[3] - covparms[3])
    }

    mae_srange_cc[i] <- mae_srange_cc[i] / nsim
    mae_trange_cc[i] <- mae_trange_cc[i] / nsim
  }

  mae_srange_gp <- c()
  mae_trange_gp <- c()
  for(i in 1:length(cand.m)) {

    mae_srange_gp[i] <- 0
    mae_trange_gp[i] <- 0
    for(j in 1:nsim) {
      mae_srange_gp[i] <- mae_srange_gp[i] + abs(fit_list[[5]][[i]][[j]]$covparms[2] - covparms[2])
      mae_trange_gp[i] <- mae_trange_gp[i] + abs(fit_list[[5]][[i]][[j]]$covparms[3] - covparms[3])
    }

    mae_srange_gp[i] <- mae_srange_gp[i] / nsim
    mae_trange_gp[i] <- mae_trange_gp[i] / nsim
  }

  mae_srange <- list()
  mae_srange[[1]] <- mae_srange_b1 ; rm(mae_srange_b1)
  mae_srange[[2]] <- mae_srange_b2 ; rm(mae_srange_b2)
  mae_srange[[3]] <- mae_srange_b3 ; rm(mae_srange_b3)
  mae_srange[[4]] <- mae_srange_cc ; rm(mae_srange_cc)
  mae_srange[[5]] <- mae_srange_gp ; rm(mae_srange_gp)

  mae_trange <- list()
  mae_trange[[1]] <- mae_trange_b1 ; rm(mae_trange_b1)
  mae_trange[[2]] <- mae_trange_b2 ; rm(mae_trange_b2)
  mae_trange[[3]] <- mae_trange_b3 ; rm(mae_trange_b3)
  mae_trange[[4]] <- mae_trange_cc ; rm(mae_trange_cc)
  mae_trange[[5]] <- mae_trange_gp ; rm(mae_trange_gp)

  ### MAD, MSD
  vislist_srange <- list()
  vislist_trange <- list()
  for(j in 1:nsim) {
    vislist_srange[[j]] <- data.frame(nsim = j, m = cand.m, est_b1 = NA, est_b2 = NA, est_b3 = NA, est_cc = NA, est_gp = NA)
    vislist_trange[[j]] <- data.frame(nsim = j, m = cand.m, est_b1 = NA, est_b2 = NA, est_b3 = NA, est_cc = NA, est_gp = NA)
    for(i in 1:length(cand.m)) {
      vislist_srange[[j]]$est_b1[i] <- fit_list[[1]][[i]][[j]]$covparms[2]
      vislist_srange[[j]]$est_b2[i] <- fit_list[[2]][[i]][[j]]$covparms[2]
      vislist_srange[[j]]$est_b3[i] <- fit_list[[3]][[i]][[j]]$covparms[2]
      vislist_srange[[j]]$est_cc[i] <- fit_list[[4]][[i]][[j]]$covparms[2]
      vislist_srange[[j]]$est_gp[i] <- fit_list[[5]][[i]][[j]]$covparms[2]

      vislist_trange[[j]]$est_b1[i] <- fit_list[[1]][[i]][[j]]$covparms[3]
      vislist_trange[[j]]$est_b2[i] <- fit_list[[2]][[i]][[j]]$covparms[3]
      vislist_trange[[j]]$est_b3[i] <- fit_list[[3]][[i]][[j]]$covparms[3]
      vislist_trange[[j]]$est_cc[i] <- fit_list[[4]][[i]][[j]]$covparms[3]
      vislist_trange[[j]]$est_gp[i] <- fit_list[[5]][[i]][[j]]$covparms[3]
    }
  }

  visdat_srange <- vislist_srange[[1]]
  visdat_trange <- vislist_trange[[1]]
  for(j in 2:nsim) {
    visdat_srange <- rbind(visdat_srange, vislist_srange[[j]])
    visdat_trange <- rbind(visdat_trange, vislist_trange[[j]])
  }

  visdat_srange$ad_b1 <- abs(visdat_srange$est_b1 - visdat_srange$est_gp)
  visdat_srange$ad_b2 <- abs(visdat_srange$est_b2 - visdat_srange$est_gp)
  visdat_srange$ad_b3 <- abs(visdat_srange$est_b3 - visdat_srange$est_gp)
  visdat_srange$ad_cc <- abs(visdat_srange$est_cc - visdat_srange$est_gp)

  visdat_trange$ad_b1 <- abs(visdat_trange$est_b1 - visdat_trange$est_gp)
  visdat_trange$ad_b2 <- abs(visdat_trange$est_b2 - visdat_trange$est_gp)
  visdat_trange$ad_b3 <- abs(visdat_trange$est_b3 - visdat_trange$est_gp)
  visdat_trange$ad_cc <- abs(visdat_trange$est_cc - visdat_trange$est_gp)

  visdat_srange$sd_b1 <- (visdat_srange$est_b1 - visdat_srange$est_gp)^2
  visdat_srange$sd_b2 <- (visdat_srange$est_b2 - visdat_srange$est_gp)^2
  visdat_srange$sd_b3 <- (visdat_srange$est_b3 - visdat_srange$est_gp)^2
  visdat_srange$sd_cc <- (visdat_srange$est_cc - visdat_srange$est_gp)^2

  visdat_trange$sd_b1 <- (visdat_trange$est_b1 - visdat_trange$est_gp)^2
  visdat_trange$sd_b2 <- (visdat_trange$est_b2 - visdat_trange$est_gp)^2
  visdat_trange$sd_b3 <- (visdat_trange$est_b3 - visdat_trange$est_gp)^2
  visdat_trange$sd_cc <- (visdat_trange$est_cc - visdat_trange$est_gp)^2

  mad_srange <- list()

  mad_srange[[1]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_srange$m == cand.m[i]
    mad_srange[[1]][i]  <- mean(visdat_srange$ad_b1[idx])
  }

  mad_srange[[2]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_srange$m == cand.m[i]
    mad_srange[[2]][i]  <- mean(visdat_srange$ad_b2[idx])
  }

  mad_srange[[3]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_srange$m == cand.m[i]
    mad_srange[[3]][i]  <- mean(visdat_srange$ad_b3[idx])
  }

  mad_srange[[4]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_srange$m == cand.m[i]
    mad_srange[[4]][i]  <- mean(visdat_srange$ad_cc[idx])
  }

  msd_srange <- list()

  msd_srange[[1]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_srange$m == cand.m[i]
    msd_srange[[1]][i]  <- mean(visdat_srange$sd_b1[idx])
  }

  msd_srange[[2]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_srange$m == cand.m[i]
    msd_srange[[2]][i]  <- mean(visdat_srange$sd_b2[idx])
  }

  msd_srange[[3]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_srange$m == cand.m[i]
    msd_srange[[3]][i]  <- mean(visdat_srange$sd_b3[idx])
  }

  msd_srange[[4]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_srange$m == cand.m[i]
    msd_srange[[4]][i]  <- mean(visdat_srange$sd_cc[idx])
  }

  mad_trange <- list()

  mad_trange[[1]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_trange$m == cand.m[i]
    mad_trange[[1]][i]  <- mean(visdat_trange$ad_b1[idx])
  }

  mad_trange[[2]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_trange$m == cand.m[i]
    mad_trange[[2]][i]  <- mean(visdat_trange$ad_b2[idx])
  }

  mad_trange[[3]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_trange$m == cand.m[i]
    mad_trange[[3]][i]  <- mean(visdat_trange$ad_b3[idx])
  }

  mad_trange[[4]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_trange$m == cand.m[i]
    mad_trange[[4]][i]  <- mean(visdat_trange$ad_cc[idx])
  }

  msd_trange <- list()

  msd_trange[[1]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_trange$m == cand.m[i]
    msd_trange[[1]][i]  <- mean(visdat_trange$sd_b1[idx])
  }

  msd_trange[[2]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_trange$m == cand.m[i]
    msd_trange[[2]][i]  <- mean(visdat_trange$sd_b2[idx])
  }

  msd_trange[[3]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_trange$m == cand.m[i]
    msd_trange[[3]][i]  <- mean(visdat_trange$sd_b3[idx])
  }

  msd_trange[[4]] <- rep(NA, length(cand.m))
  for(i in 1:length(cand.m)) {
    idx                 <- visdat_trange$m == cand.m[i]
    msd_trange[[4]][i]  <- mean(visdat_trange$sd_cc[idx])
  }

  ### return
  time.tot              <- proc.time() - time.tot

  result                <- list()

  result$names          <- c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP")
  result$fit.all        <- fit_list

  result$kls.full       <- kls.full
  result$kls.average    <- kls.average
  result$lik.full       <- lik.full
  result$lik.average    <- lik.average

  result$mse.srange     <- mse_srange
  result$mse.trange     <- mse_trange
  result$mae.srange     <- mae_srange
  result$mae.trange     <- mae_trange
  result$msd.srange     <- msd_srange
  result$msd.trange     <- msd_trange
  result$mad.srange     <- mad_srange
  result$mad.trange     <- mad_trange

  result$time.tot       <- time.tot

  return( result )
}

#' @title parSim_sptm_prediction
#'
#' @param cand.m 1
#' @param nsim 1
#' @param n 1
#' @param n.pred 1
#' @param d 1
#' @param t 1
#' @param nuggets 1
#' @param covmodel 1
#' @param covparms 1
#' @param abs.corr 1
#' @param method.locs 1
#' @param method.locs.pred 1
#' @param method.modify 1
#' @param pivot 1
#' @param tol 1
#' @param ncores 1
#'
#' @return 1
#'
#' @export
#'
#' @examples
#' 1 + 1
parSim_sptm_prediction <- function(cand.m, nsim, n, n.pred, d, t, nuggets, covmodel, covparms, abs.corr, method.locs, method.locs.pred, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
{
  ### initialization
  if(is.function(covmodel)) {

    covname <- deparse(substitute(covmodel))

  } else if(is.character(covmodel)) {

    covname <- covmodel

  } else {

    stop("Please use a function or string to specify the argument covmodel!")
  }

  if( covname %in% c("cov_expo_spacetime", "cov_expo_spacetime_cpp") ) {

    covmodel  <- cov_expo_spacetime
    covftn    <- "cov_expo_spacetime"

  } else if( covname %in% c("cov_matern_spacetime", "cov_matern_spacetime_cpp") ) {

    covmodel  <- cov_matern_spacetime
    covftn    <- "cov_matern_spacetime"

  } else if( covname %in% c("cov_matern_scaledim", "cov_matern_scaledim_cpp") ) {

    covmodel  <- cov_matern_scaledim
    covftn    <- "cov_matern_scaledim"

  } else {

    stop("This function only works for cov_expo_spacetime, cov_matern_spacetime, and cov_expo_scaledim for now.")
  }

  ### body
  time.tot <- proc.time()

  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")

  `%do%` <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`

  ### cand.all
  cand.all              <- data.frame(cand.m)
  cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
  colnames(cand.all)    <- c("index", "m")

  candid                <- data.frame(index = 1:nsim, nsim = 1)

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
  sim <- foreach::foreach(nsim = candid$nsim, .packages = c("correlationVecchia")) %dopar% simulate_spacetime_prediction(nsim = nsim, n = n, n.pred = n.pred, d = d, t = t, m = cand.all$m, nuggets = nuggets, abs.corr = abs.corr, method.locs = method.locs, method.locs.pred = method.locs.pred, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covftn, covparms = covparms)
  parallel::stopCluster(cl)

  ### return
  result                <- list()
  result[[1]]           <- sim[[1]]$setting
  result[[1]]$m         <- cand.m

  result[[2]]           <- list()
  result[[2]]           <- sim[[1]]$simout
  if(length(sim) > 1) {

    for(i in 2:length(sim)) {

      result[[2]]$b1$mspe <- rbind(result[[2]]$b1$mspe, sim[[i]]$simout$b1$mspe)
      result[[2]]$b2$mspe <- rbind(result[[2]]$b2$mspe, sim[[i]]$simout$b2$mspe)
      result[[2]]$b3$mspe <- rbind(result[[2]]$b3$mspe, sim[[i]]$simout$b3$mspe)
      result[[2]]$cc$mspe <- rbind(result[[2]]$cc$mspe, sim[[i]]$simout$cc$mspe)
      result[[2]]$eu$mspe <- rbind(result[[2]]$eu$mspe, sim[[i]]$simout$eu$mspe)

      result[[2]]$b1$logscore <- rbind(result[[2]]$b1$logscore, sim[[i]]$simout$b1$logscore)
      result[[2]]$b2$logscore <- rbind(result[[2]]$b2$logscore, sim[[i]]$simout$b2$logscore)
      result[[2]]$b3$logscore <- rbind(result[[2]]$b3$logscore, sim[[i]]$simout$b3$logscore)
      result[[2]]$cc$logscore <- rbind(result[[2]]$cc$logscore, sim[[i]]$simout$cc$logscore)
      result[[2]]$eu$logscore <- rbind(result[[2]]$eu$logscore, sim[[i]]$simout$eu$logscore)
    }
  }

  mspe.b1               <- colMeans(result[[2]]$b1$mspe)
  mspe.b2               <- colMeans(result[[2]]$b2$mspe)
  mspe.b3               <- colMeans(result[[2]]$b3$mspe)
  mspe.cc               <- colMeans(result[[2]]$cc$mspe)
  mspe.eu               <- colMeans(result[[2]]$eu$mspe)

  logs.b1               <- colMeans(result[[2]]$b1$logscore)
  logs.b2               <- colMeans(result[[2]]$b2$logscore)
  logs.b3               <- colMeans(result[[2]]$b3$logscore)
  logs.cc               <- colMeans(result[[2]]$cc$logscore)
  logs.eu               <- colMeans(result[[2]]$eu$logscore)

  result[[3]]           <- list(mspe = list(b1 = mspe.b1, b2 = mspe.b2, b3 = mspe.b3, cc = mspe.cc, eu = mspe.eu), logscore = list(b1 = logs.b1, b2 = logs.b2, b3 = logs.b3, cc = logs.cc, eu = logs.eu))

  ### return
  time.tot              <- proc.time() - time.tot
  result[[4]]           <- time.tot

  names(result) <- c("setting", "simout", "output", "time.tot")

  return( result )
}

#' @title parSim_sptm_posterior
#'
#' @param cand.m 1
#' @param target 1
#' @param ordfix 1
#' @param n 1
#' @param d 1
#' @param t 1
#' @param nuggets 1
#' @param method.locs 1
#' @param N 1
#' @param xlim 1
#' @param sdlog 1
#' @param method.modify 1
#' @param pivot 1
#' @param tol 1
#' @param verbose 1
#' @param tol.dec tol.dec
#' @param covparms 1
#' @param ncores 1
#'
#' @return 1
#'
#' @export
#'
#' @examples
#' 1 + 1
parSim_sptm_posterior <- function(cand.m, target, ordfix, n, d, t, nuggets = 0, method.locs, N = 100, xlim = c(0.05, 0.15), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = 4, covparms, ncores = NULL)
{
  time.tot  <- proc.time()

  # import::from(foreach,"%do%")
  # import::from(foreach,"%dopar%")

  `%do%`    <- foreach::`%do%`
  `%dopar%` <- foreach::`%dopar%`

  ### ncores
  if(is.null(ncores)) {
    no_cores  <- parallel::detectCores() - 2
  } else {
    no_cores  <- ncores
  }

  ### generate locs and z
  if( method.locs %in% c("random", "all.random") ) {

    process         <- generate_gp_spacetime(nsim = 1, n = n, d = d, t.len = t, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)

  } else if( method.locs %in% c("monitoring", "space.random.time.grid") ) {

    process         <- generate_gp_spacetime(nsim = 1, n = n, d = d, t.len = t, method.locs = "space.random.time.grid", covmodel = cov_matern_spacetime, covparms = covparms)

  } else if( method.locs %in% c("satellite") ) {

    process         <- generate_gp_spacetime(nsim = 1, n = n, d = d, t.len = t, method.locs = "satellite", covmodel = cov_matern_spacetime, covparms = covparms)

  } else {

    stop("Please check the argument method.locs!")
  }

  locs            <- process$sim$sim1$locs
  z               <- process$sim$sim1$y + sqrt(nuggets) * rnorm(nrow(locs))
  # fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")

  rm(process)

  ### candid
  approxs         <- c("b1", "b2", "b3", "cc", "ex") ; # approxs <- c("b1", "b2", "b3", "cc", "ex", "tr")
  candid          <- expand.grid(cand.m, approxs) ; colnames(candid) <- c("m", "approx")
  candid$approx   <- as.character(candid$approx)

  ### body
  cl              <- parallel::makeCluster(no_cores)

  if( is.null(ordfix) ) {

    doParallel::registerDoParallel(cl)
    simout <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia")) %dopar% posterior_ordfix_spacetime(target = target, approx = candid$approx[i], method.U = "simple", z = z, locs = locs, covparms = covparms, nugget = nuggets, m = candid$m[i], N = N, xlim = xlim, sdlog = sdlog, verbose = FALSE, tol.dec = tol.dec)

    simout2 <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia")) %dopar% posterior_ordfix_spacetime(target = target, approx = candid$approx[i], method.U = "ic0", z = z, locs = locs, covparms = covparms, nugget = nuggets, m = candid$m[i], N = N, xlim = xlim, sdlog = sdlog, verbose = FALSE, tol.dec = tol.dec)

    simout3 <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia")) %dopar% posterior_spacetime(target = target, approx = candid$approx[i], method.U = "simple", z = z, locs = locs, covparms = covparms, nugget = nuggets, m = candid$m[i], N = N, xlim = xlim, sdlog = sdlog, verbose = FALSE)

    simout4 <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia")) %dopar% posterior_spacetime(target = target, approx = candid$approx[i], method.U = "ic0", z = z, locs = locs, covparms = covparms, nugget = nuggets, m = candid$m[i], N = N, xlim = xlim, sdlog = sdlog, verbose = FALSE)
    parallel::stopCluster(cl)

  } else if( ordfix ) {

    doParallel::registerDoParallel(cl)
    simout <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia")) %dopar% posterior_ordfix_spacetime(target = target, approx = candid$approx[i], method.U = "simple", z = z, locs = locs, covparms = covparms, nugget = nuggets, m = candid$m[i], N = N, xlim = xlim, sdlog = sdlog, verbose = FALSE, tol.dec = tol.dec)

    simout2 <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia")) %dopar% posterior_ordfix_spacetime(target = target, approx = candid$approx[i], method.U = "ic0", z = z, locs = locs, covparms = covparms, nugget = nuggets, m = candid$m[i], N = N, xlim = xlim, sdlog = sdlog, verbose = FALSE, tol.dec = tol.dec)
    parallel::stopCluster(cl)

  } else {

    doParallel::registerDoParallel(cl)
    simout <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia")) %dopar% posterior_spacetime(target = target, approx = candid$approx[i], method.U = "simple", z = z, locs = locs, covparms = covparms, nugget = nuggets, m = candid$m[i], N = N, xlim = xlim, sdlog = sdlog, verbose = FALSE)

    simout2 <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia")) %dopar% posterior_spacetime(target = target, approx = candid$approx[i], method.U = "ic0", z = z, locs = locs, covparms = covparms, nugget = nuggets, m = candid$m[i], N = N, xlim = xlim, sdlog = sdlog, verbose = FALSE)
    parallel::stopCluster(cl)
  }

  ### return
  if( is.null(ordfix) ) {

    result          <- list()
    result[[1]]     <- list(candid = candid, target = target, ordfix = TRUE, n = n, d = d, t = t, nuggets = nuggets, method.U = "simple", method.locs = method.locs, N = N, xlim = xlim, sdlog = sdlog, covparms = covparms, locs = locs, z = z)
    result[[2]]     <- simout

    result[[3]]     <- list(candid = candid, target = target, ordfix = TRUE, n = n, d = d, t = t, nuggets = nuggets, method.U = "ic0", method.locs = method.locs, N = N, xlim = xlim, sdlog = sdlog, covparms = covparms, locs = locs, z = z)
    result[[4]]     <- simout2

    result[[5]]     <- list(candid = candid, target = target, ordfix = FALSE, n = n, d = d, t = t, nuggets = nuggets, method.U = "simple", method.locs = method.locs, N = N, xlim = xlim, sdlog = sdlog, covparms = covparms, locs = locs, z = z)
    result[[6]]     <- simout3

    result[[7]]     <- list(candid = candid, target = target, ordfix = FALSE, n = n, d = d, t = t, nuggets = nuggets, method.U = "ic0", method.locs = method.locs, N = N, xlim = xlim, sdlog = sdlog, covparms = covparms, locs = locs, z = z)
    result[[8]]     <- simout4

    time.tot        <- proc.time() - time.tot
    result[[9]]     <- time.tot

    names(result)   <- c("setting.ordfix", "simout.ordfix", "setting.ordfix.ic0", "simout.ordfix.ic0", "setting.notfix", "simout.notfix", "setting.notfix.ic0", "simout.notfix.ic0", "time.tot")

  } else {

    result          <- list()
    result[[1]]     <- list(candid = candid, target = target, ordfix = ordfix, n = n, d = d, t = t, nuggets = nuggets, method.U = "simple", method.locs = method.locs, N = N, xlim = xlim, sdlog = sdlog, covparms = covparms, locs = locs, z = z)
    result[[2]]     <- simout

    result[[3]]     <- list(candid = candid, target = target, ordfix = ordfix, n = n, d = d, t = t, nuggets = nuggets, method.U = "ic0", method.locs = method.locs, N = N, xlim = xlim, sdlog = sdlog, covparms = covparms, locs = locs, z = z)
    result[[4]]     <- simout2

    time.tot        <- proc.time() - time.tot
    result[[5]]     <- time.tot

    names(result)   <- c("setting", "simout", "setting.ic0", "simout.ic0", "time.tot")

  }

  return( result )
}
