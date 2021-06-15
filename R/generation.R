####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes R functions to generate spatial/spatio-temporal locations and corresponding observations of gaussian processes.
###
###   Contents:
###     generate_gp / generate_gp_space / generate_gp_spacetime
###
####################################################################################

#' @title Generating realizations of mean-zero Gaussian Process (GP) using a user-defined covariance function
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of locations
#' @param d Dimension of domain
#' @param domain Type of domain. It must be either "space" or "spacetime."
#' @param p A number of variates (GPs). At \code{NULL} by default
#' @param t A number of repeated measurement (= a number of different temporal locations). At \code{NULL} by default
#' @param covmodel A covariance function
#' @param method.locs If domain = space, then it must be random, overlap, or grid. If domain = spacetime, then it must be random, monitoring, grid, or satellite. 
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix.  At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param ... Covariance parameters
#'
#' @return List. Please check generate_gp_space() and generate_gp_spacetime().
#' 
#' @importFrom stats rnorm
#' @importFrom stats runif
#' 
#' @export
#'
#' @examples
#' out <- generate_gp(nsim = 1, n = 100, d = 2, domain = "space", 
#'                    p = 1, method.locs = "random", 
#'                    covmodel = cov_expo_iso, covparms = c(1, 0.1))
#' 
#' plot(as.data.frame(out$sim$sim1$locs), col = 'red', lwd = 2, main = 'Spatial locations')
#' 
#' out <- generate_gp(nsim = 1, n = 100, d = 2, domain = "spacetime", 
#'                    t = 1, method.locs = "random", 
#'                    covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0))
#' 
#' plot(as.data.frame(out$sim$sim1$locs[, 1:2]), col = 'red', lwd = 2, main = "Spatial locations")
#' plot(out$sim$sim1$locs[, 3], main = "Measurement time against index")
generate_gp <- function(nsim, n, d, domain, p = NULL, t = NULL, covmodel, method.locs, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, ...)
{
  if(missing(domain)) {
    
    stop("The argument domain must be provided!")
    
  } else if(domain %in% c("space", "s", "sp", "spatial", "spatio")) {
    
    if(is.null(p)) stop("The argument p must be provided.")
    
    if(method.locs %in% c("all.random", "random")) {
      
      output <- generate_gp_space(nsim = nsim, n = n , d = d, p = p, method.locs = "random", covmodel = covmodel, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, ...)
      
    } else if(method.locs %in% c("all.overlap", "overlap")) {
      
      output <- generate_gp_space(nsim = nsim, n = n , d = d, p = p, method.locs = "overlap", covmodel = covmodel, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, ...)
      
    } else if(method.locs %in% c("all.grid", "grid")) {
      
      output <- generate_gp_space(nsim = nsim, n = n , d = d, p = p, method.locs = "grid", covmodel = covmodel, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, ...)
      
    } else {
      
      stop("Please check the argument method.locs!")
    }
    
  } else if(domain %in% c("spacetime", "st", "sptm", "spatio-temporal")) {
    
    if(is.null(t)) stop("The argument t must be provided.")
    
    if(method.locs %in% c("all.random", "random")) {
      
      output <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t, method.locs = "all.random", covmodel = covmodel, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, ...)
      
    } else if(method.locs %in% c("space.random.time.grid", "monitoring")) {
      
      output <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t, method.locs = "space.random.time.grid", covmodel = covmodel, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, ...)
      
    } else if(method.locs %in% c("all.grid", "grid")) {
      
      output <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t, method.locs = "all.grid", covmodel = covmodel, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, ...)
      
    } else if(method.locs %in% c("satellite", "footprint", "remote")) {
      
      output <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t, method.locs = "satellite", covmodel = covmodel, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, ...)
      
    } else {
      
      stop("Please check the argument method.locs!")
    }
    
  } else {
    
    stop("Check the argument domain!")
  }
  
  names(output)[which(names(output) == "t.len")] = "t"
  return( output )
}

#' @title Generating realizations of spatial mean-zero Gaussian Process (GP) using a user-defined covariance function
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of locations
#' @param d A dimension of domain
#' @param p A number of GPs
#' @param method.locs random, overlap, or grid 
#' @param covmodel A covariance function
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix.  At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param ... Covariance parameters
#'
#' @return \code{generate_gp} returns
#'     \itemize{
#'         \item{\code{nsim}: } A number of realizations of the GP
#'         \item{\code{n}: } An actual number of locations
#'         \item{\code{d}: } A dimension of domain
#'         \item{\code{p}: } A number of GPs
#'         \item{\code{method.locs}: } The argument method.locs
#'         \item{\code{method.modify}: } The argument method.modify
#'         \item{\code{pivot}: } The argument pivot
#'         \item{\code{tol}: } The argument tol
#'         \item{\code{sim}: } Simulation results
#'     }
#'  
#' @importFrom stats rnorm
#' @importFrom stats runif
#'           
#' @export
#'
#' @examples
#' par(mfrow = c(2, 3))
#' 
#' out <- generate_gp_space(nsim = 2, n = 19, d = 2, p = 1, 
#'                          method.locs = 'random', covmodel = cov_expo_iso, 
#'                          method.modify = NULL, pivot = FALSE, 
#'                          tol = .Machine$double.eps, 
#'                          covparms = c(1, 0.1))
#' plot(as.data.frame(out$sim$sim1$locs), col = 'red', lwd = 2, 
#'      main = 'Random locs for univariate GP')
#' 
#' out <- generate_gp_space(nsim = 2, n = 19, d = 2, p = 1, 
#'                          method.locs = 'overlap', covmodel = cov_expo_iso, 
#'                          method.modify = NULL, pivot = FALSE, 
#'                          tol = .Machine$double.eps, 
#'                          covparms = c(1, 0.1))
#' plot(as.data.frame(out$sim$sim1$locs), col = 'red', lwd = 2, 
#'      main = 'Overlapped locs for univariate GP')
#' 
#' out <- generate_gp_space(nsim = 2, n = 19, d = 2, p = 1, 
#'                          method.locs = 'grid', covmodel = cov_expo_iso, 
#'                          method.modify = NULL, pivot = FALSE, 
#'                          tol = .Machine$double.eps, 
#'                          covparms = c(1, 0.1))
#' plot(as.data.frame(out$sim$sim1$locs), col = 'red', lwd = 2, 
#'      main = 'Grid locs for univariate GP')
#' 
#' out <- generate_gp_space(nsim = 2, n = 19, d = 2, p = 2, 
#'                          method.locs = 'random', covmodel = cov_flexMatern_biv, 
#'                          method.modify = NULL, pivot = FALSE, 
#'                          tol = .Machine$double.eps, 
#'                          sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                          nu.mat = matrix(0.5, 2, 2), 
#'                          alpha.mat = matrix(1, 2, 2))
#' plot(as.data.frame(out$sim$sim1$locs[[1]]), col = 'red', lwd = 2, 
#'      main = 'Random locs for bivariate GP')
#' points((as.data.frame(out$sim$sim1$locs[[2]])), col = 'blue')
#' 
#' out <- generate_gp_space(nsim = 2, n = 19, d = 2, p = 2, 
#'                          method.locs = 'overlap', covmodel = cov_flexMatern_biv, 
#'                          method.modify = NULL, pivot = FALSE, 
#'                          tol = .Machine$double.eps, 
#'                          sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                          nu.mat = matrix(0.5, 2, 2), 
#'                          alpha.mat = matrix(1, 2, 2))
#' plot(as.data.frame(out$sim$sim1$locs[[1]]), col = 'red', lwd = 2, 
#'      main = 'Overlapped locs for bivariate GP')
#' points((as.data.frame(out$sim$sim1$locs[[2]])), col = 'blue')
#' 
#' out <- generate_gp_space(nsim = 2, n = 19, d = 2, p = 2, 
#'                          method.locs = 'grid', covmodel = cov_flexMatern_biv, 
#'                          method.modify = NULL, pivot = FALSE, 
#'                          tol = .Machine$double.eps, 
#'                          sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                          nu.mat = matrix(0.5, 2, 2), 
#'                          alpha.mat = matrix(1, 2, 2))
#' plot(as.data.frame(out$sim$sim1$locs[[1]]), col = 'red', lwd = 2,
#'      main = 'Grid locs for bivariate GP')
#' points((as.data.frame(out$sim$sim1$locs[[2]])), col = 'blue')
#' 
#' par(mfrow = c(1, 1))
generate_gp_space <- function(nsim, n, d, p, method.locs, covmodel, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, ...)
{
  ### location
  if(method.locs == "random") {
    
    if(p == 1) {
      locs.full     <- matrix(stats::runif(nsim * n * d, 0, 1), nsim * n, d)
      noise.full    <- stats::rnorm(nsim * n)
    } else if(p > 1) {
      locs.full     <- matrix(stats::runif(nsim * n * p * d, 0, 1), nsim * n * p, d)
      noise.full    <- stats::rnorm(nsim * n * p)
    } else {
      stop("The argument p is invalid.")
    }
    
  } else if(method.locs == "overlap") {
    
    if(p == 1) {
      locs.full     <- matrix(stats::runif(nsim * n * d, 0, 1), nsim * n, d)
      noise.full    <- stats::rnorm(nsim * n)
    } else if(p > 1) {
      locs.full     <- matrix(stats::runif(nsim * n * d, 0, 1), nsim * n, d)
      ind           <- unlist(rep(split(seq(nsim * n), rep(seq(nsim), each = n)), each = p))
      locs.full     <- locs.full[ind, , drop = FALSE]
      noise.full    <- stats::rnorm(nsim * n * p)
    } else {
      stop("The argument p is invalid.")
    }
    
  } else if(method.locs == "grid") {
    
    if(p == 1) {
      ni            <- as.integer(n^(1/d)) ; if(ni^d != n) message(paste0("Note: the effective sample size is ", ni^d, " but not ", n, "."))
      n             <- ni^d
      
      locs.full     <- as.matrix(expand.grid(split(rep(seq(from = 0, to = 1, length.out = ni), each = d), seq(d))))
      locs.full     <- locs.full[rep(seq(nrow(locs.full)), times = nsim), , drop = FALSE]
      noise.full    <- stats::rnorm(nsim * n)
    } else if(p > 1) {
      ni            <- as.integer(n^(1/d)) ; if(ni^d != n) message(paste0("Note: the effective sample size is ", ni^d, " but not ", n, "."))
      n             <- ni^d
      
      locs.full     <- as.matrix(expand.grid(split(rep(seq(from = 0, to = 1, length.out = ni), each = d), seq(d))))
      locs.full     <- locs.full[rep(seq(nrow(locs.full)), times = nsim * p), , drop = FALSE]
      noise.full    <- stats::rnorm(nsim * n * p)
    } else {
      stop("The argument p is invalid.")
    }
    
  } else {
    stop("The argument method.locs is not valid.")
  }
  
  ### process
  sim <- list()
  if(p == 1) {
    
    for(k in 1:nsim) {
      ind           <- seq(from = 1 + (k - 1) * n, to = k * n, by = 1)
      locs          <- locs.full[ind, , drop = FALSE]
      
      covmat        <- covmodel(locs, ...)
      covfac        <- factorize(covmat = covmat, pivot = pivot, method = method.modify, tol = tol, return.err = TRUE, verbose = verbose)
      
      y             <- as.numeric(t(covfac$covfactor) %*% noise.full[ind])
      covmat        <- crossprod(covfac$covfactor)  
      
      # sim[[k]] <- ind
      sim[[k]]      <- list(locs = locs, y = y, covmat = covmat, covfac = covfac)
    }
    
  } else if(p > 1) {
    
    for(k in 1:nsim) {
      locs  <- list()
      for(i in 1:p) {
        ind           <- seq(from = 1 + (k - 1) * n * p + (i - 1) * n, to = (k - 1) * n * p + i * n, by = 1)
        
        # locs[[i]] <- ind
        locs[[i]]     <- locs.full[ind, , drop = FALSE]
      }
      names(locs)   <- paste0("locs", seq(p))
      
      covmat        <- covmodel(locs, ...)
      covfac        <- factorize(covmat = covmat, pivot = pivot, method = method.modify, tol = tol, return.err = TRUE, verbose = verbose)
      
      ind           <- seq(from = 1 + (k - 1) * n * p, to = k * n * p, by = 1)
      # y <- ind
      y             <- as.numeric(t(covfac$covfactor) %*% noise.full[ind])
      covmat        <- crossprod(covfac$covfactor)
      
      sim[[k]]      <- list(locs = locs, y = y, covmat = covmat, covfac = covfac)
    }
    
  } else {
    stop("The argument p is invalid.")
  }
  
  ### return
  names(sim)  <- paste0("sim", seq(nsim))
  return( list(nsim = nsim, n = n, d = d, p = p, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, sim = sim) )
}

#' @title Generating realizations of spatio-temporal mean-zero Gaussian Process (GP) using a user-defined covariance function
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of different spatial locations
#' @param d A dimension of space domain
#' @param t.len A number of repeated measurement (= a number of different temporal locations)
#' @param method.locs all.random, space.random.time.grid, all.grid, or satellite
#' @param covmodel A covariance function
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix.  At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param ... Covariance parameters
#' 
#' @return \code{generate_gp_spacetime} returns
#'     \itemize{
#'         \item{\code{nsim}: } A number of realizations of the GP
#'         \item{\code{size}: } A size of realizations of the GP
#'         \item{\code{n}: } An actual number of spatial locations
#'         \item{\code{d}: } A dimension of space domain
#'         \item{\code{t.len}: } A number of repeated measurement (= a number of different temporal locations)
#'         \item{\code{method.locs}: } The argument method.locs
#'         \item{\code{method.modify}: } The argument method.modify
#'         \item{\code{pivot}: } The argument pivot
#'         \item{\code{tol}: } The argument tol
#'         \item{\code{sim}: } Simulation results
#'     }
#'  
#' @importFrom stats rnorm
#' @importFrom stats runif
#' 
#' @export
#'
#' @examples
#' par(mfrow = c(4, 2))
#' 
#' out <- generate_gp_spacetime(nsim = 2, n = 300, d = 2, t.len = 3, 
#'                              method.locs = "all.random", 
#'                              covmodel = cov_expo_spacetime_nonsep, 
#'                              covparms = c(1, 0.75, 50, 25))
#' plot(as.data.frame(out$sim$sim1$locs[, 1:2]), main = "Spatial locations")
#' plot(out$sim$sim1$locs[, 3], main = "Measurement time against index")
#' 
#' out <- generate_gp_spacetime(nsim = 2, n = 300, d = 2, t.len = 3, 
#'                              method.locs = "space.random.time.grid", 
#'                              covmodel = cov_expo_spacetime_nonsep, 
#'                              covparms = c(1, 0.75, 50, 25))
#' plot(as.data.frame(out$sim$sim1$locs[, 1:2]), main = "Spatial locations")
#' plot(out$sim$sim1$locs[, 3], main = "Measurement time against index")
#' 
#' out <- generate_gp_spacetime(nsim = 2, n = 300, d = 2, t.len = 3, 
#'                              method.locs = "all.grid", 
#'                              covmodel = cov_expo_spacetime_nonsep, 
#'                              covparms = c(1, 0.75, 50, 25))
#' plot(as.data.frame(out$sim$sim1$locs[, 1:2]), main = "Spatial locations")
#' plot(out$sim$sim1$locs[, 3], main = "Measurement time against index")
#' 
#' out <- generate_gp_spacetime(nsim = 2, n = 300, d = 2, t.len = 3, 
#'                              method.locs = "satellite", 
#'                              covmodel = cov_expo_spacetime_nonsep, 
#'                              covparms = c(1, 0.75, 50, 25))
#' plot(as.data.frame(out$sim$sim1$locs[, 1:2]), main = "Spatial locations")
#' plot(out$sim$sim1$locs[, 3], main = "Measurement time against index")
#' 
#' par(mfrow = c(1, 1))
generate_gp_spacetime <- function(nsim, n, d, t.len, method.locs, covmodel, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, ...)
{
  ### location
  if(method.locs == "all.random") {
    
    message("Note: t.len (or t) is not used when method.locs is all.random (or random).")
    
    size          <- n
    locs.full     <- matrix(stats::runif(nsim * n * (d+1), 0, 1), nsim * n, d+1)
    noise.full    <- stats::rnorm(nsim * n)
    
  } else if(method.locs == "space.random.time.grid") {
    
    size          <- n * t.len
    locs.space    <- matrix(stats::runif(nsim * n * d, 0, 1), nsim * n, d)
    locs.time     <- seq(from = 0, to = 1, length.out = t.len)
    
    ind.space     <- rep(seq(nsim * n), each = t.len)
    ind.time      <- rep(seq(t.len), times = nsim * n)
    
    locs.full     <- cbind(locs.space[ind.space, ], locs.time[ind.time])
    noise.full    <- stats::rnorm(nsim * size)
    
  } else if(method.locs == "all.grid") {
    
    ni            <- as.integer(n^(1/d)) ; if(ni^d != n) message(paste0("Note: the effective sample size is ", ni^d, " but not ", n, "."))
    n             <- ni^d
    
    size          <- n * t.len
    locs.space    <- as.matrix(expand.grid(split(rep(seq(from = 0, to = 1, length.out = ni), each = d), seq(d))))
    locs.time     <- seq(from = 0, to = 1, length.out = t.len)
    
    ind.space     <- rep(seq(n), each = t.len)
    ind.time      <- rep(seq(t.len), times = n)
    
    locs.full     <- cbind(locs.space[ind.space, ], locs.time[ind.time])
    locs.full     <- locs.full[rep(seq(nrow(locs.full)), times = nsim), , drop = FALSE]
    noise.full    <- stats::rnorm(nsim * size)
    
  } else if(method.locs == "satellite") {
    
    message("Note: location matrix is hard-coded when method.locs is satellite. The number of different x coordinate values is 5 and satellite is assumed to measure each location twice.")
    
    nx <- 5 ; ny <- as.integer(n / nx / 2)
    if(nx * ny * 2 != n) message(paste0("Note: the effective sample size is ", nx * ny * 2, " but not ", n, "."))
    n             <- nx * ny * 2
    size          <- n 
    
    locs.space    <- cbind(rep(seq(from = 0.1, to = 0.9, length.out = nx), each = ny), rep(seq(from = 0.1, to = 0.9, length.out = ny), times = nx))
    locs.space    <- rbind(locs.space, locs.space)
    locs.time     <- seq(from = 0, to = 1, length.out = n)
    
    locs.full     <- cbind(locs.space, locs.time)
    locs.full     <- locs.full[rep(seq(nrow(locs.full)), times = nsim), , drop = FALSE]
    noise.full    <- stats::rnorm(nsim * n)
    
    # message("Note: location matrix is hard-coded when method.locs is satellite. The number of different x coordinate values is 10.")
    # 
    # nx <- 10 ; ny <- as.integer(n / nx)
    # if(nx * ny != n) message(paste0("Note: the effective sample size is ", nx * ny, " but not ", n, "."))
    # n             <- nx * ny
    # size          <- n
    # 
    # locs.space    <- cbind(rep(seq(from = 0.1, to = 0.9, length.out = nx), each = ny), rep(seq(from = 0.1, to = 0.9, length.out = ny), times = nx))
    # locs.time     <- seq(from = 0, to = 1, length.out = n)
    # 
    # locs.full     <- cbind(locs.space, locs.time)
    # locs.full     <- locs.full[rep(seq(nrow(locs.full)), times = nsim), , drop = FALSE]
    # noise.full    <- stats::rnorm(nsim * n)
    
  } else {
    stop("The argument method.locs is not valid.")
  }
  
  ### process
  sim <- list()
  for(k in 1:nsim) {
    ind           <- seq(from = 1 + (k - 1) * size, to = k * size, by = 1)
    locs          <- locs.full[ind, , drop = FALSE]
    
    covmat        <- covmodel(locs, ...)
    covfac        <- factorize(covmat = covmat, pivot = pivot, method = method.modify, tol = tol, return.err = TRUE, verbose = verbose)
    
    y             <- as.numeric(t(covfac$covfactor) %*% noise.full[ind])
    covmat        <- crossprod(covfac$covfactor)
    
    # sim[[k]] <- ind
    sim[[k]]      <- list(locs = locs, y = y, covmat = covmat, covfac = covfac)
  }
  
  ### return
  names(sim)  <- paste0("sim", seq(nsim))
  return( list(nsim = nsim, size = size, n = n, d = d, t.len = t.len, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, sim = sim) )
}