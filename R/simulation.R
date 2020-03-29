####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to conduct simulation.
###
###   Contents:
###
####################################################################################



#' @title Conducting simulation on vecchia approximations for univariate GPs with known covariance parameters
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of locations
#' @param d A dimension of domain
#' @param m A size of conditioning sets
#' @param method.locs random or grid 
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param covmodel Covariance function
#' @param ... Covariance parameters
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' out <- simulate_univariate_knownCovparms(nsim = 2, n = 10^2, d = 2, m = 10, 
#'                                          method.locs = 'random', method.modify = NULL, 
#'                                          pivot = FALSE, tol = .Machine$double.eps, 
#'                                          verbose = TRUE, 
#'                                          covmodel = cov_expo_aniso, covparms = c(1, 0.1, 10))
#' out$kls.average
simulate_univariate_knownCovparms <- function(nsim, n, d, m, method.locs, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel, ...)
{
  time.tot      <- proc.time()
  time.sim      <- list()
  time.gen      <- proc.time()
  
  ### generation
  if(verbose == TRUE) {
    message("")
    message("------------------------------------------------------------")
    message(paste0("-----                    ", "generation", "                    -----"))
    message("------------------------------------------------------------")
  }
  
  if(verbose == TRUE) message(paste0("System: Simulation starts. [", Sys.time(), "]"))
  
  realization   <- generate_gp_space(nsim = nsim, n = n, d = d, p = 1, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, covmodel = covmodel, ...)
  
  if(verbose == TRUE) message(paste0("System: Realizations of GPs are generated. [", Sys.time(), "]"))
  
  time.gen      <- proc.time() - time.gen
  
  ### simulation
  kls.full      <- matrix(NA, nsim, 20)
  decomp.err    <- rep(NA, nsim)
  for(k in 1:nsim) {
    
    time.sim[[k]] <- proc.time()
    
    if(verbose == TRUE) {
      message("")
      message("------------------------------------------------------------")
      message(paste0("-----                   ", "simulation ", k, "                   -----"))
      message("------------------------------------------------------------")
    }
    
    ## basic information
    locs              <- realization$sim[[k]]$locs
    covmat.modified   <- realization$sim[[k]]$covmat
    decomp.err[k]     <- realization$sim[[k]]$covfac$decomp.err
    
    ## specify vecchia approximations
    approx        <- list()
    
    approx[[1]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: E-Maxmin + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
    approx[[2]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: E-Maxmin + C-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
    approx[[3]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-Maxmin + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
    approx[[4]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-Maxmin + C-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
    
    approx[[5]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "coord", ordering.method = "euclidean", coordinate = c(1), abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: X-coord + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
    approx[[6]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "coord", ordering.method = "euclidean", coordinate = c(2), abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: Y-coord + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
    
    ## compute approximate covariance matrices
    n.approx    <- length(approx)
    covmat.hat  <- list()
    kls         <- rep(NA, n.approx)
    for(i in 1:n.approx) {
      
      covmat.ord.modified   <- covmat.modified[approx[[i]]$ord, approx[[i]]$ord]
      U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified)$U
      
      # # previous version (wrong but noticeable)
      # covmat.ord            <- covmodel(approx[[i]]$locsord, ...)
      # covmat.ord.modified   <- modify(covmat.ord, pivot = pivot, method = method, tol = tol, return.err = TRUE, verbose = verbose)
      # U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified$covmat.modified)$U
      
      revord                <- order(approx[[i]]$ord)
      covmat.hat[[i]]       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
      
      kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]])
    }
    
    kls.full[k, 1:n.approx] <- kls
    
    if(verbose == TRUE) message(paste0("System: All vecchia approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
    
    time.sim[[k]] <- proc.time() - time.sim[[k]]
  }
  
  kls.full <- data.frame(kls.full[, 1:n.approx, drop = FALSE])
  colnames(kls.full) <- paste0("approx_", 1:n.approx)
  
  kls.average <- colMeans(kls.full)
  names(kls.average) <- paste0("approx_", 1:n.approx)
  
  time.tot  <- proc.time() - time.tot
  names(time.sim) <- paste0("simulation_", 1:nsim)
  
  ### return
  result                <- list()
  
  result$nsim           <- nsim
  result$n              <- n
  result$d              <- d
  result$m              <- m
  result$method.locs    <- method.locs
  result$method.modify  <- method.modify
  result$pivot          <- pivot
  result$tol            <- tol
  
  result$process        <- realization
  result$n.approx       <- n.approx
  result$kls.average    <- kls.average
  result$kls.full       <- kls.full
  result$decomp.error   <- decomp.err
  result$time.tot       <- time.tot
  result$time.gen       <- time.gen
  result$time.sim       <- time.sim
  
  if(verbose == TRUE) {
    message("")
    message(paste0("System: Simulation is finished. [", Sys.time(), "]"))
  }
  
  return(result)
}



#' @title Conducting simulation on vecchia approximations for multivariate GPs with known covariance parameters
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of locations
#' @param d A dimension of domain
#' @param p A number of GPs
#' @param m A size of conditioning sets
#' @param method.locs "random", "overlap", or "grid" 
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param covmodel Covariance function
#' @param ... Covariance parameters
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' out <- simulate_multivariate_knownCovparms(nsim = 2, n = 10^2, d = 2, p = 2, m = 10,
#'                                            method.locs = 'random', method.modify = NULL,
#'                                            pivot = FALSE, tol = .Machine$double.eps,
#'                                            verbose = TRUE,
#'                                            covmodel = cov_bivariate_expo_latDim, 
#'                                            covparms = c(1, 0.1, 1))
#' out$kls.average
simulate_multivariate_knownCovparms <- function(nsim, n, d, p, m, method.locs, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel, ...)
{
  time.tot      <- proc.time()
  time.sim      <- list()
  time.gen      <- proc.time()
  
  ### generation
  if(verbose == TRUE) {
    message("")
    message("------------------------------------------------------------")
    message(paste0("-----                    ", "generation", "                    -----"))
    message("------------------------------------------------------------")
  }
  
  if(verbose == TRUE) message(paste0("System: Simulation starts. [", Sys.time(), "]"))
  
  realization   <- generate_gp_space(nsim = nsim, n = n, d = d, p = p, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, covmodel = covmodel, ...)
  
  if(verbose == TRUE) message(paste0("System: Realizations of GPs are generated. [", Sys.time(), "]"))
  
  time.gen      <- proc.time() - time.gen
  
  ### simulation
  kls.full      <- matrix(NA, nsim, 20)
  decomp.err    <- rep(NA, nsim)
  for(k in 1:nsim) {
    
    time.sim[[k]] <- proc.time()
    
    if(verbose == TRUE) {
      message("")
      message("------------------------------------------------------------")
      message(paste0("-----                   ", "simulation ", k, "                   -----"))
      message("------------------------------------------------------------")
    }
    
    ## basic information
    locs              <- do.call(rbind, realization$sim[[k]]$locs)
    covmat.modified   <- realization$sim[[k]]$covmat
    decomp.err[k]     <- realization$sim[[k]]$covfac$decomp.err
    
    ## specify vecchia approximations
    approx        <- list()
    
    approx[[1]]   <- baseline_1_multivariate_specify(locs = realization$sim[[k]]$locs, m = m)
    if(verbose == TRUE) message(paste0("System: The first baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[2]]   <- baseline_2_multivariate_specify(locs = realization$sim[[k]]$locs, m = m)
    if(verbose == TRUE) message(paste0("System: The second baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[3]]   <- baseline_3_multivariate_specify(locs = realization$sim[[k]]$locs, m = m)
    if(verbose == TRUE) message(paste0("System: The third baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[4]]   <- baseline_4_multivariate_specify(locs = realization$sim[[k]]$locs, m = m, covmodel = covmodel, ...)
    if(verbose == TRUE) message(paste0("System: The fourth baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[5]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-Maxmin + C-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
    
    ## compute approximate covariance matrices
    n.approx    <- length(approx)
    covmat.hat  <- list()
    kls         <- rep(NA, n.approx)
    for(i in 1:n.approx) {
      
      locsord <- list()
      for(j in 1:p) {
        locsord[[j]] <- approx[[i]]$locsord[seq(from = 1 + (j-1) * n, to = j * n, by = 1), , drop = FALSE]
      }
      
      covmat.ord.modified   <- covmat.modified[approx[[i]]$ord, approx[[i]]$ord]
      U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified)$U
      
      # # previous version (wrong but noticeable)
      # covmat.ord            <- covmodel(locsord, ...)
      # covmat.ord.modified   <- modify(covmat.ord, pivot = pivot, method = method, tol = tol, return.err = TRUE, verbose = verbose)
      # U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified$covmat.modified)$U
      
      revord                <- order(approx[[i]]$ord)
      covmat.hat[[i]]       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
      
      kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]])
    }
    
    kls.full[k, 1:n.approx] <- kls
    
    if(verbose == TRUE) message(paste0("System: All vecchia approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
    
    time.sim[[k]] <- proc.time() - time.sim[[k]]
  }
  
  kls.full <- data.frame(kls.full[, 1:n.approx, drop = FALSE])
  colnames(kls.full) <- paste0("approx_", 1:n.approx)
  
  kls.average <- colMeans(kls.full)
  names(kls.average) <- paste0("approx_", 1:n.approx)
  
  time.tot  <- proc.time() - time.tot
  names(time.sim) <- paste0("simulation_", 1:nsim)
  
  ### return
  result                <- list()
  
  result$nsim           <- nsim
  result$n              <- n
  result$d              <- d
  result$p              <- p
  result$m              <- m
  result$method.locs    <- method.locs
  result$method.modify  <- method.modify
  result$pivot          <- pivot
  result$tol            <- tol
  
  result$process        <- realization
  result$n.approx       <- n.approx
  result$kls.average    <- kls.average
  result$kls.full       <- kls.full
  result$decomp.error   <- decomp.err
  result$time.tot       <- time.tot
  result$time.gen       <- time.gen
  result$time.sim       <- time.sim
  
  if(verbose == TRUE) {
    message("")
    message(paste0("System: Simulation is finished. [", Sys.time(), "]"))
  }
  
  return(result)
}



#' @title Conducting simulation on vecchia approximations for spacetime GPs with known covariance parameters
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of locations
#' @param d A dimension of domain
#' @param t.len A number of repeated measurement (= a number of different temporal locations)
#' @param m A size of conditioning sets
#' @param method.locs "all.random", "space.random.time.grid", "all.grid", or "satellite"
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param covmodel Covariance function
#' @param ... Covariance parameters 
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' out <- simulate_spacetime_knownCovparms(nsim = 2, n = 10^2, d = 2, 
#'                                         t.len = 2, m = 5, 
#'                                         method.locs = 'all.random', 
#'                                         covmodel = cov_spacetime_expo, 
#'                                         covparms = c(1, 0.75, 50, 25))
#' out$kls.average
simulate_spacetime_knownCovparms <- function(nsim, n, d, t.len, m, method.locs, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel, ...)
{
  time.tot      <- proc.time()
  time.sim      <- list()
  time.gen      <- proc.time()
  
  ### generation
  if(verbose == TRUE) {
    message("")
    message("------------------------------------------------------------")
    message(paste0("-----                    ", "generation", "                    -----"))
    message("------------------------------------------------------------")
  }
  
  if(verbose == TRUE) message(paste0("System: Simulation starts. [", Sys.time(), "]"))
  
  realization   <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t.len, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, covmodel = covmodel, ...)
  
  if(verbose == TRUE) message(paste0("System: Realizations of GPs are generated. [", Sys.time(), "]"))
  
  time.gen      <- proc.time() - time.gen
  
  ### simulation
  kls.full      <- matrix(NA, nsim, 20)
  decomp.err    <- rep(NA, nsim)
  for(k in 1:nsim) {
    
    time.sim[[k]] <- proc.time()
    
    if(verbose == TRUE) {
      message("")
      message("------------------------------------------------------------")
      message(paste0("-----                   ", "simulation ", k, "                   -----"))
      message("------------------------------------------------------------")
    }
    
    ## basic information
    locs              <- realization$sim[[k]]$locs
    covmat.modified   <- realization$sim[[k]]$covmat
    decomp.err[k]     <- realization$sim[[k]]$covfac$decomp.err
    
    ## specify vecchia approximations
    approx        <- list()
    
    approx[[1]]   <- baseline_1_spacetime_specify(locs = locs, m = m, coordinate = NULL)
    if(verbose == TRUE) message(paste0("System: The first baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[2]]   <- baseline_2_spacetime_specify(locs = locs, m = m, coordinate = NULL, theta = 1)
    if(verbose == TRUE) message(paste0("System: The second baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[3]]   <- baseline_3_spacetime_specify(locs = locs, m = m, coordinate = NULL, covmodel = covmodel, ...)
    if(verbose == TRUE) message(paste0("System: The third baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[4]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-Maxmin + C-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
    
    ## compute approximate covariance matrices
    n.approx    <- length(approx)
    covmat.hat  <- list()
    kls         <- rep(NA, n.approx)
    for(i in 1:n.approx) {
      
      covmat.ord.modified   <- covmat.modified[approx[[i]]$ord, approx[[i]]$ord]
      U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified)$U
      
      # # previous version (wrong but noticeable)
      # covmat.ord            <- covmodel(approx[[i]]$locsord, ...)
      # covmat.ord.modified   <- modify(covmat.ord, pivot = pivot, method = method, tol = tol, return.err = TRUE, verbose = verbose)
      # U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified$covmat.modified)$U
      
      revord                <- order(approx[[i]]$ord)
      covmat.hat[[i]]       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
      
      kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]])
    }
    
    kls.full[k, 1:n.approx] <- kls
    
    if(verbose == TRUE) message(paste0("System: All vecchia approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
    
    time.sim[[k]] <- proc.time() - time.sim[[k]]
  }
  
  kls.full <- data.frame(kls.full[, 1:n.approx, drop = FALSE])
  colnames(kls.full) <- paste0("approx_", 1:n.approx)
  
  kls.average <- colMeans(kls.full)
  names(kls.average) <- paste0("approx_", 1:n.approx)
  
  time.tot  <- proc.time() - time.tot
  names(time.sim) <- paste0("simulation_", 1:nsim)
  
  ### return
  result                <- list()
  
  result$nsim           <- nsim
  result$n              <- n
  result$d              <- d
  result$t.len          <- t.len
  result$m              <- m
  result$method.locs    <- method.locs
  result$method.modify  <- method.modify
  result$pivot          <- pivot
  result$tol            <- tol
  
  result$process        <- realization
  result$n.approx       <- n.approx
  result$kls.average    <- kls.average
  result$kls.full       <- kls.full
  result$decomp.error   <- decomp.err
  result$time.tot       <- time.tot
  result$time.gen       <- time.gen
  result$time.sim       <- time.sim
  
  if(verbose == TRUE) {
    message("")
    message(paste0("System: Simulation is finished. [", Sys.time(), "]"))
  }
  
  return(result)
}



#' @title Conducting simulation on vecchia approximations for (GP, dGP/dx, dGP/dy) with known covariance parameters
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of locations
#' @param d A dimension of domain. It must be 2 for now
#' @param m A size of conditioning sets
#' @param method.locs random or grid 
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param covmodel Covariance function
#' @param ... Covariance parameters 
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' par(mfrow = c(2, 1))
#' 
#' locs <- matrix(runif(10), 5, 2)
#' covmat <- cov_derivative_matern_2.5_2d(locs = locs, covparms = c(1, 0.1))
#' fields::image.plot(1 - covmat / max(diag(covmat)), main = "distance matrix with range = 0.1")
#' 
#' locs <- matrix(runif(10), 5, 2)
#' covmat <- cov_derivative_matern_2.5_2d(locs = locs, covparms = c(1, 4))
#' fields::image.plot(1 - covmat / max(diag(covmat)), main = "distance matrix with range = 4")
#' 
#' par(mfrow = c(1, 1))
#' 
#' out <- simulate_derivative_knownCovparms(nsim = 1, n = 10^2, d = 2, m = 10, 
#'                                          method.locs = 'random', 
#'                                          covmodel = cov_derivative_matern_2.5_2d, 
#'                                          covparms = c(1, 4), pivot = FALSE, 
#'                                          method.modify = "eigen-I", tol = 1e-6)
#' 
#' out$kls.average
simulate_derivative_knownCovparms <- function(nsim, n, d, m, method.locs, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel, ...)
{
  if(d != 2) stop("The dimension of the domain must be 2 for now.")
  
  time.tot      <- proc.time()
  time.sim      <- list()
  time.gen      <- proc.time()
  
  ### generation
  if(verbose == TRUE) {
    message("")
    message("------------------------------------------------------------")
    message(paste0("-----                    ", "generation", "                    -----"))
    message("------------------------------------------------------------")
  }
  
  if(verbose == TRUE) message(paste0("System: Simulation starts. [", Sys.time(), "]"))
  
  realization   <- generate_gp_space(nsim = nsim, n = n, d = d, p = 1, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, covmodel = cov_matern_2.5, ...)
  
  if(verbose == TRUE) message(paste0("System: Realizations of GPs are generated. [", Sys.time(), "]"))
  
  time.gen      <- proc.time() - time.gen
  
  ### simulation
  kls.full      <- matrix(NA, nsim, 20)
  decomp.err    <- rep(NA, nsim)
  for(k in 1:nsim) {
    
    time.sim[[k]] <- proc.time()
    
    if(verbose == TRUE) {
      message("")
      message("------------------------------------------------------------")
      message(paste0("-----                   ", "simulation ", k, "                   -----"))
      message("------------------------------------------------------------")
    }
    
    ## basic information
    locs              <- list(locs1 = realization$sim[[k]]$locs, locs2 = realization$sim[[k]]$locs, locs3 = realization$sim[[k]]$locs)
    locs.all          <- do.call(rbind, locs)
    covmat.modified   <- covmodel(locs = realization$sim[[k]]$locs, ...)
    covmat.modified   <- modify(covmat.modified, pivot = pivot, method = method.modify, tol = tol, return.err = TRUE, verbose = verbose)
    covmat.modified   <- covmat.modified$covmat.modified
      
    ## specify vecchia approximations
    approx        <- list()

    approx[[1]] <- baseline_1_multivariate_specify(locs = locs, m = m)
    if(verbose == TRUE) message(paste0("System: The first baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[2]] <- baseline_2_multivariate_specify(locs = locs, m = m)
    if(verbose == TRUE) message(paste0("System: The second baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[3]] <- baseline_3_multivariate_specify(locs = locs, m = m)
    if(verbose == TRUE) message(paste0("System: The third baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[4]] <- baseline_4_multivariate_specify(locs = locs, m = m, covmodel = covmodel, ...)
    if(verbose == TRUE) message(paste0("System: The fourth baseline approximation is accomplished. [", Sys.time(), "]"))
    approx[[5]] <- corrvecchia_specify_knownCovparms(locs = locs.all, m = m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat.modified, covparms = max(diag(covmat.modified)))
    if(verbose == TRUE) message(paste0("System: C-Maxmin + C-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
    
    ## compute approximate covariance matrices
    n.approx    <- length(approx)
    covmat.hat  <- list()
    kls         <- rep(NA, n.approx)
    for(i in 1:n.approx) {
      
      covmat.ord.modified   <- covmat.modified[approx[[i]]$ord, approx[[i]]$ord]
      U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified)$U
      
      revord                <- order(approx[[i]]$ord)
      covmat.hat[[i]]       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
      
      kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]])
    }
    
    kls.full[k, 1:n.approx] <- kls
    
    if(verbose == TRUE) message(paste0("System: All vecchia approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
    
    time.sim[[k]] <- proc.time() - time.sim[[k]]
  }
  
  kls.full <- data.frame(kls.full[, 1:n.approx, drop = FALSE])
  colnames(kls.full) <- paste0("approx_", 1:n.approx)
  
  kls.average <- colMeans(kls.full)
  names(kls.average) <- paste0("approx_", 1:n.approx)
  
  time.tot  <- proc.time() - time.tot
  names(time.sim) <- paste0("simulation_", 1:nsim)
  
  ### return
  result                <- list()
  
  result$nsim           <- nsim
  result$n              <- n
  result$d              <- d
  result$m              <- m
  result$method.locs    <- method.locs
  result$method.modify  <- method.modify
  result$pivot          <- pivot
  result$tol            <- tol
  
  result$process        <- realization
  result$n.approx       <- n.approx
  result$kls.average    <- kls.average
  result$kls.full       <- kls.full
  
  result$time.tot       <- time.tot
  result$time.gen       <- time.gen
  result$time.sim       <- time.sim
  
  if(verbose == TRUE) {
    message("")
    message(paste0("System: Simulation is finished. [", Sys.time(), "]"))
  }
  
  return(result)
}



