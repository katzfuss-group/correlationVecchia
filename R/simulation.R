####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to perform simulation.
###
###   Contents: 
###     simulate_univ / simulate_mulv / simulate_spacetime / simulate_deriv
###
####################################################################################

#' @title Conducting simulation on vecchia approximations for univariate GPs with known covariance parameters
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of locations
#' @param d A dimension of domain
#' @param m A size of conditioning sets
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
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param covmodel Covariance model
#' @param ... Covariance parameters
#'
#' @return Simulation output (list)
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' output <- simulate_univ(nsim = 2, n = 20^2, d = 2, m = 20, 
#'                         abs.corr = FALSE, method.locs = "random",
#'                         method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps,
#'                         verbose = TRUE, 
#'                         covmodel = cov_expo_aniso, covparms = c(1, 0.1, 10))
#' 
#' barplot(log10(output$kls.average), names.arg = output$names, 
#'         main = "log10-scale KL divergences")
#' }
simulate_univ <- function(nsim, n, d, m, abs.corr, method.locs, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel, ...)
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
  
  realization   <- generate_gp(nsim = nsim, n = n, d = d, domain = "space", p = 1, t = NULL, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, covmodel = covmodel, ...)
  
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
    
    approx[[1]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "euc", conditioning = "NN", conditioning.method = "euc", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: E-MM + E-NN is accomplished. [", Sys.time(), "]"))
    approx[[2]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "euc", conditioning = "NN", conditioning.method = "cor", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: E-MM + C-NN is accomplished. [", Sys.time(), "]"))
    approx[[3]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "euc", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-MM + E-NN is accomplished. [", Sys.time(), "]"))
    approx[[4]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-MM + C-NN is accomplished. [", Sys.time(), "]"))
    approx[[5]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = c(1), ordering = "coord", ordering.method = "euc", conditioning = "NN", conditioning.method = "euc", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: X-ord + E-NN is accomplished. [", Sys.time(), "]"))
    approx[[6]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = c(2), ordering = "coord", ordering.method = "euc", conditioning = "NN", conditioning.method = "euc", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: Y-ord + E-NN is accomplished. [", Sys.time(), "]"))
    
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
      
      kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]], method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps)
    }
    
    kls.full[k, 1:n.approx] <- kls
    
    if(verbose == TRUE) message(paste0("System: All the approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
    
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
  result$abs.corr       <- abs.corr
  result$method.locs    <- method.locs
  result$method.modify  <- method.modify
  result$pivot          <- pivot
  result$tol            <- tol
  
  result$names          <- c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN")
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

#' @title Conducting simulation on vecchia approximations for univariate GPs with known covariance parameters (without coord ordering)
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of locations
#' @param d A dimension of domain
#' @param m A size of conditioning sets
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
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param covmodel Covariance model
#' @param ... Covariance parameters
#'
#' @return Simulation output (list)
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' output <- simulate_univ_cpt(nsim = 2, n = 20^2, d = 2, m = 20, 
#'                         abs.corr = FALSE, method.locs = "random",
#'                         method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps,
#'                         verbose = TRUE, 
#'                         covmodel = cov_expo_aniso, covparms = c(1, 0.1, 10))
#' 
#' barplot(log10(output$kls.average), names.arg = output$names, 
#'         main = "log10-scale KL divergences")
#' }
simulate_univ_cpt <- function(nsim, n, d, m, abs.corr, method.locs, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel, ...)
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
  
  realization   <- generate_gp(nsim = nsim, n = n, d = d, domain = "space", p = 1, t = NULL, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, covmodel = covmodel, ...)
  
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
    
    approx[[1]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "euc", conditioning = "NN", conditioning.method = "euc", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: E-MM + E-NN is accomplished. [", Sys.time(), "]"))
    approx[[2]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "euc", conditioning = "NN", conditioning.method = "cor", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: E-MM + C-NN is accomplished. [", Sys.time(), "]"))
    approx[[3]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "euc", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-MM + E-NN is accomplished. [", Sys.time(), "]"))
    approx[[4]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-MM + C-NN is accomplished. [", Sys.time(), "]"))
    
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
      
      kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]], method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps)
    }
    
    kls.full[k, 1:n.approx] <- kls
    
    if(verbose == TRUE) message(paste0("System: All the approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
    
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
  result$abs.corr       <- abs.corr
  result$method.locs    <- method.locs
  result$method.modify  <- method.modify
  result$pivot          <- pivot
  result$tol            <- tol
  
  result$names          <- c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN")
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
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param covmodel Covariance function
#' @param ... Covariance parameters
#'
#' @return Simulation result (list)
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' output <- simulate_mulv(nsim = 2, n = 20^2, d = 2, p = 2, m = 20,
#'                         method.locs = 'random',
#'                         method.modify = NULL, abs.corr = FALSE,
#'                         pivot = FALSE, tol = .Machine$double.eps,
#'                         verbose = TRUE,
#'                         covmodel = cov_latentDim_biv,
#'                         covparms = c(1, 0.1, 0.1))
#' 
#' barplot(log10(output$kls.average), names.arg = output$names, 
#'         main = "log10-scale KL divergences")
#' }
simulate_mulv <- function(nsim, n, d, p, m, abs.corr, method.locs, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel, ...)
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
  
  realization   <- generate_gp(nsim = nsim, n = n, d = d, domain = "space", p = p, t = NULL, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, covmodel = covmodel, ...)
  
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
    covmat.modified   <- realization$sim[[k]]$covmat ; realization$sim[[k]]$covmat <- NULL
    decomp.err[k]     <- realization$sim[[k]]$covfac$decomp.err
    
    ## specify vecchia approximations
    approx        <- list()
    
    approx[[1]]   <- baseline_mulv_specify(approx = 1, locs = realization$sim[[k]]$locs, m = m)
    if(verbose == TRUE) message(paste0("System: S-E-MM + HH-E-NN is accomplished. [", Sys.time(), "]"))
    approx[[2]]   <- baseline_mulv_specify(approx = 2, locs = realization$sim[[k]]$locs, m = m)
    if(verbose == TRUE) message(paste0("System: S-E-MM + J-E-NN is accomplished. [", Sys.time(), "]"))
    approx[[3]]   <- baseline_mulv_specify(approx = 3, locs = realization$sim[[k]]$locs, m = m)
    if(verbose == TRUE) message(paste0("System: S-E-MM + S-E-NN is accomplished. [", Sys.time(), "]"))
    approx[[4]]   <- baseline_mulv_specify(approx = 4, locs = realization$sim[[k]]$locs, m = m, abs.corr = abs.corr, covmodel = covmodel, ...)
    if(verbose == TRUE) message(paste0("System: S-E-MM + C-NN is accomplished. [", Sys.time(), "]"))
    approx[[5]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-MM + C-NN is accomplished. [", Sys.time(), "]"))
    
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
      
      kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]], method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps)
    }
    
    kls.full[k, 1:n.approx] <- kls
    
    if(verbose == TRUE) message(paste0("System: All the approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
    
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
  result$abs.corr       <- abs.corr
  result$method.locs    <- method.locs
  result$method.modify  <- method.modify
  result$pivot          <- pivot
  result$tol            <- tol
  
  result$names          <- c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN")
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
#' @param t A number of repeated measurement (= a number of different temporal locations)
#' @param m A size of conditioning sets
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
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param covmodel Covariance function
#' @param ... Covariance parameters 
#'
#' @return Simulation result (list)
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' output <- simulate_spacetime(nsim = 2, n = 20^2, d = 2, t = 1, m = 20, abs.corr = FALSE,
#'                              method.locs = "random", method.modify = NULL, pivot = FALSE, 
#'                              tol = .Machine$double.eps, verbose = TRUE, 
#'                              covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0))
#' 
#' barplot(log10(output$kls.average), names.arg = output$names,
#'         main = "log10-scale KL divergences")
#' }
simulate_spacetime <- function(nsim, n, d, t, m, abs.corr, method.locs, method.modify = NULL, pivot, tol = .Machine$double.eps, verbose = TRUE, covmodel, ...)
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
  
  realization   <- generate_gp(nsim = nsim, n = n, d = d, domain = "spacetime", p = NULL, t = t, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, covmodel = covmodel, ...)
  
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
    
    approx[[1]]   <- baseline_spacetime_specify(approx = 1, locs = locs, m = m, coordinate = NULL)
    if(verbose == TRUE) message(paste0("System: T-ord + T-NN is accomplished. [", Sys.time(), "]"))
    approx[[2]]   <- baseline_spacetime_specify(approx = 2, locs = locs, m = m, coordinate = NULL)
    if(verbose == TRUE) message(paste0("System: T-ord + E-NN is accomplished. [", Sys.time(), "]"))
    approx[[3]]   <- baseline_spacetime_specify(approx = 3, locs = locs, m = m, coordinate = NULL, abs.corr = abs.corr, covmodel = covmodel, ...)
    if(verbose == TRUE) message(paste0("System: T-ord + C-NN is accomplished. [", Sys.time(), "]"))
    approx[[4]]   <- cvecchia_specify(locs = locs, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-MM + C-NN is accomplished. [", Sys.time(), "]"))
    
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
      
      kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]], method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps)
    }
    
    kls.full[k, 1:n.approx] <- kls
    
    if(verbose == TRUE) message(paste0("System: All the approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
    
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
  result$t              <- t
  result$m              <- m
  result$abs.corr       <- abs.corr
  result$method.locs    <- method.locs
  result$method.modify  <- method.modify
  result$pivot          <- pivot
  result$tol            <- tol
  
  result$names          <- c("T-ord + T-NN", "T-ord + T-NN", "T-ord + C-NN", "C-MM + C-NN")
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

#' @title Conducting simulation on vecchia approximations for (GP, dGP/dx1 dGP/dx2, ...) with known covariance parameters
#'
#' @param nsim A number of realizations of the GP
#' @param n A number of locations
#' @param d A dimension of domain. It must be 2 for now
#' @param m A size of conditioning sets
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
#' @param verbose Logical at \code{TRUE} by default. If verbose is \code{TRUE}, then this function prints out all messages
#' @param covmodel Covariance function
#' @param ... Covariance parameters 
#'
#' @return Simulation result (list)
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(05222021)
#' 
#' output <- simulate_deriv(nsim = 2, n = 20^2, d = 2, m = 20, abs.corr = FALSE,
#'                          method.locs = "random", method.modify = "eigen-I", pivot = FALSE,
#'                          tol = 1e-4, verbose = TRUE,
#'                          covmodel = cov_sqexpo_deriv, covparms = c(1, 0.1))
#' 
#' barplot(log10(output$kls.average), names.arg = output$names,
#'         main = "log10-scale KL divergences")
#' }
simulate_deriv <- function(nsim, n, d, m, abs.corr, method.locs, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel, ...)
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
  
  realization   <- generate_gp(nsim = nsim, n = n, d = d, domain = "space", p = 1, t = NULL, method.locs = method.locs, method.modify = method.modify, pivot = pivot, tol = tol, verbose = verbose, covmodel = cov_expo_iso, covparms = c(1, 0.1)) # CAUTION: This cov_expo_iso is not used.
  
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
    if(d == 1) locs <- list(locs1 = realization$sim[[k]]$locs, locs2 = realization$sim[[k]]$locs)
    if(d == 2) locs <- list(locs1 = realization$sim[[k]]$locs, locs2 = realization$sim[[k]]$locs, locs3 = realization$sim[[k]]$locs)
    if(d  > 2) stop("Not yet!")
    
    locs.all          <- do.call(rbind, locs)
    
    covmat.modified   <- covmodel(locs = realization$sim[[k]]$locs, ...)
    covmat.modified   <- modify(covmat.modified, pivot = pivot, method = method.modify, tol = tol, return.err = TRUE, verbose = verbose)
    covmat.modified   <- covmat.modified$covmat.modified
    
    ## specify vecchia approximations
    approx        <- list()
    
    approx[[1]] <- baseline_mulv_specify(approx = 1, locs = locs, m = m)
    if(verbose == TRUE) message(paste0("System: S-E-MM + HH-E-NN is accomplished. [", Sys.time(), "]"))
    approx[[2]] <- baseline_mulv_specify(approx = 2, locs = locs, m = m)
    if(verbose == TRUE) message(paste0("System: S-E-MM + J-E-NN is accomplished. [", Sys.time(), "]"))
    approx[[3]] <- baseline_mulv_specify(approx = 3, locs = locs, m = m)
    if(verbose == TRUE) message(paste0("System: S-E-MM + S-E-NN is accomplished. [", Sys.time(), "]"))
    approx[[4]] <- baseline_mulv_specify(approx = 4, locs = locs, m = m, abs.corr = abs.corr, covmodel = covmodel, ...)
    if(verbose == TRUE) message(paste0("System: S-E-MM + C-NN is accomplished. [", Sys.time(), "]"))
    approx[[5]] <- cvecchia_specify(locs = locs.all, m = m, rho = NULL, initial.pt = NULL, abs.corr = abs.corr, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor", covmodel = covmat.modified, covparms = c(1))
    if(verbose == TRUE) message(paste0("System: C-MM + C-NN is accomplished. [", Sys.time(), "]"))

    ## compute approximate covariance matrices
    n.approx    <- length(approx)
    covmat.hat  <- list()
    kls         <- rep(NA, n.approx)
    for(i in 1:n.approx) {
      
      covmat.ord.modified   <- covmat.modified[approx[[i]]$ord, approx[[i]]$ord]
      U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified)$U
      
      revord                <- order(approx[[i]]$ord)
      covmat.hat[[i]]       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
      
      kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]], method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps)
    }
    
    kls.full[k, 1:n.approx] <- kls
    
    if(verbose == TRUE) message(paste0("System: All the approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
    
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
  result$abs.corr       <- abs.corr
  result$method.locs    <- method.locs
  result$method.modify  <- method.modify
  result$pivot          <- pivot
  result$tol            <- tol
  
  result$names          <- c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN")
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

#' @title simulate_spacetime_prediction
#'
#' @param nsim nsim
#' @param n n
#' @param n.pred n.pred 
#' @param d d
#' @param t t 
#' @param m m
#' @param nuggets nuggets 
#' @param abs.corr logical
#' @param method.locs random, monitoring, or satellite
#' @param method.locs.pred random or subset
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix. 
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param verbose logical
#' @param covmodel covmodel
#' @param covparms covparms
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' 1 + 1
simulate_spacetime_prediction <- function(nsim, n, n.pred, d, t, m, nuggets = 0, abs.corr = FALSE, method.locs, method.locs.pred, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, covmodel = cov_matern_spacetime, covparms)
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
  
  ### Generate locs and locs.pred
  if( method.locs %in% c("random", "all.random") ) {
    
    if( method.locs.pred %in% c("random", "separate") ) {
      
      process         <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t, method.locs = "all.random", covmodel = covmodel, covparms = covparms)
      process.pred    <- generate_gp_spacetime(nsim = nsim, n = n.pred, d = d, t.len = t, method.locs = "all.random", covmodel = covmodel, covparms = covparms)

      locs <- list()
      for(i in 1:nsim) locs[[i]] <- process$sim[[i]]$locs

      locs.pred <- list()
      for(i in 1:nsim) locs.pred[[i]] <- process.pred$sim[[i]]$locs
      
      rm(process, process.pred)
      
    } else if( method.locs.pred %in% c("subset", "systematic") ) {
      
      process         <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t, method.locs = "all.random", covmodel = covmodel, covparms = covparms)
      
      locs            <- list()
      locs.pred       <- list()
      for(i in 1:nsim) {
        
        locs[[i]]       <- process$sim[[i]]$locs
        
        ind.unobs       <- sort(sample(seq(n), n.pred))
        
        locs.pred[[i]]  <- locs[[i]][ind.unobs, , drop = FALSE]
        locs[[i]]       <- locs[[i]][-ind.unobs, , drop = FALSE]
      }
      
      rm(process)
      
      # plot(locs[[1]][, 1:2], col = "red", main = "Random locations")
      # points(locs.pred[[1]][, 1:2], col = "blue")
      
    } else {
      
      stop("Please check the argument method.locs.pred!")
    }
    
  } else if( method.locs %in% c("monitoring", "space.random.time.grid") ) {
    
    if( method.locs.pred %in% c("random", "separate") ) {
      
      process         <- generate_gp_spacetime(nsim = nsim, n = n, d = 2, t.len = t, method.locs = "space.random.time.grid", covmodel = covmodel, covparms = covparms)
      process.pred    <- generate_gp_spacetime(nsim = nsim, n = n.pred, d = 2, t.len = 1, method.locs = "all.random", covmodel = covmodel, covparms = covparms)

      locs <- list()
      for(i in 1:nsim) locs[[i]] <- process$sim[[i]]$locs

      locs.pred <- list()
      for(i in 1:nsim) locs.pred[[i]] <- process.pred$sim[[i]]$locs
      
      rm(process, process.pred)
      
    } else if( method.locs.pred %in% c("subset", "systematic") ) {
      
      process         <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t, method.locs = "space.random.time.grid", covmodel = covmodel, covparms = covparms)

      locs            <- list()
      locs.pred       <- list()
      for(i in 1:nsim) {
        
        locs[[i]]       <- as.data.frame(process$sim[[i]]$locs)
        
        locs.spatial    <- dplyr::distinct(locs[[i]][, 1:2])
        locs.spatial    <- cbind(locs.spatial, seq(nrow(locs.spatial)))
        
        locs[[i]]       <- dplyr::left_join(locs[[i]], locs.spatial, by = c("V1", "V2"))
        colnames(locs[[i]]) <- c("coord1", "coord2", "time", "index")
        
        ind.unobs       <- sort(sample(seq(n), n.pred/t))
        bool.unobs      <- locs[[i]]$index %in% ind.unobs
        
        locs.pred[[i]]  <- locs[[i]][bool.unobs, , drop = FALSE]
        locs[[i]]       <- locs[[i]][!bool.unobs, , drop = FALSE]
        
        locs.pred[[i]]  <- as.matrix(locs.pred[[i]][, c("coord1", "coord2", "time")])
        locs[[i]]       <- as.matrix(locs[[i]][, c("coord1", "coord2", "time")])
      }
      
      rm(process)
      
      # plot(locs[[1]][, 1:2], col = "red", main = "Monitoring stations")
      # points(locs.pred[[1]][, 1:2], col = "blue")
      
    } else {
      
      stop("Please check the argument method.locs.pred!")
    }
    
  } else if( method.locs %in% c("satellite") ) {
    
    if( method.locs.pred %in% c("random", "separate") ) {
      
      process         <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t, method.locs = "satellite", covmodel = covmodel, covparms = covparms)
      process.pred    <- generate_gp_spacetime(nsim = nsim, n = n.pred, d = d, t.len = t, method.locs = "all.random", covmodel = covmodel, covparms = covparms)

      locs <- list()
      for(i in 1:nsim) locs[[i]] <- process$sim[[i]]$locs

      locs.pred <- list()
      for(i in 1:nsim) locs.pred[[i]] <- process.pred$sim[[i]]$locs
      
      rm(process, process.pred)
      
    } else if( method.locs.pred %in% c("subset", "systematic") ) {
     
      process         <- generate_gp_spacetime(nsim = nsim, n = n, d = d, t.len = t, method.locs = "satellite", covmodel = covmodel, covparms = covparms)
      
      locs            <- list()
      locs.pred       <- list()
      for(i in 1:nsim) {
        
        locs[[i]]       <- process$sim[[i]]$locs
        
        ind.unobs       <- sample(seq(n), 1)
        ind.unobs       <- seq(from = ind.unobs, by = 1, length.out = n.pred)
        
        ind.unobs[ind.unobs > n] <- ind.unobs[ind.unobs > n] - n
        ind.unobs       <- sort(ind.unobs)
        
        locs.pred[[i]]  <- locs[[i]][ind.unobs, , drop = FALSE]
        locs[[i]]       <- locs[[i]][-ind.unobs, , drop = FALSE]
      }
      
      rm(process)
      
      # plot(locs[[1]][, 1:2], col = "red", main = "Remote sensing")
      # points(locs.pred[[1]][, 1:2], col = "blue")
       
    } else {
      
      stop("Please check the argument method.locs.pred!")
    }
    
  } else {
    
    stop("Please check the argument method.locs!")
  }
  
  ### Generate z and z.pred
  z <- list() ; z.pred <- list()
  for(i in 1:nsim) {
    
    covmat          <- covmodel(locs = rbind(locs[[i]], locs.pred[[i]]), covparms = covparms) # + nuggets * diag(nrow(locs[[i]]) + nrow(locs.pred[[i]]))
    
    z.all           <- as.numeric(t(chol(x = covmat)) %*% as.matrix(rnorm(nrow(locs[[i]]) + nrow(locs.pred[[i]])))) + sqrt(nuggets) * rnorm(nrow(locs[[i]]) + nrow(locs.pred[[i]]))
    z[[i]]          <- z.all[seq(nrow(locs[[i]]))]
    z.pred[[i]]     <- z.all[seq(nrow(locs.pred[[i]])) + nrow(locs[[i]])]
  }
  
  ### body
  result       <- list()
  result[[1]]  <- list(n = n, n.pred = n.pred, d = d, t = t, m = m, covparms = covparms, nuggets = nuggets, locs = locs, locs.pred = locs.pred, z = z, z.pred = z.pred)
  
  simout <- list()
  for(i in 1:5) simout[[i]] <- list(mspe = matrix(NA, nrow = nsim, ncol = length(m)), logscore = matrix(NA, nrow = nsim, ncol = length(m)))
  
  names(simout) <- c("b1", "b2", "b3", "cc", "eu")
  
  for(i in 1:nsim) {
    
    for(j in 1:length(m)) {
      
      out.baseline1   <- prediction_spacetime_baseline(approx = 1, z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covmodel = covmodel, covparms = covparms, nuggets = nuggets, predcond.method = "general", var.exact = TRUE, return.values = "all")
      out.baseline2   <- prediction_spacetime_baseline(approx = 2, z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covmodel = covmodel, covparms = covparms, nuggets = nuggets, predcond.method = "general", var.exact = TRUE, return.values = "all")
      out.baseline3   <- prediction_spacetime_baseline(approx = 3, z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covmodel = covmodel, covparms = covparms, nuggets = nuggets, predcond.method = "general", var.exact = TRUE, return.values = "all")
      out.euclidean   <- prediction_spacetime_cvecchia(method = "euc", z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], initial.pt = NULL, covmodel = covftn, covparms = covparms, nuggets = nuggets, predcond.method = "general", var.exact = TRUE, return.values = "all")
      out.corrvecchia <- prediction_spacetime_cvecchia(method = "cor", z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], initial.pt = NULL, covmodel = covftn, covparms = covparms, nuggets = nuggets, predcond.method = "general", var.exact = TRUE, return.values = "all")
      
      simout$b1$mspe[i, j] <- mean((out.baseline1$predict$mu.pred - z.pred[[i]])^2)
      simout$b2$mspe[i, j] <- mean((out.baseline2$predict$mu.pred - z.pred[[i]])^2)
      simout$b3$mspe[i, j] <- mean((out.baseline3$predict$mu.pred - z.pred[[i]])^2)
      simout$cc$mspe[i, j] <- mean((out.corrvecchia$predict$mu.pred - z.pred[[i]])^2)
      simout$eu$mspe[i, j] <- mean((out.euclidean$predict$mu.pred - z.pred[[i]])^2)
      
      simout$b1$logscore[i, j] <- mean(logscore(out.baseline1$predict$mu.pred, out.baseline1$predict$var.pred, z.pred[[i]]))
      simout$b2$logscore[i, j] <- mean(logscore(out.baseline2$predict$mu.pred, out.baseline2$predict$var.pred, z.pred[[i]]))
      simout$b3$logscore[i, j] <- mean(logscore(out.baseline3$predict$mu.pred, out.baseline3$predict$var.pred, z.pred[[i]]))
      simout$cc$logscore[i, j] <- mean(logscore(out.corrvecchia$predict$mu.pred, out.corrvecchia$predict$var.pred, z.pred[[i]]))
      simout$eu$logscore[i, j] <- mean(logscore(out.euclidean$predict$mu.pred, out.euclidean$predict$var.pred, z.pred[[i]]))
    }
    
    if(verbose == TRUE) { print(i) ; print(Sys.time()) }
  }
  
  result[[2]]  <- simout
  
  mspe.b1 <- colMeans(result[[2]]$b1$mspe)
  mspe.b2 <- colMeans(result[[2]]$b2$mspe)
  mspe.b3 <- colMeans(result[[2]]$b3$mspe)
  mspe.cc <- colMeans(result[[2]]$cc$mspe)
  mspe.eu <- colMeans(result[[2]]$eu$mspe)
  
  logs.b1 <- colMeans(result[[2]]$b1$logscore)
  logs.b2 <- colMeans(result[[2]]$b2$logscore)
  logs.b3 <- colMeans(result[[2]]$b3$logscore)
  logs.cc <- colMeans(result[[2]]$cc$logscore)
  logs.eu <- colMeans(result[[2]]$eu$logscore)
  
  result[[3]] <- list(mspe = list(b1 = mspe.b1, b2 = mspe.b2, b3 = mspe.b3, cc = mspe.cc, eu = mspe.eu), logscore = list(b1 = logs.b1, b2 = logs.b2, b3 = logs.b3, cc = logs.cc, eu = logs.eu))
  
  ### return
  names(result) <- c("setting", "simout", "output")
  
  return(result)
}

