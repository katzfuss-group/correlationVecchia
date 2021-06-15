####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several small functions to evaluate performance of vecchia approximations.
###
###   Contents: 
###     kldiv / performance
###     loop_kldiv / loop_meankldiv / loop_loglik
###
####################################################################################

#' @title Multivariate normal KL divergence between true distribution N(mu0, covmat0) and approximate distribution N(mu1, covmat1)
#'
#' @param covmat0 A covariance matrix of the true distribution
#' @param covmat1 A covariance matrix of the approximate distribution
#' @param mu0 A mean vector of the true distribution assigned to zero vector by default
#' @param mu1 A mean vector of the approximate distribution assigned to zero vector by default
#' @param chol1 if cholesky factor of covmat1 has already been computed, provide chol1 instead of covmat1
#' @param method.modify An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                      If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                      If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                      If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                      Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'                      Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix. 
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix. At \code{FALSE} by default
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#'
#' @return KL divergence between true and approx distributions
#' 
#' @export
#'
#' @examples
#' kldiv(diag(1, 4), diag(2, 4))
kldiv <- function(covmat0, covmat1, mu0 = rep(0, nrow(covmat0)), mu1 = rep(0, nrow(covmat0)), chol1 = NULL, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps)
{
  n           <- nrow(covmat0)
  
  chol0       <- tryCatch(t( chol(covmat0) ), error = function(e) "Oops")
  if(is.character(chol0)) {
    
    chol0       <- factorize(covmat0, pivot = pivot, method = method.modify, tol = tol, return.err = FALSE, verbose = FALSE)
    chol0       <- t(chol0$covfactor)
  }

  if(is.null(chol1)) {
  
    chol1       <- tryCatch(t( chol(covmat1) ), error = function(e) "Oops")
    if(is.character(chol1)) {
      
      chol1       <- factorize(covmat1, pivot = pivot, method = method.modify, tol = tol, return.err = FALSE, verbose = FALSE)
      chol1       <- t(chol1$covfactor)
    }
  }
  
  # trace term
  M           <- solve( chol1, chol0 )
  traceterm   <- sum( M^2 ) - n
  
  # log det term
  logdetterm  <- 2 * sum( log(diag(chol1)) ) - 2 * sum( log(diag(chol0)) )
  
  # mean term (equal to zero if means are zero)
  meandiff    <- mu1 - mu0
  
  # t(meandiff) %*% solve(covmat1, meandiff)
  temp        <- forwardsolve( chol1, meandiff )
  meanterm    <- sum( temp^2 )
  
  # kldiv
  kld         <- 1/2 * ( traceterm + logdetterm + meanterm )
  return(kld)
}

#' @title Evaluating performance of Vecchia approximations in terms of KL divergence
#'
#' @param vecchia.approx object returned by specifying functions
#' @param locs A matrix of locations
#' @param covmodel A covariance model
#' @param covparms A numeric vector of covariance parameters
#' @param nuggets Nuggets at 0 by default
#'
#' @return KL divergence between true and approximate covariance matrices
#' 
#' @export
#'
#' @examples
#' kldiv(diag(1, 4), diag(2, 4))
performance <- function(vecchia.approx, locs, covmodel, covparms, nuggets = 0)
{
  Sigma           <- covmodel(locs, covparms)
  
  Sigma.ord       <- covmodel(vecchia.approx$locsord, covparms) # true cov in appropriate ordering
  U               <- GPvecchia::createU(vecchia.approx, covparms, nuggets, covmodel = Sigma.ord)$U # vecchia.approx, covparms, nuggets, covmodel
  revord          <- order(vecchia.approx$ord)
  Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
  
  kld             <- kldiv(Sigma, Sigma.hat)
  
  return(kld)
}

#################################################################################### For Fisher-scoring approach

.loop_kldiv <- function(nsim, ms, approx, locs.ls, covmat.ls, fit.ls, verbose = TRUE)
{
  covmodel    <- function(locs, covparms) GpGp::matern_spacetime(covparms = covparms, locs = locs)
  
  result      <- list()
  for(j in 1:length(ms)) {
    
    if(verbose == TRUE) print(j)
    
    result[[j]] <- list()
    for(i in 1:nsim) {
      
      if( approx %in% c("b1", "T-ord + T-NN") ) {
        
        out_specify <- baseline_1_spacetime_specify(locs = locs.ls[[i]], m = ms[j])
        
      } else if( approx %in% c("b2", "T-ord + E-NN") ) {
        
        out_specify <- baseline_2_spacetime_specify(locs = locs.ls[[i]], m = ms[j])
        
      } else if( approx %in% c("b3", "T-ord + C-NN") ) {
        
        out_specify <- baseline_3_spacetime_specify(locs = locs.ls[[i]], m = ms[j], covmodel = covmodel, covparms = fit.ls[[j]][[i]]$covparms)
        
      } else if( approx %in% c("cc", "C-MM + C-NN") ) {
        
        out_specify <- cvecchia_m_specify(locs = locs.ls[[i]], m = ms[j], ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor", covmodel = covmodel, covparms = fit.ls[[j]][[i]]$covparms)
        
      } else {
        
        stop("The argument approx is invalid!")
        
      }
      
      indmat <- matrix(0, nrow = nrow(locs.ls[[i]]), ncol = nrow(locs.ls[[i]]))
      for(k in 1:nrow(locs.ls[[i]])) indmat[k, out_specify$U.prep$revNNarray[k, ]] <- 1
      indmat <- indmat + t(indmat) ; diag(indmat) <- 1
      # fields::image.plot(indmat)
      
      # ?GpGp::matern_spacetime
      covmat_est  <- covmodel(locs = locs.ls[[i]], covparms = fit.ls[[j]][[i]]$covparms)
      mu_est      <- rep(fit.ls[[j]][[i]]$betahat, nrow(locs.ls[[i]]))
      
      revord      <- order(out_specify$ord)
      covmat_ord  <- covmat_est[out_specify$ord, out_specify$ord]
      
      U           <- GPvecchia::createU(out_specify, c(1), 0, covmodel = covmat_ord)$U
      covmat_est  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
      
      # ?correlationVecchia::kldiv
      result[[j]][[i]] <- kldiv(covmat0 = covmat.ls[[i]], covmat1 = covmat_est, mu0 = rep(0, length(mu_est)), mu1 = mu_est)
    }
  }
  
  return(result)
}

.loop_meankldiv <- function(nsim, m, kldiv.ls)
{
  result <- rep(0, length(m))
  
  for(j in 1:length(m)) {
    for(i in 1:nsim) {
      result[j] <- result[j] + kldiv.ls[[j]][[i]]
    }
  }
  
  return(result / length(m))
}

.loop_loglik <- function(y, locs, fit, approx, covmodel)
{
  covmat_est  <- covmodel(locs = locs, covparms = fit$covparms)
  mu_est      <- rep(fit$betahat, nrow(locs))
  
  revord      <- order(approx$ord)
  covmat_ord  <- covmat_est[approx$ord, approx$ord]
  
  U           <- GPvecchia::createU(approx, c(1), 0, covmodel = covmat_ord)$U
  covmat_est  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
  
  return(mvtnorm::dmvnorm(as.numeric(y), mu_est, covmat_est, log = TRUE))
}

#################################################################################### For prediction

#' Simple function for calculating logarithmic prediction score
#'
#' @param mu.pred predicted mean vector
#' @param var.pred predicted variance vector
#' @param z observations of the process of interest
#'
#' @return A numeric value
#' 
#' @export
#'
#' @examples
#' 1 + 1
logscore <- function(mu.pred, var.pred, z)
{
  S <- dnorm(x = z, mean = mu.pred, sd = sqrt(var.pred), log = TRUE)
  
  return(S)
}