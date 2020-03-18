####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several small functions to evaluate performance of vecchia approximations.
###
###   Contents:
###       FTN - kldiv()
###       FTN - performance()
###
####################################################################################



#' @title Multivariate normal KL divergence between true distribution N(mu0, covmat0) and approximate distribution N(mu1, covmat1)
#'
#' @param covmat0 A covariance matrix of the true distribution
#' @param covmat1 A covariance matrix of the approximate distribution
#' @param mu0 A mean vector of the true distribution assigned to zero vector by default
#' @param mu1 A mean vector of the approximate distribution assigned to zero vector by default
#' @param chol1 if cholesky factor of covmat1 has already been computed, provide chol1 instead of covmat1
#'
#' @return KL divergence between true and approx distributions
#' 
#' @export
#'
#' @examples
#' kldiv(diag(1, 4), diag(2, 4))
kldiv <- function(covmat0, covmat1, mu0 = rep(0, nrow(covmat0)), mu1 = rep(0, nrow(covmat0)), chol1 = NULL)
{
  n       <- nrow(covmat0)
  chol0   <- t(chol(covmat0))
  if(is.null(chol1)) chol1 <- t(chol(covmat1))
  
  # trace term
  M           <- solve( chol1, chol0 )
  traceterm   <- sum( M^2 ) - n
  
  # log det term
  logdetterm  <- 2 * sum(log(diag(chol1))) - 2 * sum(log(diag(chol0)))
  
  # mean term (equal to zero if means are zero)
  meandiff    <- mu1 - mu0
  
  # t(meandiff) %*% solve(covmat1,meandiff)
  temp        <- forwardsolve(chol1, meandiff)
  meanterm    <- sum(temp^2)
  
  # kldiv
  kldiv <- 1/2 * ( traceterm + logdetterm + meanterm )
  return(kldiv)
}



#' Evaluating performance of Vecchia approximations in terms of KL divergence
#'
#' @param out vecchia_specify
#' @param locs a matrix of locations
#' @param covmodel covariance function
#' @param covparms a numeric vector of covariance parameters
#'
#' @return KL divergence between true and approx distributions
#' @export
#'
performance <- function(out, locs, covmodel, covparms)
{
  Sigma           <- covmodel(locs, covparms)
  Sigma.ord       <- covmodel(out$locsord, covparms) # true cov in appropriate ordering
  U               <- GPvecchia::createU(out, covparms, 0, covmodel = Sigma.ord)$U
  revord          <- order(out$ord)
  Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
  kls             <- kldiv(Sigma, Sigma.hat)
  
  return(kls)
}

