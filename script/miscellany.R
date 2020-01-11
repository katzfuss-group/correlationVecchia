####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to evaluate performance of vecchia approximation approaches.
###
####################################################################################



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
