####################################################################
####
####  Author: Myeongjong (MJ) Kang (kmj.stat@gmail.com)
####
####  Description: corrvecchia_knownCovparms(), corrvecchia_unknownCovparms()
####
####################################################################

####################################################################
#### load required packages
####################################################################

require(spatstat)
require(scanstatistics)
require(dplyr)

####################################################################
#### corrvecchia_knownCovparms()
####################################################################

#' @param dat: the observed data
#' @param locs: nxd matrix of observed locs
#' @param m: Number of nearby points to condition on
#' 
#' @param ordering: 'maxmin'
#' @param conditioning: 'NN' (nearest neighbor)
#' 
#' @param covmodel: covariance function
#' @param covparms: covariance parameters as a vector
#' 
corrvecchia_knownCovparms <- function(dat, locs, m, ordering = "maxmin", conditioning = "NN", covmodel, covparms)
{
  
  p     <- ncol(locs)
  n     <- nrow(locs)
  
  ord       <- order_maxmin(locs)
  locsord   <- locs[ord, , drop=FALSE]
  
  dist.matrix   <- distance_correlation(locsord, covmodel, covparms)
  cond.sets     <- conditioning_nn(dist.matrix)
  
  # sparsity structure
  
  vecchia.approx <- generate_vecchia()
  return(vecchia.approx)
  
}

order_maxmin <- function(locs)
{
  
}

distance_correlation <- function(locsord, covmodel, covparms)
{
  
}

conditioning_nn <- function(dist.matrix)
{
  
  n     <- nrow(dist.matrix)
  NN    <- matrix(rep(NA, n^2), nrow = n, ncol = n)
  for(i in 1:n){
    dist.matrix.nested <- dist.matrix[seq(i), seq(i)]
  }
}

### Test validity of dist_to_knn() ###

locs <- matrix(runif(100, 0, 1), 50, 2)
locsord <- locs

nn.spatstat <- nnwhich(locs, k = 1:10)
head(nn.spatstat)

nn.scanstat <- dist_to_knn(dist(locs, diag = T, upper = T), k = 10)
head(nn.scanstat)





