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
require(fields)

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
  cond.sets     <- conditioning_nn(m, dist.matrix)
  
  # sparsity structure
  
  vecchia.approx <- generate_vecchia()
  return(vecchia.approx)
  
}

order_maxmin <- function(locs)
{
  
}

distance_correlation <- function(locsord, covmodel, covparms)
{
  covparms[1]   <- 1 # correlation function
  dist.matrix   <- 1 - covmodel(locsord, covparms) # 1-rho

  return(dist.matrix)
}

# # Test code for distance_correlation()
# locs        <- matrix(runif(20, 0, 1), 10, 2)
# cov.iso     <- function(locs, covparms) covparms[1] * exp(-rdist(locs) / covparms[2])
# cov.aniso   <- function(locs, covparms) covparms[1] * exp(-rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
# covparms    <- c(1, 1, 5)
# 
# distance_correlation(locs, cov.iso, covparms[1:2])
# distance_correlation(locs, cov.aniso, covparms)


conditioning_nn <- function(m, dist.matrix)
{
  # initialize an output object NN which is a n*n matrix
  n     <- nrow(dist.matrix)
  NN    <- matrix(rep(NA, n^2), nrow = n, ncol = n) ; NN[1, 1] <- 1
  
  # Find the nearest neighbor conditioning set for each i-th location using the 'dist_to_knn()' function 
  for(i in 2:n) {
    k               <- min(i, m) # the number of neighbors of the i-th observation
    NN[i, seq(k)]   <- scanstatistics::dist_to_knn(dist.matrix[seq(i), seq(i)], k = k)[i, seq(k)]
  }
  
  return(NN)
}

# # Test code for conditioning_nn()
# locs <- matrix(runif(20, 0, 1), 10, 2)
# dist.matrix <- as.matrix(dist(locs))
# m <- 3
# conditioning_nn(m, dist.matrix)
# 
# # Test validity of dist_to_knn()
# locs <- matrix(runif(100, 0, 1), 50, 2)
# locsord <- locs
#
# nn.spatstat <- nnwhich(locs, k = 1:10)
# head(nn.spatstat)
#
# nn.scanstat <- dist_to_knn(dist(locs, diag = T, upper = T), k = 10)
# head(nn.scanstat)





