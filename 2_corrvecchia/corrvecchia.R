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

require(Rfast)

require(spatstat)
require(scanstatistics)

require(fields)
require(FNN)

####################################################################
#### corrvecchia_knownCovparms()
####################################################################

#' @param locs: nxd matrix of observed locs
#' @param m: Number of nearby points to condition on
#' 
#' @param ordering: 'maxmin'
#' @param conditioning: 'NN' (nearest neighbor)
#' 
#' @param covmodel: covariance function
#' @param covparms: covariance parameters as a vector (variance, range, degree of anisotropy)
#' 
corrvecchia_knownCovparms <- function(locs, m, ordering = "maxmin", conditioning = "NN", covmodel, covparms)
{
  
  p     <- ncol(locs)
  n     <- nrow(locs)
  
  ord           <- order_maxmin_euclidean(locs)
  locsord       <- locs[ord, ]
  
  dist.matrix   <- distance_correlation(locsord, covmodel, covparms)
  cond.sets     <- conditioning_nn(m, dist.matrix)
    
  vecchia.approx <- list(locsord = locsord, ord = ord, ord.pred='general', cond.yz = 'false', conditioning = 'NN')
  return(vecchia.approx)
  
}

# locs          <- matrix(runif(15^2 * 2, 0, 1), 15^2, 2)
# n             <- 15^2
# m             <- 15
# # covmodel
# cov.iso       <- function(locs, covparms) covparms[1] * exp(-rdist(locs) / covparms[2])
# cov.aniso     <- function(locs, covparms) covparms[1] * exp(-rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
# covparms      <- c(1, 0.1, 10)
# 
# # true cov matrix
# Sigma <- cov.aniso(locs, covparms)
# 
# # Visualize the process
# y <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
# quilt.plot(locs[,1], locs[,2], y)
# 
# sim.iso     <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", conditioning = "NN", covmodel = cov.iso, covparms = covparms)
# sim.aniso   <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", conditioning = "NN", covmodel = cov.aniso, covparms = covparms)


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


order_maxmin_euclidean <- function(locs)
{
  n     <- nrow(locs)
  p     <- ncol(locs)
  ord   <- rep(NA, n)
  cen   <- t(as.matrix(colMeans(locs)))
  
  ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  cand.argmax   <- seq(n)[seq(n) != ord[1]]
  
  cdist         <- fields::rdist(locs[cand.argmax, ], matrix(locs[ord[1], ], nrow = 1, ncol = 2))
  ord[2]        <- cand.argmax[which.max(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]
  
  for(j in 3:(n-1)){
    cdist         <- fields::rdist(locs[cand.argmax, ], locs[ord[seq(j-1)], ])
    cdist         <- Rfast::rowMins(cdist, value = T)
    ord[j]        <- cand.argmax[which.max(cdist)]
    cand.argmax   <- cand.argmax[cand.argmax != ord[j]]
  } 
  
  ord[n]        <- cand.argmax
  
  return(ord)
}

# # Test code
# locs <- matrix(runif(100, 0, 1), 50, 2)
# identical(order_maxmin_euclidean(locs), GPvecchia::order_maxmin_exact(locs))
# 
# # Comparison in performance between two different approaches to find minimum
# locs <- matrix(runif(1000, 0, 1), 500, 2)
# library(microbenchmark)
# microbenchmark(which.min(diag(tcrossprod(locs)) -2 * as.numeric(tcrossprod(locs, cen))), which.min(rowSums((locs - matrix(as.numeric(cen), nrow = 500, ncol = 2, byrow = T))^2)))


conditioning_nn <- function(m, dist.matrix)
{
  # initialize an output object NN which is a n*n matrix
  n     <- nrow(dist.matrix)
  NN    <- matrix(rep(NA, n^2), nrow = n, ncol = m + 1) ; NN[1, 1] <- 1
  
  # Find the nearest neighbor conditioning set for each i-th location using the 'dist_to_knn()' function 
  for(i in 2:n) {
    k               <- min(i, m + 1) # the number of neighbors of the i-th observation
    NN[i, seq(k)]   <- scanstatistics::dist_to_knn(dist.matrix[seq(i), seq(i)], k = k)[i, seq(k)]
  }
  
  return(NN)
}

# # Test code for conditioning_nn()
# locs <- matrix(runif(20, 0, 1), 10, 2)
# dist.matrix <- as.matrix(dist(locs))
# m <- 3
# conditioning_nn(m, dist.matrix)
# GpGp::find_ordered_nn(locs, m)
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
# 
# # Comparison in performance between get.knn and dist_to_knn
# library(microbenchmark)
# microbenchmark(FNN::get.knn(locs, k = m), scanstatistics::dist_to_knn(dist(locs, diag = T, upper = T), k = 3))




