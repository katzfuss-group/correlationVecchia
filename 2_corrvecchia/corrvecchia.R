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

# require(Rfast)
# 
# require(spatstat)
# require(scanstatistics)
# 
# require(fields)
# require(FNN)
# 
# library(GPvecchia)

####################################################################
#### corrvecchia_knownCovparms()
####################################################################

#' @param locs: nxd matrix of observed locs
#' @param m: Number of nearby points to condition on
#' 
#' @param ordering: 'maxmin'
#' @param ordering.method: 'euclidean', 'correlation'
#' @param conditioning: 'NN' (nearest neighbor)
#' 
#' @param covmodel: covariance function (or matrix)
#' @param covparms: covariance parameters as a vector (variance, range, degree of anisotropy). The first element must be its variance.
#' 
corrvecchia_knownCovparms <- function(locs, m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", covmodel, covparms)
{
  
  p     <- ncol(locs)
  n     <- nrow(locs)
  
  ord           <- if(ordering.method == "euclidean") order_maxmin_euclidean(locs) else order_maxmin_correlation(locs, covmodel, covparms)
  locsord       <- locs[ord, ]
  
  dist.matrix   <- distance_correlation(locsord, covmodel, covparms)
  cond.sets     <- conditioning_nn(m, dist.matrix)
  
  Cond          <- matrix(NA, nrow(cond.sets),ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    
  obs           <- rep(TRUE,n)
  U.prep        <- U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', conditioning = 'NN')
  return(vecchia.approx)
}

# locs          <- matrix(runif(15^2 * 2, 0, 1), 15^2, 2)
# n             <- 15^2
# m             <- 15
# # covmodel
# cov.iso       <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
# cov.aniso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
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
# 
# Sigma.ord       <- cov.aniso(sim.iso$locsord, covparms) # true cov in appropriate ordering
# U               <- createU(sim.iso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
# revord          <- order(sim.iso$ord)
# Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
# source("2_corrvecchia/kldiv.R")
# kls             <- kldiv(Sigma, Sigma.hat)
# kls
# 
# sim.aniso   <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", conditioning = "NN", covmodel = cov.aniso, covparms = covparms)
# 
# Sigma.ord       <- cov.aniso(sim.aniso$locsord, covparms) # true cov in appropriate ordering
# U               <- createU(sim.aniso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
# revord          <- order(sim.aniso$ord)
# Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
# source("2_corrvecchia/kldiv.R")
# kls             <- kldiv(Sigma, Sigma.hat)
# kls

### CAUTION: This function can cause numerical issue. Please use the 'correlation()' function, instead. ###
distance_correlation <- function(locs, covmodel, covparms)
{
  if(is.function(covmodel)) {
    covparms[1]   <- 1 # correlation function
    dist.matrix   <- 1 - covmodel(locs, covparms) # 1-rho
  } else if(is.matrix(covmodel)) {
    dist.matrix   <- 1 - covmodel / covparms[1]
  } else {
    dist.matrix   <- 1 - diag(1, nrow = nrow(locs), ncol = nrow(locs))
  }

  return(dist.matrix)
}

correlation <- function(locs, covmodel, covparms)
{
  if(is.function(covmodel)) {
    covparms[1]   <- 1 # correlation function
    corr.matrix   <- covmodel(locs, covparms) # 1-rho
  } else if(is.matrix(covmodel)) {
    corr.matrix   <- covmodel / covparms[1]
  } else {
    corr.matrix   <- diag(1, nrow = nrow(locs), ncol = nrow(locs))
  }
  
  return(corr.matrix)
}

# # Test code for distance_correlation()
# locs        <- matrix(runif(20, 0, 1), 10, 2)
# cov.iso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
# cov.aniso   <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
# covparms    <- c(1, 1, 5)
# 
# distance_correlation(locs, cov.iso, covparms[1:2])
# distance_correlation(locs, cov.aniso, covparms)
# 
# M <- cov.iso(locs, covparms)
# distance_correlation(locs, M, covparms[1:2])
# 
# M <- cov.aniso(locs, covparms)
# distance_correlation(locs, cov.aniso, covparms)
# 
# # Test code for correlation()
# locs        <- matrix(runif(20, 0, 1), 10, 2)
# cov.iso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
# cov.aniso   <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
# covparms    <- c(1, 1, 5)
# 
# correlation(locs, cov.iso, covparms[1:2])
# correlation(locs, cov.aniso, covparms)
# 
# M <- cov.iso(locs, covparms)
# correlation(locs, M, covparms[1:2])
# 
# M <- cov.aniso(locs, covparms)
# correlation(locs, cov.aniso, covparms)


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


### Caution: old version is literally matched with the algorithm but cause numerical issues. ###
order_maxmin_correlation_old <- function(locs, covmodel, covparms, initial.pt = NULL)
{
  n     <- nrow(locs)
  p     <- ncol(locs)
  ord   <- rep(NA, n)
  d     <- distance_correlation(locs, covmodel, covparms)
  
  if( is.null(initial.pt) ){
    ord[1]        <-  which.min(rowSums(d))
  } else if( is.numeric(initial.pt) ) {
    ord[1]        <- initial.pt
  } else if( initial.pt == 'center' ) {
    cen           <- t(as.matrix(colMeans(locs)))
    ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  } else { # at random
    ord[1]        <- sample(1:n, 1)
  }
  cand.argmax   <- seq(n)[seq(n) != ord[1]]
  
  ind           <- as.matrix(expand.grid(cand.argmax, ord[1]))
  cdist         <- d[ind]
  ord[2]        <- cand.argmax[which.max(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]
  
  for(j in 3:(n-1)){
    ind           <- as.matrix(expand.grid(cand.argmax, ord[seq(j-1)]))
    cdist         <- matrix(d[ind], nrow = length(cand.argmax), ncol = j-1, byrow = F)
    cdist         <- Rfast::rowMins(cdist, value = T)
    ord[j]        <- cand.argmax[which.max(cdist)]
    cand.argmax   <- cand.argmax[cand.argmax != ord[j]]
  } 
  
  ord[n]        <- cand.argmax
  
  return(ord)
}

order_maxmin_correlation <- function(locs, covmodel, covparms, initial.pt = NULL)
{
  n     <- nrow(locs)
  p     <- ncol(locs)
  ord   <- rep(NA, n)
  dinv  <- correlation(locs, covmodel, covparms)
  
  if( is.null(initial.pt) ){
    ord[1]        <-  which.max(rowSums(dinv))
  } else if( is.numeric(initial.pt) ) {
    ord[1]        <- initial.pt
  } else if( initial.pt == 'center' ) {
    cen           <- t(as.matrix(colMeans(locs)))
    ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  } else { # at random
    ord[1]        <- sample(1:n, 1)
  }
  cand.argmax   <- seq(n)[seq(n) != ord[1]]
  
  ind           <- as.matrix(expand.grid(cand.argmax, ord[1]))
  cdist         <- dinv[ind]
  ord[2]        <- cand.argmax[which.min(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]
  
  for(j in 3:(n-1)){
    ind           <- as.matrix(expand.grid(cand.argmax, ord[seq(j-1)]))
    cdist         <- matrix(dinv[ind], nrow = length(cand.argmax), ncol = j-1, byrow = F)
    cdist         <- Rfast::rowMaxs(cdist, value = T)
    ord[j]        <- cand.argmax[which.min(cdist)]
    cand.argmax   <- cand.argmax[cand.argmax != ord[j]]
  } 
  
  ord[n]        <- cand.argmax
  
  return(ord)
}

# # Test code
# locs <- matrix(runif(100, 0, 1), 50, 2)
# cov.iso       <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
# cov.aniso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
# covparms      <- c(1, 0.1, 10)
# covmodel      <- cov.iso
# 
# order_maxmin_euclidean(locs)
# order_maxmin_correlation_old(locs, cov.iso, covparms, initial.pt = 'center')
# order_maxmin_correlation(locs, cov.iso, covparms, initial.pt = 'center')
# GPvecchia::order_maxmin_exact(locs)
# 
# locs.trans <- cbind(locs[,1]*covparms[3],locs[,2]) / covparms[2]
# covmodel      <- cov.aniso
# 
# order_maxmin_euclidean(locs.trans)
# GPvecchia::order_maxmin_exact(locs.trans)
# order_maxmin_correlation_old(locs, cov.aniso, covparms, initial.pt = 2)
# order_maxmin_correlation(locs, cov.aniso, covparms, initial.pt = 2)


conditioning_nn <- function(m, dist.matrix)
{
  # initialize an output object NN which is a n*n matrix
  n     <- nrow(dist.matrix)
  NN    <- matrix(rep(NA, n * (m + 1)), nrow = n, ncol = m + 1) ; NN[1, 1] <- 1
  
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




