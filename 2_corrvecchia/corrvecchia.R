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
#' @param def.dist: either 'NULL' or 'abs'
#' @param ordering.method: either 'euclidean' or 'correlation'
#' @param conditioning: 'NN' (nearest neighbor)
#' @param conditioning.method: either 'euclidean' or 'correlation'
#' @param initial.pt: NULL = which.min(rowMeans(d)), center = euclidean-based center, integer = specify the first obs, 'random' = at random, and else = which.min(rowMeans(d))
#' 
#' @param covmodel: covariance function (or matrix)
#' @param covparms: covariance parameters as a vector (variance, range, degree of anisotropy). The first element must be its variance.
#' 
corrvecchia_knownCovparms <- function(locs, m, ordering = "maxmin", def.dist = NULL, ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel, covparms)
{
  
  p     <- ncol(locs)
  n     <- nrow(locs)
  
  if(ordering.method == "correlation" | conditioning.method == "correlation") {
    rho         <- correlation(locs = locs, covmodel = covmodel, covparms = covparms, def.dist = def.dist)
  }
  
  # ordering
  if(ordering.method == "euclidean") {
    ord         <- order_maxmin_euclidean(locs = locs)
  } else if(ordering.method == "correlation") {
    ord         <- order_maxmin_correlation(locs = locs, dinv = rho, covmodel = covmodel, covparms = covparms, initial.pt = initial.pt)
  } else {
    stop("Please check the ordering method.")
  }

  locsord       <- locs[ord, ]
  
  if(is.matrix(covmodel)) {
    covmodel    <- covmodel[ord, ord]
  }
  
  # conditioning
  if(conditioning.method == "euclidean") {
    cond.sets   <- GpGp::find_ordered_nn(locs = locsord, m = m)
  } else if(conditioning.method == "correlation") {
    rho         <- rho[ord, ord]
    cond.sets   <- conditioning_nn(m = m, dist.matrix = 1 - rho)
  } else {
    stop("Please check the conditioning method.")
  }
  
  Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- GPvecchia::U_sparsity(locsord, cond.sets, obs, Cond)
  
  # output
  vecchia.approx <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = ordering, ordering.method = ordering.method, conditioning = 'NN', conditioning.method = conditioning.method)
  return(vecchia.approx)
}

# # Test code
# # Test code
# library(GPvecchia)
# 
# locs          <- matrix(runif(15^2 * 2, 0, 1), 15^2, 2)
# n             <- 15^2
# m             <- 15
# 
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
# fields::quilt.plot(locs[,1], locs[,2], y)
# 
# out1 <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = cov.aniso, covparms = covparms)
# out2 <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = cov.aniso, covparms = covparms)
# out3 <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = cov.aniso, covparms = covparms)
# out4 <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = cov.aniso, covparms = covparms)
# 
# source("2_corrvecchia/kldiv.R")
# 
# Sigma.ord       <- cov.aniso(out1$locsord, covparms) # true cov in appropriate ordering
# U               <- createU(out1, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
# revord          <- order(out1$ord)
# Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
# kls             <- kldiv(Sigma, Sigma.hat)
# kls
# 
# source("2_corrvecchia/vecchia_specify_adjusted.R")
# outref <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
# 
# Sigma.ord       <- cov.aniso(outref$locsord, covparms) # true cov in appropriate ordering
# U               <- createU(outref, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
# revord          <- order(outref$ord)
# Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
# kls             <- kldiv(Sigma, Sigma.hat)
# kls
# 
# Sigma.ord       <- cov.aniso(out2$locsord, covparms) # true cov in appropriate ordering
# U               <- createU(out2, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
# revord          <- order(out2$ord)
# Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
# kls             <- kldiv(Sigma, Sigma.hat)
# kls
# 
# Sigma.ord       <- cov.aniso(out3$locsord, covparms) # true cov in appropriate ordering
# U               <- createU(out3, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
# revord          <- order(out3$ord)
# Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
# kls             <- kldiv(Sigma, Sigma.hat)
# kls
# 
# Sigma.ord       <- cov.aniso(out4$locsord, covparms) # true cov in appropriate ordering
# U               <- createU(out4, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
# revord          <- order(out4$ord)
# Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
# kls             <- kldiv(Sigma, Sigma.hat)
# kls


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


### CAUTION: This function can cause numerical issue. Please use the 'correlation()' function, instead. ###
distance_correlation <- function(locs, covmodel, covparms, def.dist)
{
  if(is.function(covmodel)) {
    covparms[1]   <- 1 # correlation function
    dist.matrix   <- if(is.null(def.dist)) 1 - covmodel(locs, covparms) else 1 - abs(covmodel(locs, covparms)) # 1-rho
  } else if(is.matrix(covmodel)) {
    dist.matrix   <- if(is.null(def.dist)) 1 - covmodel / covparms[1] else 1 - abs(covmodel) / covparms[1]
  } else {
    dist.matrix   <- 1 - diag(1, nrow = nrow(locs), ncol = nrow(locs))
  }
  
  return(dist.matrix)
}

correlation <- function(locs, covmodel, covparms, def.dist)
{
  if(is.function(covmodel)) {
    covparms[1]   <- 1 # correlation function
    corr.matrix   <- if(is.null(def.dist)) covmodel(locs, covparms) else abs(covmodel(locs, covparms)) # 1-rho
  } else if(is.matrix(covmodel)) {
    corr.matrix   <- if(is.null(def.dist)) covmodel / covparms[1] else abs(covmodel) / covparms[1]
  } else {
    corr.matrix   <- diag(1, nrow = nrow(locs), ncol = nrow(locs))
  }
  
  return(corr.matrix)
}


### Caution: old version is literally matched with the algorithm but cause numerical issues. ###
order_maxmin_correlation_old <- function(locs, d, covmodel, covparms, initial.pt = NULL)
{
  n     <- nrow(locs)
  p     <- ncol(locs)
  ord   <- rep(NA, n)
  
  if( is.null(initial.pt) ){
    ord[1]        <-  which.min(rowSums(d))
  } else if( is.numeric(initial.pt) ) {
    ord[1]        <- initial.pt
  } else if( initial.pt == 'center' ) {
    cen           <- t(as.matrix(colMeans(locs)))
    ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  } else if( initial.pt == 'random') { # at random
    ord[1]        <- sample(1:n, 1)
  } else {
    ord[1]        <-  which.min(rowSums(d))
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

order_maxmin_correlation <- function(locs, dinv, covmodel, covparms, initial.pt = NULL)
{
  n     <- nrow(locs)
  p     <- ncol(locs)
  ord   <- rep(NA, n)
  
  if( is.null(initial.pt) ){
    ord[1]        <-  which.max(rowSums(dinv))
  } else if( is.numeric(initial.pt) ) {
    ord[1]        <- initial.pt
  } else if( initial.pt == 'center' ) {
    cen           <- t(as.matrix(colMeans(locs)))
    ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  } else if( initial.pt == 'random') { # at random
    ord[1]        <- sample(1:n, 1)
  } else {
    ord[1]        <-  which.max(rowSums(dinv))
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
# 
# covmodel      <- cov.iso
# 
# order_maxmin_euclidean(locs)
# GPvecchia::order_maxmin_exact(locs)
# 
# rho <- correlation(locs, cov.iso, covparms, NULL)
# order_maxmin_correlation(locs, rho, cov.iso, covparms, initial.pt = GPvecchia::order_maxmin_exact(locs)[1])
# order_maxmin_correlation_old(locs, 1 - rho, cov.iso, covparms, initial.pt = GPvecchia::order_maxmin_exact(locs)[1])
# 
# locs.trans <- cbind(locs[ ,1] * covparms[3], locs[,2]) / covparms[2]
# covmodel <- cov.aniso
# 
# order_maxmin_euclidean(locs.trans)
# GPvecchia::order_maxmin_exact(locs.trans)
# 
# rho <- correlation(locs, cov.aniso, covparms, NULL)
# order_maxmin_correlation(locs, rho, cov.aniso, covparms, initial.pt = GPvecchia::order_maxmin_exact(locs.trans)[1])
# order_maxmin_correlation_old(locs, 1 - rho, cov.aniso, covparms, initial.pt = GPvecchia::order_maxmin_exact(locs.trans)[1])


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




