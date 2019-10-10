####################################################################
####
####  Author: Myeongjong (MJ) Kang (kmj.stat@gmail.com)
####
####  Description:
####
####################################################################

gc()
rm(list = ls())

library(GPvecchia)
library(fields)

source("1_Pilot_Study/2_vecchia_specify_adjusted.R")
source("2_corrvecchia/corrvecchia.R")
source("2_corrvecchia/kldiv.R")

set.seed(12345)

####################################################################
#### example 1 : degree of anisotropy = 1
####################################################################

n             <- 15^2
m             <- 10
locs          <- matrix(runif(n * 2, 0, 1), n, 2)

# covmodel
cov.iso       <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
cov.aniso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
covparms      <- c(1, 0.1, 1)

locs.trans    <- cbind(locs[, 1] * covparms[3], locs[, 2])

# true cov matrix
Sigma <- cov.aniso(locs, covparms)

# Visualize the process
y <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
fields::quilt.plot(locs[,1], locs[,2], y)


### specify vecchia approximations
approx <- list()

# standard vecchia
approx[[1]]           <- vecchia_specify_adjusted(locs, m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
# pilot correlation-based vecchia
approx[[2]]           <- vecchia_specify_adjusted(locs.trans, m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
approx[[2]]$locsord   <- locs[approx[[2]]$ord,]
# standard vecchia with the corrvecchia function
approx[[3]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", covmodel = cov.iso, covparms = covparms)
# correlation-based vecchia with the corrvecchia function
approx[[4]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", covmodel = cov.aniso, covparms = covparms)


### compute approximate covariance matrices
Sigma.hat   <- list()
kls         <- c()
for(i in 1:4){
  
  Sigma.ord       <- cov.aniso(approx[[i]]$locsord, covparms) # true cov in appropriate ordering
  
  U               <- createU(approx[[i]], c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
  revord          <- order(approx[[i]]$ord)
  Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
  
  kls[i]          <- kldiv(Sigma, Sigma.hat[[i]])
}

kls
# standard vecchia / pilot correlation-based vecchia / standard vecchia with the corrvecchia function / correlation-based vecchia with the corrvecchia function
# 0.05371653 0.05371653 0.05371653 0.05308003

####################################################################
#### example 2 : degree of anisotropy = 10
####################################################################

n             <- 15^2
m             <- 10
locs          <- matrix(runif(n * 2, 0, 1), n, 2)

# covmodel
cov.iso       <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
cov.aniso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
covparms      <- c(1, 0.1, 10)

locs.trans    <- cbind(locs[, 1] * covparms[3], locs[, 2])

# true cov matrix
Sigma <- cov.aniso(locs, covparms)

# Visualize the process
y <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
fields::quilt.plot(locs[,1], locs[,2], y)


### specify vecchia approximations
approx <- list()

# standard vecchia
approx[[1]]           <- vecchia_specify_adjusted(locs, m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
# pilot correlation-based vecchia
approx[[2]]           <- vecchia_specify_adjusted(locs.trans, m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
approx[[2]]$locsord   <- locs[approx[[2]]$ord,]
# standard vecchia with the corrvecchia function
approx[[3]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", covmodel = cov.iso, covparms = covparms)
# correlation-based vecchia with the corrvecchia function
approx[[4]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", covmodel = cov.aniso, covparms = covparms)


### compute approximate covariance matrices
Sigma.hat   <- list()
kls         <- c()
for(i in 1:4){
  
  Sigma.ord       <- cov.aniso(approx[[i]]$locsord, covparms) # true cov in appropriate ordering
  
  U               <- createU(approx[[i]], c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
  revord          <- order(approx[[i]]$ord)
  Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
  
  kls[i]          <- kldiv(Sigma, Sigma.hat[[i]])
}

kls
# standard vecchia / pilot correlation-based vecchia / standard vecchia with the corrvecchia function / correlation-based vecchia with the corrvecchia function
# 0.4090684634 0.0002972307 0.4090684634 0.0003749438

####################################################################
#### example 3: basic nonstaionary covariance with function-type covmodel
####################################################################

covparms <- c(1)

# spatially-varying standard deviation
sigma <- function(loc) 1
# spatially-varying local anisotropy (controlling both the range and direction of dependence)
aniso_mat<- function(loc) {
  
  eta <- 0
  rot.mat <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
  
  range <- c(1, 1)
  diag.mat <- diag(range, nrow = length(loc))
  
  aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
  
  return(aniso.mat)
}
# matern's smoothness
smoothness <- function(loc) 1

# v <- seq(0, 6, by = 0.05)
# plot(v, Matern(v, smoothness = 0.2), type = 'l', col = 'red', lwd = 3)
# lines(v, Matern(v, smoothness = 0.9), col = 'blue', lwd = 3)
  
# Risser MD, Calder CA (2015)
matern_ns <- function(locs1, locs2 = NULL) {
  
  if(is.null(locs2)) locs2 = locs1
  
  n1 <- nrow(locs1) ; n2 <- nrow(locs2) ; d <- ncol(locs1)
  
  mat.cov <- matrix(NA, nrow = n1, ncol = n2)
  
  for(i in 1:n1){
    for(j in 1:n2){
      
      sigma.i <- sigma(locs1[i, ]) ; sigma.j <- sigma(locs2[j, ])
      SIGMA.i <- aniso_mat(locs1[i, ]) ; SIGMA.j <- aniso_mat(locs2[j, ])
      sigma.ij <- sigma.i * sigma.j
      norm.ij <- (SIGMA.i + SIGMA.j) / 2 
      smooth.ij <- ( smoothness(locs1[i ])+smoothness(locs2[j, ]) ) / 2
      
      q.ij <- crossprod( locs1[i, ] - locs2[j, ], solve(norm.ij, locs1[i, ] - locs2[j, ]) )
  
      mat.cov[i,j] <- sigma.ij * Matern( sqrt(q.ij), nu = smooth.ij) / sqrt( determinant(norm.ij, logarithm = F)[[1]][1] )
    }
  }
  
  return(mat.cov)
}

n             <- 15^2
m             <- 10
locs          <- matrix(runif(n * 2, 0, 1), n, 2)

# true cov matrix
Sigma <- matern_ns(locs)

# Visualize the process
y <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
fields::quilt.plot(locs[,1], locs[,2], y)


### specify vecchia approximations
approx <- list()
covmodel <- function(locs, covparms) matern_ns(locs) / covparms[1]

# standard vecchia
approx[[1]]           <- vecchia_specify_adjusted(locs, m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
# standard vecchia with the corrvecchia function
approx[[2]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", covmodel = covmodel, covparms = covparms)
# correlation-based vecchia with the corrvecchia function
approx[[3]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", covmodel = covmodel, covparms = covparms)

### compute approximate covariance matrices
Sigma.hat   <- list()
kls         <- c()
for(i in 1:3){
  
  Sigma.ord       <- matern_ns(approx[[i]]$locsord) # true cov in appropriate ordering
  
  U               <- createU(approx[[i]], c(1, 1, 1), 0, covmodel = Sigma.ord)$U
  revord          <- order(approx[[i]]$ord)
  Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
  
  kls[i]          <- kldiv(Sigma, Sigma.hat[[i]])
}

kls
# standard vecchia / standard vecchia with the corrvecchia function / correlation-based vecchia with the corrvecchia function
# 2.183058 2.178872 2.178872

####################################################################
#### example 4: basic nonstaionary covariance with matrix-type covmodel
####################################################################

covparms <- c(1)

# spatially-varying standard deviation
sigma <- function(loc) 1
# spatially-varying local anisotropy (controlling both the range and direction of dependence)
aniso_mat<- function(loc) {
  
  eta <- 0
  rot.mat <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
  
  range <- c(1, 1)
  diag.mat <- diag(range, nrow = length(loc))
  
  aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
  
  return(aniso.mat)
}
# matern's smoothness
smoothness <- function(loc) 1

# v <- seq(0, 6, by = 0.05)
# plot(v, Matern(v, smoothness = 0.2), type = 'l', col = 'red', lwd = 3)
# lines(v, Matern(v, smoothness = 0.9), col = 'blue', lwd = 3)

# Risser MD, Calder CA (2015)
matern_ns <- function(locs1, locs2 = NULL) {
  
  if(is.null(locs2)) locs2 = locs1
  
  n1 <- nrow(locs1) ; n2 <- nrow(locs2) ; d <- ncol(locs1)
  
  mat.cov <- matrix(NA, nrow = n1, ncol = n2)
  
  for(i in 1:n1){
    for(j in 1:n2){
      
      sigma.i <- sigma(locs1[i, ]) ; sigma.j <- sigma(locs2[j, ])
      SIGMA.i <- aniso_mat(locs1[i, ]) ; SIGMA.j <- aniso_mat(locs2[j, ])
      sigma.ij <- sigma.i * sigma.j
      norm.ij <- (SIGMA.i + SIGMA.j) / 2 
      smooth.ij <- ( smoothness(locs1[i ])+smoothness(locs2[j, ]) ) / 2
      
      q.ij <- crossprod( locs1[i, ] - locs2[j, ], solve(norm.ij, locs1[i, ] - locs2[j, ]) )
      
      mat.cov[i,j] <- sigma.ij * Matern( sqrt(q.ij), nu = smooth.ij) / sqrt( determinant(norm.ij, logarithm = F)[[1]][1] )
    }
  }
  
  return(mat.cov)
}

n             <- 15^2
m             <- 10
locs          <- matrix(runif(n * 2, 0, 1), n, 2)

# true cov matrix
Sigma <- matern_ns(locs)

# Visualize the process
y <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
fields::quilt.plot(locs[,1], locs[,2], y)


### specify vecchia approximations
approx <- list()

# standard vecchia
approx[[1]]           <- vecchia_specify_adjusted(locs, m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
# standard vecchia with the corrvecchia function
approx[[2]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", covmodel = Sigma, covparms = covparms)
# correlation-based vecchia with the corrvecchia function
approx[[3]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", covmodel = Sigma, covparms = covparms)

### compute approximate covariance matrices
Sigma.hat   <- list()
kls         <- c()
for(i in 1:3){
  
  Sigma.ord       <- matern_ns(approx[[i]]$locsord) # true cov in appropriate ordering
  
  U               <- createU(approx[[i]], c(1, 1, 1), 0, covmodel = Sigma.ord)$U
  revord          <- order(approx[[i]]$ord)
  Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
  
  kls[i]          <- kldiv(Sigma, Sigma.hat[[i]])
}

kls
# standard vecchia / standard vecchia with the corrvecchia function / correlation-based vecchia with the corrvecchia function
# 2.398728 2.398728 2.398728