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

source("1_Pilot_Study/2_vecchia_specify_adjusted.R")
source("2_corrvecchia/corrvecchia.R")
source("2_corrvecchia/kldiv.R")

set.seed(12345)

####################################################################
#### example 1 : degree of anisotropy = 1
####################################################################

locs          <- matrix(runif(15^2 * 2, 0, 1), 15^2, 2)
n             <- 15^2
m             <- 10

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

locs          <- matrix(runif(15^2 * 2, 0, 1), 15^2, 2)
n             <- 15^2
m             <- 10

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