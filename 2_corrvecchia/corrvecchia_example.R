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

source("2_corrvecchia/corrvecchia.R")
source("2_corrvecchia/kldiv.R")

set.seed(12345)

####################################################################
#### example 1 : degree of anisotropy = 1
####################################################################

locs          <- matrix(runif(15^2 * 2, 0, 1), 15^2, 2)
n             <- 15^2
m             <- 15
# covmodel
cov.iso       <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
cov.aniso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
covparms      <- c(1, 0.1, 1)

# true cov matrix
Sigma <- cov.aniso(locs, covparms)

# Visualize the process
y <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
fields::quilt.plot(locs[,1], locs[,2], y)

sim.iso     <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", covmodel = cov.iso, covparms = covparms)

Sigma.ord       <- cov.aniso(sim.iso$locsord, covparms) # true cov in appropriate ordering
U               <- createU(sim.iso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
revord          <- order(sim.iso$ord)
Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
kls             <- kldiv(Sigma, Sigma.hat)
kls

sim.iso     <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", conditioning = "NN", covmodel = cov.iso, covparms = covparms)

Sigma.ord       <- cov.aniso(sim.iso$locsord, covparms) # true cov in appropriate ordering
U               <- createU(sim.iso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
revord          <- order(sim.iso$ord)
Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
kls             <- kldiv(Sigma, Sigma.hat)
kls

sim.aniso   <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", covmodel = cov.aniso, covparms = covparms)

Sigma.ord       <- cov.aniso(sim.aniso$locsord, covparms) # true cov in appropriate ordering
U               <- createU(sim.aniso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
revord          <- order(sim.aniso$ord)
Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
kls             <- kldiv(Sigma, Sigma.hat)
kls

sim.aniso   <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", conditioning = "NN", covmodel = cov.aniso, covparms = covparms)

Sigma.ord       <- cov.aniso(sim.aniso$locsord, covparms) # true cov in appropriate ordering
U               <- createU(sim.aniso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
revord          <- order(sim.aniso$ord)
Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
kls             <- kldiv(Sigma, Sigma.hat)
kls

####################################################################
#### example 1 : degree of anisotropy = 10
####################################################################

locs          <- matrix(runif(15^2 * 2, 0, 1), 15^2, 2)
n             <- 15^2
m             <- 15
# covmodel
cov.iso       <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
cov.aniso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
covparms      <- c(1, 0.1, 10)

# true cov matrix
Sigma <- cov.aniso(locs, covparms)

# Visualize the process
y <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
fields::quilt.plot(locs[,1], locs[,2], y)

sim.iso     <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", covmodel = cov.iso, covparms = covparms)

Sigma.ord       <- cov.aniso(sim.iso$locsord, covparms) # true cov in appropriate ordering
U               <- createU(sim.iso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
revord          <- order(sim.iso$ord)
Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
kls             <- kldiv(Sigma, Sigma.hat)
kls

sim.iso     <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", conditioning = "NN", covmodel = cov.iso, covparms = covparms)

Sigma.ord       <- cov.aniso(sim.iso$locsord, covparms) # true cov in appropriate ordering
U               <- createU(sim.iso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
revord          <- order(sim.iso$ord)
Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
kls             <- kldiv(Sigma, Sigma.hat)
kls

sim.aniso   <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", covmodel = cov.aniso, covparms = covparms)

Sigma.ord       <- cov.aniso(sim.aniso$locsord, covparms) # true cov in appropriate ordering
U               <- createU(sim.aniso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
revord          <- order(sim.aniso$ord)
Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
kls             <- kldiv(Sigma, Sigma.hat)
kls

sim.aniso   <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", conditioning = "NN", covmodel = cov.aniso, covparms = covparms)

Sigma.ord       <- cov.aniso(sim.aniso$locsord, covparms) # true cov in appropriate ordering
U               <- createU(sim.aniso, c(1, 0.1, 0.5), 0, covmodel = Sigma.ord)$U
revord          <- order(sim.aniso$ord)
Sigma.hat       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
kls             <- kldiv(Sigma, Sigma.hat)
kls
