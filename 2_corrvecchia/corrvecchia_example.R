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

set.seed(10102019)

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
# 0.05825593 0.05825593 0.05825593 0.05564751

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
# 0.3152489353 0.0003521205 0.3152489353 0.0002679062


####################################################################
#### Nonstaionary Matern covariance model (Risser MD, Calder CA (2015))
####################################################################

# Katzfuss M. 2013. Bayesian nonstationary spatial modeling for very large datasets. Environmetrics 24(3):189–200.
# Stein ML. 2005. Nonstationary spatial covariance functions, Technical Report, University of Chicago, Department of Statistics.

matern_ns <- function(locs1, locs2 = NULL) {
  
  if(is.null(locs2)) locs2 = locs1
  
  n1 <- nrow(locs1) ; n2 <- nrow(locs2) ; d <- ncol(locs1)
  
  mat.cov <- matrix(NA, nrow = n1, ncol = n2)
  
  for(i in 1:n1){
    for(j in 1:n2){
      sigma.ij    <- sigma(locs1[i, ]) * sigma(locs2[j, ])
      kernel.ij   <- ( aniso_mat(locs1[i, ]) + aniso_mat(locs2[j, ]) ) / 2 
      smooth.ij   <- ( smoothness(locs1[i, ]) + smoothness(locs2[j, ]) ) / 2
      
      q.ij <- as.numeric(crossprod( locs1[i, ] - locs2[j, ], solve(kernel.ij, locs1[i, ] - locs2[j, ]) ))
      
      mat.cov[i,j] <- sigma.ij * Matern( sqrt(q.ij), nu = smooth.ij) / sqrt( determinant(kernel.ij, logarithm = F)[[1]][1] )
    }
  }
  
  return(mat.cov)
}


# ##### Test #####
# 
# # spatially-varying standard deviation
# sigma <- function(loc) determinant(aniso_mat(loc), logarithm = F)[[1]][1]^0.25
# # spatially-varying local anisotropy (controlling both the range and direction of dependence)
# aniso_mat<- function(loc) {
#   
#   eta <- 0
#   rot.mat <- matrix(c(cos(eta), -sin(eta), sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
#   
#   range <- c(1, 1)
#   diag.mat <- diag(range, nrow = length(loc))
#   
#   aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
#   
#   return(aniso.mat)
# }
# # matern's smoothness
# smoothness <- function(loc) 0.5
# 
# n             <- 15
# locs          <- matrix(runif(n * 2, 0, 1), n, 2)
# 
# Sigma.new <- matern_ns(locs)
# 
# cov.iso       <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
# covparms      <- c(1, 1, 10)
# 
# Sigma.old <- cov.iso(locs, covparms)
# 
# sum((Sigma.old - Sigma.new)^2)
# 
# ##### Test #####
# 
# # spatially-varying standard deviation
# sigma <- function(loc) determinant(aniso_mat(loc), logarithm = F)[[1]][1]^0.25
# # spatially-varying local anisotropy (controlling both the range and direction of dependence)
# aniso_mat<- function(loc) {
#   
#   eta <- 0
#   rot.mat <- matrix(c(cos(eta), -sin(eta), sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
#   
#   range <- c(0.1^2, 1)
#   diag.mat <- diag(range, nrow = length(loc))
#   
#   aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
#   
#   return(aniso.mat)
# }
# # matern's smoothness
# smoothness <- function(loc) 0.5
# 
# n             <- 15
# locs          <- matrix(runif(n * 2, 0, 1), n, 2)
# 
# Sigma.new <- matern_ns(locs)
# 
# cov.iso       <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
# cov.aniso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])
# covparms      <- c(1, 1, 10)
# 
# Sigma.old <- cov.aniso(locs, covparms)
# 
# sum((Sigma.old - Sigma.new)^2)

####################################################################
#### example 3: basic nonstaionary covariance with function-type covmodel
####################################################################

covparms <- c(1)

# spatially-varying standard deviation
sigma <- function(loc) determinant(aniso_mat(loc), logarithm = F)[[1]][1]^0.25
# spatially-varying local anisotropy (controlling both the range and direction of dependence)
aniso_mat<- function(loc) {
  
  eta <- 0
  rot.mat <- matrix(c(cos(eta), -sin(eta), sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
  
  range <- c(1, 1)
  diag.mat <- diag(range, nrow = length(loc))
  
  aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
  
  return(aniso.mat)
}
# matern's smoothness
smoothness <- function(loc) 0.5

# v <- seq(0, 6, by = 0.05)
# plot(v, Matern(v, smoothness = 0.2), type = 'l', col = 'red', lwd = 3)
# lines(v, Matern(v, smoothness = 0.9), col = 'blue', lwd = 3)
  
n             <- 15^2
m             <- 10
locs          <- matrix(runif(n * 2, 0, 1), n, 2)

# true cov matrix
Sigma <- matern_ns(locs)

matrixcalc::is.positive.definite(aniso_mat(c(1,1)))
temp <- eigen(aniso_mat(c(1,1)))
head(sort(temp$values))

matrixcalc::is.positive.definite(Sigma)
temp <- eigen(Sigma)
head(sort(temp$values))

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
# 0.3970591 0.3984529 0.3921881


####################################################################
#### example 4: basic nonstaionary covariance with matrix-type covmodel
####################################################################

covparms <- c(1)

# spatially-varying standard deviation
sigma <- function(loc) determinant(aniso_mat(loc), logarithm = F)[[1]][1]^0.25
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
smoothness <- function(loc) 0.5

# v <- seq(0, 6, by = 0.05)
# plot(v, Matern(v, smoothness = 0.2), type = 'l', col = 'red', lwd = 3)
# lines(v, Matern(v, smoothness = 0.9), col = 'blue', lwd = 3)

n             <- 15^2
m             <- 10
locs          <- matrix(runif(n * 2, 0, 1), n, 2)

# true cov matrix
Sigma <- matern_ns(locs)

matrixcalc::is.positive.definite(aniso_mat(c(1,1)))
temp <- eigen(aniso_mat(c(1,1)))
head(sort(temp$values))

matrixcalc::is.positive.definite(Sigma)
temp <- eigen(Sigma)
head(sort(temp$values))

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
# 0.3523416 0.3514426 0.3453771


####################################################################
#### example 5: nonstaionary covariance with matrix-type covmodel (non-anisotropy + nonstationarity)
####################################################################

covparms <- c(1)

# spatially-varying standard deviation
sigma <- function(loc) determinant(aniso_mat(loc), logarithm = F)[[1]][1]^0.25
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
smoothness <- function(loc) 0.2 * exp(loc[1])

# v <- seq(0, 6, by = 0.05)
# plot(v, Matern(v, smoothness = 0.2), type = 'l', col = 'red', lwd = 3)
# lines(v, Matern(v, smoothness = 0.9), col = 'blue', lwd = 3)

n             <- 15^2
m             <- 10
locs          <- matrix(runif(n * 2, 0, 1), n, 2)

# true cov matrix
Sigma <- matern_ns(locs)

matrixcalc::is.positive.definite(aniso_mat(c(1,1)))
temp <- eigen(aniso_mat(c(1,1)))
head(sort(temp$values))

matrixcalc::is.positive.definite(Sigma)
temp <- eigen(Sigma)
head(sort(temp$values))

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
# 0.1830565 0.2259558 0.2162229
# 0.1802892 0.1941291 0.1616930


####################################################################
#### example 6: nonstaionary covariance with matrix-type covmodel (anisotropy + nonstaionarity)
####################################################################

covparms <- c(1)

# spatially-varying standard deviation
sigma <- function(loc) determinant(aniso_mat(loc), logarithm = F)[[1]][1]^0.25
# spatially-varying local anisotropy (controlling both the range and direction of dependence)
aniso_mat<- function(loc) {
  
  eta <- pi/4
  rot.mat <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
  
  range <- c(0.1^2, 1)
  diag.mat <- diag(range, nrow = length(loc))
  
  aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
  
  return(aniso.mat)
}
# matern's smoothness
smoothness <- function(loc) 0.2 * exp(loc[1])

# v <- seq(0, 6, by = 0.05)
# plot(v, Matern(v, smoothness = 0.2), type = 'l', col = 'red', lwd = 3)
# lines(v, Matern(v, smoothness = 0.9), col = 'blue', lwd = 3)

n             <- 15^2
m             <- 10
locs          <- matrix(runif(n * 2, 0, 1), n, 2)

# true cov matrix
Sigma <- matern_ns(locs)

matrixcalc::is.positive.definite(aniso_mat(c(1,1)))
temp <- eigen(aniso_mat(c(1,1)))
head(sort(temp$values))

matrixcalc::is.positive.definite(Sigma)
temp <- eigen(Sigma)
head(sort(temp$values))

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
# 5.75016585 0.14757086 0.05779561
# 5.55710955 0.16784603 0.07602673
