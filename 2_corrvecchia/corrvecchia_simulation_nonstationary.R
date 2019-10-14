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

library(foreach)

source("2_corrvecchia/vecchia_specify_adjusted.R")
source("2_corrvecchia/corrvecchia.R")
source("2_corrvecchia/kldiv.R")

set.seed(10102019)

####################################################################
#### Nonstaionary Matern covariance model (Risser MD, Calder CA (2015))
####################################################################

covparms <- c(1)

a <- function(loc) 0.47 * loc[1] + 0.03
b <- function(loc) 1
angle <- function(loc) 0

# spatially-varying standard deviation
sigma <- function(loc) determinant(aniso_mat(loc), logarithm = F)[[1]][1]^0.25
# spatially-varying local anisotropy (controlling both the range and direction of dependence)
aniso_mat<- function(loc) {
  
  eta <- angle(loc)
  rot.mat <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
  
  range <- c(a(loc)^(-2), b(loc)^(-2))
  diag.mat <- diag(range, nrow = length(loc))
  
  aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
  
  return(aniso.mat)
}
# matern's smoothness
smoothness <- function(loc) 0.2 * exp(loc[1])

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
      
      mat.cov[i,j] <- sigma.ij * fields::Matern( sqrt(q.ij), nu = smooth.ij) / sqrt( determinant(kernel.ij, logarithm = F)[[1]][1] )
    }
  }
  
  return(mat.cov)
}

####################################################################
#### simulation function
####################################################################

simulation <- function(n = 15^2, m = 10, covparms = c(1)) {

  locs      <- matrix(runif(n * 2, 0, 1), n, 2)
  Sigma     <- matern_ns(locs)
  y         <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
  
  ### specify vecchia approximations
  approx <- list()
  
  # standard vecchia with maxmin ordering
  approx[[1]]           <- vecchia_specify_adjusted(locs, m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
  # standard vecchia with x coord ordering
  approx[[2]]           <- vecchia_specify_adjusted(locs, m, ordering = "coord", which.coord = 1, cond.yz='y', conditioning = "NN")
  # standard vecchia with y coord ordering
  approx[[3]]           <- vecchia_specify_adjusted(locs, m, ordering = "coord", which.coord = 2, cond.yz='y', conditioning = "NN")
  # correlation-based vecchia with the corrvecchia function
  approx[[4]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", covmodel = Sigma, covparms = covparms)
  
  ### compute approximate covariance matrices
  Sigma.hat   <- list()
  kls         <- c()
  for(i in 1:4){
    
    Sigma.ord       <- matern_ns(approx[[i]]$locsord) # true cov in appropriate ordering
    
    U               <- createU(approx[[i]], c(1, 1, 1), 0, covmodel = Sigma.ord)$U
    revord          <- order(approx[[i]]$ord)
    Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    kls[i]          <- kldiv(Sigma, Sigma.hat[[i]])
  }
  
  result                  <- list()
  result$n                <- n
  result$m                <- m
  result$covparms         <- covparms
  result$locs             <- locs
  result$approx           <- approx
  result$kls              <- kls
  result$Sigma            <- Sigma
  result$Sigma.hat        <- Sigma.hat
  
  return(result)  
}


####################################################################
#### simulation 1: smoothness
####################################################################

cand.m    <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45) ; n.cand.m <- length(cand.m)
sim1      <- list()

# # small case
# cand.m    <- c(1, 10, 20, 30, 40) ; n.cand.m <- length(cand.m)
# sim1      <- list()

a <- function(loc) 1
b <- function(loc) 1
angle <- function(loc) 0
smoothness <- function(loc) 0.2 + 1.3 * loc[1]

# plot(seq(0, 5, 0.01), fields::Matern(seq(0, 5, 0.01), smoothness = 1), type = "l", col = 1)
# lines(seq(0, 5, 0.01), fields::Matern(seq(0, 5, 0.01), smoothness = 0.2), type = "l", col = 2)
# lines(seq(0, 5, 0.01), fields::Matern(seq(0, 5, 0.01), smoothness = 1.5), type = "l", col = 3)
# lines(seq(0, 5, 0.01), fields::Matern(seq(0, 5, 0.01), smoothness = 2), type = "l", col = 4)

no_cores            <- parallel::detectCores() - 2
cl                  <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)
sim1 <- foreach(m = cand.m, .export = c("a", "b", "aniso_mat", "conditioning_nn", "correlation", "corrvecchia_knownCovparms", "distance_correlation", "kldiv", "matern_ns", "order_maxmin_correlation", "order_maxmin_correlation_old", "order_maxmin_euclidean", "simulation", "smoothness", "vecchia_specify_adjusted"), .packages='GPvecchia') %dopar% simulation(30^2, m = m, covparms = c(1))
parallel::stopCluster(cl)

kls.maxmin.euclidean    <- rep(NA, n.cand.m)
kls.maxmin.corr         <- rep(NA, n.cand.m)
kls.xcoord.euclidean    <- rep(NA, n.cand.m)
kls.ycoord.euclidean    <- rep(NA, n.cand.m)
for(i in 1:n.cand.m) {
  kls.maxmin.euclidean[i]    <- sim1[[i]]$kls[1]
  kls.maxmin.corr[i]         <- sim1[[i]]$kls[4]
  kls.xcoord.euclidean[i]    <- sim1[[i]]$kls[2]
  kls.ycoord.euclidean[i]    <- sim1[[i]]$kls[3]
}

vis.dat1 <- data.frame(kls.maxmin.euclidean, kls.maxmin.corr, kls.xcoord.euclidean, kls.ycoord.euclidean)
vis.dat1 <- vis.dat1[, order(colnames(vis.dat1))]
head(vis.dat1)

plot(cand.m, log10(vis.dat1$kls.maxmin.euclidean), type = "o", col = 1, lty = 1, lwd = 3,
     ylim = c(min(log10(vis.dat1)), max(log10(vis.dat1))), xlab = "m", ylab = "log10(KL)", main = NULL)
lines(cand.m, log10(vis.dat1$kls.maxmin.corr), type = "o", col = 2, lty = 2, lwd = 3)
lines(cand.m, log10(vis.dat1$kls.xcoord.euclidean), type = "o", col = 3, lty = 3, lwd = 3)
lines(cand.m, log10(vis.dat1$kls.ycoord.euclidean), type = "o", col = 4, lty = 4, lwd = 3)
legend("topright", legend=c("maxmin + Euclidean", "maxmin + correlation", "x-coord + Euclidean", "y-coord + Euclidean"), col=1:4, lty=1:4, lwd = 3, cex=1)

# save(sim1, cand.m, vis.dat1, file='2_corrvecchia/sim_nonstationarity_1.RData')
# load(file='2_corrvecchia/sim_nonstationarity_1.RData')


####################################################################
#### simulation 2: range
####################################################################

cand.m    <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45) ; n.cand.m <- length(cand.m)
sim2      <- list()

# # small case
# cand.m    <- c(1, 10, 20, 30, 40) ; n.cand.m <- length(cand.m)
# sim2      <- list()

a <- function(loc) 1/(0.03 + 0.97 * loc[1])
b <- function(loc) 1/(0.03 + 0.97 * loc[1])
angle <- function(loc) 0
smoothness <- function(loc) 0.5

no_cores            <- parallel::detectCores() - 2
cl                  <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)
sim2 <- foreach(m = cand.m, .export = c("a", "b", "aniso_mat", "conditioning_nn", "correlation", "corrvecchia_knownCovparms", "distance_correlation", "kldiv", "matern_ns", "order_maxmin_correlation", "order_maxmin_correlation_old", "order_maxmin_euclidean", "simulation", "smoothness", "vecchia_specify_adjusted"), .packages='GPvecchia') %dopar% simulation(30^2, m = m, covparms = c(1))
parallel::stopCluster(cl)

kls.maxmin.euclidean    <- rep(NA, n.cand.m)
kls.maxmin.corr         <- rep(NA, n.cand.m)
kls.xcoord.euclidean    <- rep(NA, n.cand.m)
kls.ycoord.euclidean    <- rep(NA, n.cand.m)
for(i in 1:n.cand.m) {
  kls.maxmin.euclidean[i]    <- sim2[[i]]$kls[1]
  kls.maxmin.corr[i]         <- sim2[[i]]$kls[4]
  kls.xcoord.euclidean[i]    <- sim2[[i]]$kls[2]
  kls.ycoord.euclidean[i]    <- sim2[[i]]$kls[3]
}

vis.dat2 <- data.frame(kls.maxmin.euclidean, kls.maxmin.corr, kls.xcoord.euclidean, kls.ycoord.euclidean)
vis.dat2 <- vis.dat2[, order(colnames(vis.dat2))]
head(vis.dat2)

plot(cand.m, log10(vis.dat2$kls.maxmin.euclidean), type = "o", col = 1, lty = 1, lwd = 3,
     ylim = c(min(log10(vis.dat2)), max(log10(vis.dat2))), xlab = "m", ylab = "log10(KL)", main = NULL)
lines(cand.m, log10(vis.dat2$kls.maxmin.corr), type = "o", col = 2, lty = 2, lwd = 3)
lines(cand.m, log10(vis.dat2$kls.xcoord.euclidean), type = "o", col = 3, lty = 3, lwd = 3)
lines(cand.m, log10(vis.dat2$kls.ycoord.euclidean), type = "o", col = 4, lty = 4, lwd = 3)
legend("topright", legend=c("maxmin + Euclidean", "maxmin + correlation", "x-coord + Euclidean", "y-coord + Euclidean"), col=1:4, lty=1:4, lwd = 3, cex=1)

# save(sim2, cand.m, vis.dat2, file='2_corrvecchia/sim_nonstationarity_2.RData')
# load(file='2_corrvecchia/sim_nonstationarity_2.RData')