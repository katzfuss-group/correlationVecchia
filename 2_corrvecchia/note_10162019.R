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

## To visualize results
library(tidyr) ; library(ggplot2) ; library(gridExtra) ; library(RColorBrewer)

source("2_corrvecchia/vecchia_specify_adjusted.R")
source("2_corrvecchia/corrvecchia.R")
# source("2_corrvecchia/corrvecchia_with_ab_corr_dist.R")
source("2_corrvecchia/kldiv.R")

set.seed(10102019)


####################################################################
#### Wave covariance model (Risser MD, Calder CA (2015))
####################################################################

# cov.aniso <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])

# covparms = c(sigma, period)
covmodel_dampedsine <- function(locs, covparms) {
  h       <- fields::rdist(locs)
  
  ind     <- which(h < 1e-8, arr.ind = T)
  c       <- covparms[1] * sin(h / covparms[2])* (covparms[2] / h)
  c[ind]  <- covparms[1]
  
  c
}

# covparms = c(sigma, range, period)
# 0 <= covparms[2] <= covparms[3] for R^2
covmodel_dampedcosine <- function(locs, covparms) {
  h <- fields::rdist(locs)
  
  covparms[1] * exp(- h / covparms[2]) * cos(h / covparms[3])
}

# covparms = c(sigma, nu, period)
covmodel_besselJ <- function(locs, covparms) {
  h <- fields::rdist(locs)
  c <- gamma(covparms[2] + 1) * (2*covparms[3]/c(h))^covparms[2] * Bessel::BesselJ(c(h) / covparms[3], covparms[2])
  c[which(is.nan(c))] <- 1
  
  covparms[1] * matrix(c, nrow = nrow(h), ncol = ncol(h))
}

covmodel_wave <- function(locs, covparms, covm = "Dampedsine") {
  if(covm == "Dampedsine") {
    return(covmodel_dampedsine(locs, covparms))
  } else if(covm == "Dampedcosine") {
    return(covmodel_dampedcosine(locs, covparms))
  } else if(covm == "BesselJ") {
    return(covmodel_besselJ(locs, covparms))
  } else {
    stop("Select one of the following covariance models: Dampedsine, Dampedcosine, or besselJ.")
  }
}


####################################################################
#### pd function 
####################################################################

positive_def <- function(Sigma, tol){
  eig.decomp  <- eigen(Sigma)
  diagvec     <- ifelse(eig.decomp$values < tol, tol, eig.decomp$values)
  
  Sigma.modified <- eig.decomp$vectors %*% diag(diagvec) %*% t(eig.decomp$vectors)
  return(Sigma.modified)
}


####################################################################
#### visualization: 1-rho
####################################################################

# covparms = c(sigma, period)
covparms <- c(1, 1/10)

n         <- 10^2
locs      <- matrix(runif(n * 2, 0, 1), n, 2)
Sigma     <- covmodel_wave(locs = locs, covparms = covparms, covm = "Dampedsine")

Sigma.modified <- positive_def(Sigma, tol = 1e-5)

par(mfrow = c(2, 2))

y         <- as.numeric(t(chol(Sigma.modified)) %*% rnorm(n))
fields::quilt.plot(locs[,1], locs[,2], y, main = 'process')

h         <- fields::rdist(locs)
ord       <- order(unlist(h))
hord      <- unlist(h)[ord]
c         <- correlation(locs = locs, covmodel = Sigma.modified, covparms = covparms, def.dist = NULL)
cord      <- unlist(c)[ord]
plot(hord, cord, type = 'l', main = 'covariance model', xlab = NA, ylab = NA)
abline(h = 0, col = 'red')

d         <- distance_correlation(locs = locs, covmodel = Sigma.modified, covparms = covparms, def.dist = NULL)
fields::quilt.plot(locs[,1], locs[,2], d[, 1], main = 'distance')

approx <- corrvecchia_knownCovparms(locs = locs, m = 10, ordering = "maxmin", def.dist = NULL, ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = covparms)
ord.maxmin <- approx$ord 
k <- 10
each.vec <- rep(k, as.integer(n / k))
each.vec[length(each.vec)] <- each.vec[length(each.vec)] + n %% k
o <- rep(seq(as.integer(n / k)), each.vec) ; ord.o <- order(o) ; o <- paste0('gp', o)
vis.dat <- cbind(locs, o) 
vis.dat <- vis.dat[ord.o, ]
fac <- factor(vis.dat[, 3], levels = unique(vis.dat[, 3]))
plot(vis.dat[, 1], vis.dat[, 2], col = brewer.pal(as.integer(n / k), "Spectral")[fac], pch = 16, xlim = c(0, 1.2), ylim = c(0, 1), main = 'maxmin order', xlab = NA, ylab = NA)
abline(v = 1)
legend("topright", legend = levels(fac), col = brewer.pal(as.integer(n / k), "Spectral"), pch = 16, cex = 0.5)

par(mfrow = c(1, 1))


####################################################################
#### visualization: 1-|rho|
####################################################################

# covparms = c(sigma, period)
covparms <- c(1, 1/10)

n         <- 10^2
locs      <- matrix(runif(n * 2, 0, 1), n, 2)
Sigma     <- covmodel_wave(locs = locs, covparms = covparms, covm = "Dampedsine")

Sigma.modified <- positive_def(Sigma, tol = 1e-5)

par(mfrow = c(2, 2))

y         <- as.numeric(t(chol(Sigma.modified)) %*% rnorm(n))
fields::quilt.plot(locs[,1], locs[,2], y, main = 'process')

h         <- fields::rdist(locs)
ord       <- order(unlist(h))
hord      <- unlist(h)[ord]
c         <- correlation(locs = locs, covmodel = Sigma.modified, covparms = covparms, def.dist = NULL)
cord      <- unlist(c)[ord]
plot(hord, abs(cord), type = 'l', main = 'covariance model', xlab = NA, ylab = NA)
abline(h = 0, col = 'red')

d         <- distance_correlation(locs = locs, covmodel = Sigma.modified, covparms = covparms, def.dist = "abs")
fields::quilt.plot(locs[,1], locs[,2], d[, 1], main = 'distance')

approx <- corrvecchia_knownCovparms(locs = locs, m = 10, ordering = "maxmin", def.dist = "abs", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = covparms)
ord.maxmin <- approx$ord 
k <- 10
each.vec <- rep(k, as.integer(n / k))
each.vec[length(each.vec)] <- each.vec[length(each.vec)] + n %% k
o <- rep(seq(as.integer(n / k)), each.vec) ; ord.o <- order(o) ; o <- paste0('gp', o)
vis.dat <- cbind(locs, o) 
vis.dat <- vis.dat[ord.o, ]
fac <- factor(vis.dat[, 3], levels = unique(vis.dat[, 3]))
plot(vis.dat[, 1], vis.dat[, 2], col = brewer.pal(as.integer(n / k), "Spectral")[fac], pch = 16, xlim = c(0, 1.2), ylim = c(0, 1), main = 'maxmin order', xlab = NA, ylab = NA)
abline(v = 1)
legend("topright", legend = levels(fac), col = brewer.pal(as.integer(n / k), "Spectral"), pch = 16, cex = 0.5)

par(mfrow = c(1, 1))

