
gc()
rm(list = ls())

library(GPvecchia)

library(foreach)

## To visualize results
library(tidyr) ; library(ggplot2) ; library(gridExtra) ; library(RColorBrewer)

source("2_corrvecchia/vecchia_specify_adjusted.R")
# source("2_corrvecchia/corrvecchia.R")
source("2_corrvecchia/corrvecchia_with_ab_corr_dist.R")
source("2_corrvecchia/kldiv.R")

set.seed(10102019)

# covparms = c(sigma, range)
covparms <- c(1, 1/10)

covmodel_periodic <- function(locs, covparms) {
  h       <- fields::rdist(locs)
  
  ind     <- which(h < 1e-8, arr.ind = T)
  c       <- covparms[1] * sin(h / covparms[2])* (covparms[2] / h)
  c[ind]  <- covparms[1]
  
  c
}

positive_def <- function(Sigma, tol){
  eig.decomp  <- eigen(Sigma)
  diagvec     <- ifelse(eig.decomp$values < tol, tol, eig.decomp$values)
  
  Sigma.modified <- eig.decomp$vectors %*% diag(diagvec) %*% t(eig.decomp$vectors)
  return(Sigma.modified)
}

n         <- 10^2
locs      <- matrix(runif(n * 2, 0, 1), n, 2)
Sigma     <- covmodel_periodic(locs, covparms)

Sigma.modified <- positive_def(Sigma, tol = 1e-5)

par(mfrow = c(2, 2))

y         <- as.numeric(t(chol(Sigma.modified)) %*% rnorm(n))
fields::quilt.plot(locs[,1], locs[,2], y, main = 'process')

h         <- fields::rdist(locs)
c         <- correlation(locs, covmodel_periodic, covparms)
ord       <- order(unlist(h))
hord      <- unlist(h)[ord]
cord      <- unlist(c)[ord]
plot(hord, cord, type = 'l', main = 'covariance model', xlab = NULL, ylab = NULL)
abline(h = 0, col = 'red')

d         <- distance_correlation(locs, covmodel_periodic, covparms)
fields::quilt.plot(locs[,1], locs[,2], d[, 1], main = 'distance')

approx <- corrvecchia_knownCovparms(locs = locs, m = 10, ordering = 'maxmin', ordering.method = 'correlation', initial.pt = NULL, conditioning = 'NN', covmodel = covmodel_periodic, covparms = covparms)
ord.maxmin <- approx$ord 
k <- 10
each.vec <- rep(k, as.integer(n / k))
each.vec[length(each.vec)] <- each.vec[length(each.vec)] + n %% k
o <- rep(seq(as.integer(n / k)), each.vec) ; ord.o <- order(o) ; o <- paste0('gp', o)
vis.dat <- cbind(locs, o) 
vis.dat <- vis.dat[ord.o, ]
plot(vis.dat[, 1], vis.dat[, 2], col = factor(vis.dat[, 3], levels = unique(vis.dat[, 3])), pch = 16, xlim = c(0, 1.2), ylim = c(0, 1), main = 'maxmin order', xlab = NULL, ylab = NULL)
abline(v = 1)
legend("topright", legend = levels(factor(vis.dat[, 3], levels = unique(vis.dat[, 3]))), col = c(1:length(unique(vis.dat[, 3]))), pch = 16)
par(mfrow = c(1, 1))

