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

## To visualize results
library(tidyr) ; library(ggplot2) ; library(gridExtra) ; library(RColorBrewer)

source("1_Pilot_Study/2_vecchia_specify_adjusted.R")
source("2_corrvecchia/corrvecchia.R")
source("2_corrvecchia/kldiv.R")
source("1_Pilot_Study/5_visualization.R")

set.seed(10102019)

####################################################################
#### Nonstaionary Matern covariance model (Risser MD, Calder CA (2015))
####################################################################

covparms <- c(1)

a <- function(loc) 10

# spatially-varying standard deviation
sigma <- function(loc) determinant(aniso_mat(loc), logarithm = F)[[1]][1]^0.25
# spatially-varying local anisotropy (controlling both the range and direction of dependence)
aniso_mat<- function(loc) {
  
  eta <- 0
  rot.mat <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
  
  range <- c(a(loc)^(-2), 1)
  diag.mat <- diag(range, nrow = length(loc))
  
  aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
  
  return(aniso.mat)
}
# matern's smoothness
smoothness <- function(loc) 0.5

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

# n         <- 15^2
# locs      <- matrix(runif(n * 2, 0, 1), n, 2)
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
#### simulation 1
####################################################################

cand.m    <- c(5, 20, 45) ; n.cand.m <- length(cand.m)
cand.scale  <- c(1, 10, 25) ; n.cand.scale <- length(cand.scale)
sim1      <- list()

cand.all <- expand.grid(cand.m, cand.scale)
cand.all <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all) <- c('index', 'm', 'scale')
n.cand.all <- n.cand.m * n.cand.scale

for(i in cand.all$index){
  a <- function(loc) cand.all[i, 3]
  sim1[[i]] <- simulation(10^2, m = cand.all[i, 2], covparms = c(1))
}

kls.maxmin.euclidean    <- rep(NA, n.cand.all)
kls.maxmin.corr         <- rep(NA, n.cand.all)
kls.xcoord.euclidean    <- rep(NA, n.cand.all)
kls.ycoord.euclidean    <- rep(NA, n.cand.all)
for(i in 1:n.cand.all) {
  kls.maxmin.euclidean[i]    <- sim1[[i]]$kls[1]
  kls.maxmin.corr[i]         <- sim1[[i]]$kls[4]
  kls.xcoord.euclidean[i]    <- sim1[[i]]$kls[2]
  kls.ycoord.euclidean[i]    <- sim1[[i]]$kls[3]
}

set.scale <- 10
ind <- cand.all$scale == set.scale
vis.dat1 <- data.frame(kls.maxmin.euclidean[ind], kls.maxmin.corr[ind], kls.xcoord.euclidean[ind], kls.ycoord.euclidean[ind])
vis.dat1 <- vis.dat1[, order(colnames(vis.dat1))]
vis.dat1 <- cbind(rep(cand.m, times = 4), tidyr::gather(vis.dat1))
colnames(vis.dat1) <- c("m", "method", "KL")
head(vis.dat1)

set.m <- 20
ind <- cand.all$m == set.m
vis.dat2 <- data.frame(kls.maxmin.euclidean[ind], kls.maxmin.corr[ind], kls.xcoord.euclidean[ind], kls.ycoord.euclidean[ind])
vis.dat2 <- vis.dat2[, order(colnames(vis.dat2))]
vis.dat2 <- cbind(rep(cand.scale, times = 4), tidyr::gather(vis.dat2))
colnames(vis.dat2) <- c("scale", "method", "KL")
head(vis.dat2)

kls.legend <- c("Correlation + Maxmin     ", "Euclidean + Maxmin     ", "Euclidean + x-coord     ", "Euclidean + y-coord")
vis_arrange(vdat1 = vis.dat1, vdat2 = vis.dat2, combined.legend = kls.legend, color.pal = brewer.pal(4, "Set1"), alpha.value = 0.7, size.legend = 16, size.lab = 16, size.text = 12)





