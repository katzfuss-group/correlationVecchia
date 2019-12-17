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
source("2_corrvecchia/kldiv.R")

set.seed(10102019)


####################################################################
#### Flexible multivariate Matern covariance model 
#### (Genton and Kleiber, 2015; Apanasovich, Genton, and Sun, 2012; Gneiting, Kleiber, and Schlather, 2010)
####################################################################

cov_flexMatern_bruteforce <- function(locs, sigma.mat, nu.mat, alpha.mat) {
  
  n <- nrow(locs)
  p <- nrow(sigma.mat)
  
  cov.mat <- matrix(NA, n * p, n * p)
  for(r in 1:n) {
    for(s in 1:n) {
      
      block.rs <- matrix(NA, p, p)
      for(i in 1:p) {
        for(j in 1:p) {
          block.rs[i, j] <- c_ij_fields_2(loc1 = locs[r, ], loc2 = locs[s, ], sigma = sigma.mat[i, j], nu = nu.mat[i, j], alpha = alpha.mat[i, j])
        }
      }
      cov.mat[seq(from = 1 + p * (r - 1), to = p * r, by = 1), seq(from = 1 + p * (s - 1), to = p * s, by = 1)] <- block.rs
      
    }
  }
  
  return(cov.mat)
}

c_ij_bruteforce <- function(loc1, loc2, sigma, nu, alpha) { # returns NaN when h = 0
  
  h <- sqrt(sum((loc1 - loc2)^2))
  
  return(  sigma * (1 / 2^(nu-1) / gamma(nu)) * (alpha * h)^(nu) * besselK(x = alpha * h, nu = nu)  )
}

c_ij_fields_1 <- function(loc1, loc2, sigma, nu, alpha) { # returns 1 when h = 0
  
  h <- sqrt(sum((loc1 - loc2)^2))
  
  return(  sigma * fields::Matern(h, range = 1/alpha, smoothness = nu)  )
}

c_ij_fields_2 <- function(loc1, loc2, sigma, nu, alpha) { # is equivalent to the c_ij_fields_1() function
  
  h <- sqrt(sum((loc1 - loc2)^2))
  
  return(  sigma * fields::Matern(h, alpha = alpha, nu = nu)  )
}

# locs <- matrix(runif(4), 2, 2)
# c_ij_bruteforce(loc1 = locs[1, ], loc2 = locs[2, ], sigma = 1, nu = 1, alpha = 1)
# c_ij_fields_1(loc1 = locs[1, ], loc2 = locs[2, ], sigma = 1, nu = 1, alpha = 1)
# c_ij_fields_2(loc1 = locs[1, ], loc2 = locs[2, ], sigma = 1, nu = 1, alpha = 1)
# 
# n <- 5 ; d <- 2 ; p <- 1 
# locs <- matrix(runif(n * d), n, d)
# sigma.mat <- diag(p) 
# nu.mat <- matrix(0.5, p, p)
# alpha.mat <- matrix(1, p, p)
# covmat <- cov_flexMatern_bruteforce(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat)
# covmat
# exp(-fields::rdist(locs))
# isSymmetric(covmat)
# matrixcalc::is.positive.definite(covmat)

cov_flexMatern_old <- function(locs, sigma.mat, nu.mat, alpha.mat) {
  
  n <- nrow(locs)
  p <- nrow(sigma.mat)
  
  cov.mat <- matrix(NA, n*p, n*p)
  for(i in 1:p) {
    for(j in 1:p) {
      ind.row     <- seq(from = i, by = p, length.out = n)
      ind.col     <- seq(from = j, by = p, length.out = n)
      
      cov.mat[ind.row, ind.col]   <- c_ij(locs, sigma = sigma.mat[i, j], nu = nu.mat[i, j], alpha = alpha.mat[i, j])
    }
  }
  
  return(cov.mat)
}

c_ij <- function(locs, sigma, nu, alpha) {
  return(  sigma * fields::Matern(fields::rdist(locs), alpha = alpha, nu = nu)  )
}

cov_flexMatern <- function(locs, locs2 = NULL, p = NULL, rho = NULL, sigma.mat = NULL, nu.vec = NULL, nu.mat = NULL, alpha = NULL, alpha.mat = NULL) {
  
  n         <- nrow(locs)
  m         <- ifelse(is.null(locs2), n, nrow(locs2))
  
  if( !is.null(sigma.mat) & !is.null(rho) ) {
    stop("Please use only one argument specifying the collocated covariance coefficient sigma.")
  } else if( is.null(sigma.mat) & is.null(rho) ) {
    stop("Please specify the collocated covariance coefficient sigma.")
  } else if( is.null(sigma.mat) & !is.null(rho) ) {
    if(is.null(p)) stop("Please specify the number of GPs.")
    sigma.mat <- rho^fields::rdist(seq(p))
  }
  
  p         <- nrow(sigma.mat)
  
  if( !is.null(nu.mat) & !is.null(nu.vec) ) {
    stop("Please use only one argument specifying the smoothness parameter nu.")
  } else if( is.null(nu.mat) & is.null(nu.vec) ) {
    stop("Please specify the smoothness parameter nu.")
  } else if( is.null(nu.mat) & !is.null(nu.vec) ) {
    nu.mat <- (matrix(nu.vec, p, p, byrow = T) + matrix(nu.vec, p, p, byrow = F)) / 2
  }
  
  if( !is.null(alpha.mat) & !is.null(alpha) ) {
    stop("Please use only one argument specifying the scale parameter alpha.")
  } else if( is.null(alpha.mat) & is.null(alpha) ) {
    stop("Please specify the scale parameter alpha.")
  } else if( is.null(alpha.mat) & !is.null(alpha) ) {
    alpha.mat <- matrix(alpha, p, p)
  }
  
  cov.mat   <- matrix(NA, n * p, m * p)
  for(i in 1:p) {
    for(j in 1:p) {
      ind.row   <- seq(from = i, by = p, length.out = n)
      ind.col   <- seq(from = j, by = p, length.out = m)
      
      cov.mat[ind.row, ind.col] <- sigma.mat[i, j] * fields::Matern(fields::rdist(x1 = locs, x2 = locs2), alpha = alpha.mat[i, j], nu = nu.mat[i, j])
    }
  }
  
  return(cov.mat)
}

# n <- 30 ; d <- 2 ; p <- 5
# locs <- matrix(runif(n * d), n, d)
# sigma.mat <- diag(p) ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)
# 
# covmat_bruteforce <- cov_flexMatern_bruteforce(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat)
# covmat_old <- cov_flexMatern_old(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat)
# covmat_new <- cov_flexMatern(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat)
# covmat_rho <- cov_flexMatern(locs = locs, p = 5, rho = 0, sigma.mat = NULL, nu.mat = nu.mat, alpha.mat = alpha.mat)
# sqrt(sum((covmat_old - covmat_bruteforce)^2))
# sqrt(sum((covmat_new - covmat_bruteforce)^2))
# sqrt(sum((covmat_rho - covmat_bruteforce)^2))
# 
# microbenchmark::microbenchmark(
#   cov_flexMatern_bruteforce(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat),
#   cov_flexMatern_old(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat),
#   cov_flexMatern(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat),
#   times = 10
# )
# 
# n1 <- 3 ; n2 <- 4 ; d <- 2 ; p <- 2
# locs1 <- matrix(runif(n1 * d), n1, d) ; locs2 <- matrix(runif(n2 * d), n2, d) ; locs <- rbind(locs1, locs2)
# sigma.mat <- diag(0.5, 2) + 0.5 ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)
# 
# round(cov_flexMatern(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat), 3)[seq(1, n1*p), seq(n1*p+1, (n1+n2)*p)]
# round(cov_flexMatern(locs = locs1, locs2 = locs2, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat), 3)


####################################################################
#### how to generate and plot it
####################################################################

# ### full bivariate Matern model (p = 2)
# n <- 20^2 ; p <- 2
# 
# locs        <- matrix(runif(n*2, 0, 1), n, 2)
# sigma.mat   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
# nu.vec      <- c(1.5, 0.5)
# alpha       <- 1
# 
# Sigma       <- cov_flexMatern(locs = locs, locs2 = NULL, sigma.mat = sigma.mat, nu.vec = nu.vec, nu.mat = NULL, alpha = 1, alpha.mat = NULL)
# # min(eigen(Sigma)$values)
# 
# y.full      <- as.numeric(t(chol(Sigma)) %*% rnorm(n * p))
# y1          <- y.full[seq(from = 1, to = n*p, by = p)]
# y2          <- y.full[seq(from = 2, to = n*p, by = p)]
# 
# par(mfrow = c(1, 2))
# fields::quilt.plot(locs[, 1], locs[, 2], y1, main = 'GP 1') # smooth
# fields::quilt.plot(locs[, 1], locs[, 2], y2, main = 'GP 2') # rough
# par(mfrow = c(1, 1))
# 
# ### parsimonious multivariate Matern covariance model (p = 4)
# n <- 20^2 ; p <- 4
# 
# locs        <- matrix(runif(n*2, 0, 1), n, 2)
# rho         <- 0.5
# nu.vec      <- c(1.5, 1, 0.5, 0.1)
# alpha       <- 1
# 
# Sigma       <- cov_flexMatern(locs = locs, locs2 = NULL, p = p, rho = rho, nu.vec = nu.vec, nu.mat = NULL, alpha = 1, alpha.mat = NULL)
# # min(eigen(Sigma)$values)
# 
# y.full      <- as.numeric(t(chol(Sigma)) %*% rnorm(n * p))
# y1          <- y.full[seq(from = 1, to = n*p, by = p)]
# y2          <- y.full[seq(from = 2, to = n*p, by = p)]
# y3          <- y.full[seq(from = 3, to = n*p, by = p)]
# y4          <- y.full[seq(from = 4, to = n*p, by = p)]
# 
# par(mfrow = c(2, 2))
# fields::quilt.plot(locs[, 1], locs[, 2], y1, main = 'GP 1')
# fields::quilt.plot(locs[, 1], locs[, 2], y2, main = 'GP 2') 
# fields::quilt.plot(locs[, 1], locs[, 2], y3, main = 'GP 3') 
# fields::quilt.plot(locs[, 1], locs[, 2], y4, main = 'GP 4') 
# par(mfrow = c(1, 1))


####################################################################
#### simulation function
####################################################################

positive_def <- function(Sigma, tol){
  eig.decomp  <- eigen(Sigma)
  diagvec     <- ifelse(eig.decomp$values < tol, tol, eig.decomp$values)

  Sigma.modified <- eig.decomp$vectors %*% diag(diagvec) %*% t(eig.decomp$vectors)
  return(Sigma.modified)
}

simulation <- function(n = 15^2, m = 10, p = 2, rho = 0.5, sigma.mat = NULL, nu.vec = 0.5, nu.mat = NULL, alpha = 1, alpha.mat = NULL, tol = 1e-6) {

  locs      <- matrix(runif(n * 2, 0, 1), n, 2)
  Sigma     <- cov_flexMatern(locs = locs, locs2 = NULL, p = p, rho = rho, sigma.mat = sigma.mat, nu.vec = nu.vec, nu.mat = nu.mat, alpha = alpha, alpha.mat = alpha.mat)

  Sigma.modified <- positive_def(Sigma, tol)

  y         <- as.numeric(t(chol(Sigma.modified)) %*% rnorm(n * p))

  # locs
  locs <- locs[rep(1:n, times = p), ]
  
  ### specify vecchia approximations
  approx <- list()
  
  # standard vecchia with maxmin ordering
  approx[[1]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
  # standard vecchia with x coord ordering
  approx[[2]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "coord", which.coord = 1, cond.yz='y', conditioning = "NN")
  # standard vecchia with y coord ordering
  approx[[3]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "coord", which.coord = 2, cond.yz='y', conditioning = "NN")
  # euclidean-based ordering + euclidean-based NN conditioning
  approx[[4]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma.modified, covparms = c(1))
  # euclidean-based ordering + correlation-based NN conditioning
  approx[[5]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = c(1))
  # correlation-based ordering + euclidean-based NN conditioning
  approx[[6]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma.modified, covparms = c(1))
  # correlation-based ordering + correlation-based NN conditioning
  approx[[7]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = c(1))
  # x-coordinate-based ordering + correlation-based NN conditioning
  approx[[8]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(1), def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = c(1))
  # y-coordinate-based ordering + correlation-based NN conditioning
  approx[[9]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(2), def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = c(1))
  # (x+y)-coordinate-based ordering + correlation-based NN conditioning
  approx[[10]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(1, 2), def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = c(1))
  # (x+y)-coordinate-based ordering + euclidean-based NN conditioning
  approx[[11]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(1, 2), def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma.modified, covparms = c(1))
  # standard vecchia with coord ordering
  approx[[12]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "coord", which.coord = c(1, 2), cond.yz='y', conditioning = "NN")
  
  ### compute approximate covariance matrices
  n.approx    <- length(approx)
  Sigma.hat   <- list()
  kls         <- rep(NA, n.approx)
  for(i in 1:n.approx){
    
    # true cov in appropriate ordering
    Sigma.ord       <- cov_flexMatern(locs = approx[[i]]$locsord, locs2 = NULL, p = p, rho = rho, sigma.mat = sigma.mat, nu.vec = nu.vec, nu.mat = nu.mat, alpha = alpha, alpha.mat = alpha.mat)

    Sigma.ord.modified <- positive_def(Sigma.ord, tol)
    
    U               <- createU(approx[[i]], c(1, 1, 1), 0, covmodel = Sigma.ord.modified)$U
    revord          <- order(approx[[i]]$ord)
    Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    kls[i]          <- kldiv(Sigma.modified, Sigma.hat[[i]])
  }
  
  result                  <- list()
  result$n                <- n
  result$m                <- m
  result$p                <- p
  result$rho              <- rho
  result$sigma.mat        <- sigma.mat
  result$nu.vec           <- nu.vec
  result$nu.mat           <- nu.mat
  result$alpha            <- alpha
  result$alpha.mat        <- alpha.mat
  result$tol              <- tol
  
  result$locs             <- locs[1:n, ]
  result$Sigma            <- Sigma
  result$Sigma.modified   <- Sigma.modified
  result$approx           <- approx
  result$Sigma.hat        <- Sigma.hat
  result$kls              <- kls
  
  return(result)  
}

out <- simulation(n = 15^2, m = 10, p = 2, rho = 0.5, sigma.mat = NULL, nu.vec = 0.5, nu.mat = NULL, alpha = 1, alpha.mat = NULL, tol = 1e-6)
plot(log(out$kls), col = c(3, 1, 1, 3, 3, 2, 2, 1, 1, 1, 1, 1), lwd = 2, main = 'Log-scale KL divergence', ylab = 'KL divergence', xlab = 'method')
# green = Euclidean-based ordering (does not make sense but work much better)
# red = correlation-based ordering
# black = others (does not make sense)

# 01 = standard vecchia with maxmin ordering
# 02 = standard vecchia with x coord ordering
# 03 = standard vecchia with y coord ordering
# 04 = euclidean-based ordering + euclidean-based NN conditioning
# 05 = euclidean-based ordering + correlation-based NN conditioning
# 06 = correlation-based ordering + euclidean-based NN conditioning
# 07 = correlation-based ordering + correlation-based NN conditioning
# 08 = x-coordinate-based ordering + correlation-based NN conditioning
# 09 = y-coordinate-based ordering + correlation-based NN conditioning
# 10 = (x+y)-coordinate-based ordering + correlation-based NN conditioning
# 11 = (x+y)-coordinate-based ordering + euclidean-based NN conditioning
# 12 = standard vecchia with coord ordering


####################################################################
#### visualization
####################################################################


####################################################################
#### simulation 1:
####################################################################


####################################################################
#### simulation 2:
####################################################################
