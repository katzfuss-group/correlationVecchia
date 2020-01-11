####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is the main script!
###
####################################################################################

rm(list = ls())

####################################################################################
### Initialization
####################################################################################

##### Import R packages #####

# GPvecchia
library(GPvecchia)

# Methodology
require(fields)
require(Bessel)
require(mvLSW)
require(Matrix)
require(GpGp)
require(Rfast)
require(scanstatistics)

# Parallel 
require(parallel)
require(foreach)
require(doParallel)

# Visualization
require(dplyr)
require(tidyr)
require(ggplot2)
require(RColorBrewer)
require(gridExtra)

##### Import R scripts #####

source("script/covmodels.R")
source("script/factorization.R")
source("script/vecchia.R")
source("script/simulation.R")
source("script/miscellany.R")
source("script/visualization.R")


####################################################################################
### Anisotropic case
####################################################################################

# out1 <- simulate_corrvecchia_knownCovparms(nsim = 3, n = 15^2, d = 2, p = 1,
#                                            same.locs = FALSE, m = 10,
#                                            covmodel = cov_expo_iso,
#                                            pivot = FALSE, method = NULL,
#                                            tol = .Machine$double.eps, verbose = TRUE,
#                                            covparms = c(1, 0.1))
# 
# rm(out1)


####################################################################################
### Nonstationary case
####################################################################################

# kernel <- function(loc) {
# 
#   d         <- length(loc)
# 
#   a         <- function(loc) 10
#   b         <- function(loc) 10
#   angle     <- function(loc) 0
# 
#   eta       <- angle(loc)
#   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
# 
#   range     <- c(a(loc)^(-2), b(loc)^(-2))
# 
#   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
# }
# 
# sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25
# 
# smoothness <- function(loc) 0.2 + 1.3 * loc[1]
# 
# out2 <- simulate_corrvecchia_knownCovparms(nsim = 3, n = 15^2, d = 2, p = 1,
#                                            same.locs = FALSE, m = 10,
#                                            covmodel = cov_matern_ns_bruteforce,
#                                            pivot = FALSE, method = "eigen-I",
#                                            tol = 1e-4, verbose = TRUE,
#                                            sigma = sigma,
#                                            smoothness = smoothness,
#                                            kernel = kernel)
#   
# rm(out2)


####################################################################################
### Multivariate case
####################################################################################

# out3 <- simulate_corrvecchia_knownCovparms(nsim = 3, n = 15^2, d = 2, p = 2,
#                                            same.locs = FALSE, m = 10,
#                                            covmodel = cov_multivariate_flexMatern,
#                                            pivot = FALSE, method = "eigen-I",
#                                            tol = 1e-4, verbose = TRUE,
#                                            sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2),
#                                            nu.mat = matrix(c(1.5, 1, 1, 0.5), 2, 2),
#                                            alpha.mat = matrix(0.01, 2, 2))
# 
# rm(out3)


####################################################################################
### Anisotropic cases
####################################################################################

# cand.m = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)
# cand.a = c(1, 5, 10, 15, 20, 25)
# nsim = 10
# n = 30^2
# d = 2
# covmodel = cov_expo_aniso
# covparms = c(1, 0.1)
# pivot = FALSE
# method = NULL
# tol = .Machine$double.eps

cand.m = c(10, 20, 30)
cand.a = c(1, 5, 10)
nsim = 3
n = 15^2
d = 2
covmodel = cov_expo_aniso
covparms = c(1, 0.1)
pivot = FALSE
method = NULL
tol = .Machine$double.eps

cand.all              <- expand.grid(cand.m, cand.a)
cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all)    <- c("index", "m", "a")

sim                   <- list()

no_cores              <- parallel::detectCores() - 2
cl                    <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)

sim <- foreach::foreach(m = cand.all$m, a = cand.all$a, 
                        .export = NULL, 
                        .packages = c("GPvecchia")) %dopar% simulate_corrvecchia_knownCovparms( nsim = nsim, n = n, d = d, p = 1, same.locs = FALSE, m = m, covmodel = covmodel, pivot = pivot, method = method, tol = tol, verbose = FALSE, covparms = covparms, a = c(a, rep(1, d - 1)) )

parallel::stopCluster(cl)

# KL divergence
n.approx        <- sim[[1]]$n.approx
kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
kls[, 1]        <- seq(nrow(cand.all))
for(k in 1:length(sim)) {
  kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
}
kls             <- as.data.frame(kls)
colnames(kls)   <- c("index", names(sim[[1]]$kls.average))

out <- list(vars = cand.all, kldiv = kls)

vdat1 <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(a == 10) %>% select(-a)
vdat2 <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)

vis <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2)


####################################################################################
### Nonstationary cases
####################################################################################

##### spatially varying smoothness #####

# cand.m = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)
# nsim = 10
# n = 30^2
# d = 2
# covmodel = cov_matern_ns_bruteforce
# pivot = FALSE
# method = NULL
# tol = .Machine$double.eps

cand.m = c(10, 20, 30)
nsim = 3
n = 15^2
d = 2
covmodel = cov_matern_ns_bruteforce
pivot = FALSE
method = NULL
tol = .Machine$double.eps

kernel <- function(loc) {

  d         <- length(loc)

  a         <- function(loc) 10
  b         <- function(loc) 10
  angle     <- function(loc) 0

  eta       <- angle(loc)
  rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)

  range     <- c(a(loc)^(-2), b(loc)^(-2))

  return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
}

sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25

smoothness <- function(loc) 0.2 + 1.3 * loc[1]

cand.all              <- data.frame(cand.m)
cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all)    <- c("index", "m")

sim                   <- list()

no_cores              <- parallel::detectCores() - 2
cl                    <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)

sim <- foreach::foreach(m = cand.all$m,
                        .export = NULL,
                        .packages = c("GPvecchia")) %dopar% simulate_corrvecchia_knownCovparms( nsim = nsim, n = n, d = d, p = 1, same.locs = FALSE, m = m, covmodel = covmodel, pivot = pivot, method = method, tol = tol, verbose = FALSE, sigma = sigma, smoothness = smoothness, kernel = kernel )

parallel::stopCluster(cl)

# KL divergence
n.approx        <- sim[[1]]$n.approx
kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
kls[, 1]        <- seq(nrow(cand.all))
for(k in 1:length(sim)) {
  kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
}
kls             <- as.data.frame(kls)
colnames(kls)   <- c("index", names(sim[[1]]$kls.average))

out1 <- list(vars = cand.all, kldiv = kls)

##### rotating anisotropic kernel #####

# cand.m = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)
# nsim = 10
# n = 30^2
# d = 2
# covmodel = cov_matern_ns_bruteforce
# pivot = FALSE
# method = NULL
# tol = .Machine$double.eps

cand.m = c(10, 20, 30)
nsim = 3
n = 15^2
d = 2
covmodel = cov_matern_ns_bruteforce
pivot = FALSE
method = NULL
tol = .Machine$double.eps

kernel <- function(loc) {

  d         <- length(loc)

  a         <- function(loc) 10 * 10
  b         <- function(loc) 10
  angle     <- function(loc) pi * loc[1] / 2

  eta       <- angle(loc)
  rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)

  range     <- c(a(loc)^(-2), b(loc)^(-2))

  return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
}

sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25

smoothness <- function(loc) 0.5

cand.all              <- data.frame(cand.m)
cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all)    <- c("index", "m")

sim                   <- list()

no_cores              <- parallel::detectCores() - 2
cl                    <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)

sim <- foreach::foreach(m = cand.all$m,
                        .export = NULL,
                        .packages = c("GPvecchia")) %dopar% simulate_corrvecchia_knownCovparms( nsim = nsim, n = n, d = d, p = 1, same.locs = FALSE, m = m, covmodel = covmodel, pivot = pivot, method = method, tol = tol, verbose = FALSE, sigma = sigma, smoothness = smoothness, kernel = kernel )

parallel::stopCluster(cl)

# KL divergence
n.approx        <- sim[[1]]$n.approx
kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
kls[, 1]        <- seq(nrow(cand.all))
for(k in 1:length(sim)) {
  kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
}
kls             <- as.data.frame(kls)
colnames(kls)   <- c("index", names(sim[[1]]$kls.average))

out2 <- list(vars = cand.all, kldiv = kls)

vdat1 <- out1$vars %>% left_join(out1$kldiv, by = "index")
vdat2 <- out2$vars %>% left_join(out2$kldiv, by = "index")

vis <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2)

####################################################################################
### Multivariate cases
####################################################################################

##### constant a #####

# cand.m = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)
# nsim = 10
# n = 30^2
# d = 2
# p = 2
# covmodel = cov_multivariate_flexMatern
# pivot = FALSE
# method = "eigen-I"
# tol = 1e-4

cand.m = c(10, 20, 30)
nsim = 3
n = 15^2
d = 2
p = 2
covmodel = cov_multivariate_flexMatern
pivot = FALSE
method = "eigen-I"
tol = 1e-4

sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2)
nu.mat = matrix(c(1.5, 1, 1, 0.5), 2, 2)
alpha.mat = matrix(0.01, 2, 2)

cand.all              <- data.frame(cand.m)
cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all)    <- c("index", "m")

sim                   <- list()

no_cores              <- parallel::detectCores() - 2
cl                    <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)

sim <- foreach::foreach(m = cand.all$m, # .export = NULL
                        .packages = c("vecchia", "GPvecchia")) %dopar% simulate_corrvecchia_knownCovparms( nsim = nsim, n = n, d = d, p = p, same.locs = FALSE, m = m, covmodel = covmodel, pivot = pivot, method = method, tol = tol, verbose = FALSE, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat )

parallel::stopCluster(cl)

# KL divergence
n.approx        <- sim[[1]]$n.approx
kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
kls[, 1]        <- seq(nrow(cand.all))
for(k in 1:length(sim)) {
  kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
}
kls             <- as.data.frame(kls)
colnames(kls)   <- c("index", names(sim[[1]]$kls.average))

out1 <- list(vars = cand.all, kldiv = kls)

##### non-constant a #####

# cand.m = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)
# nsim = 10
# n = 30^2
# d = 2
# p = 2
# covmodel = cov_multivariate_flexMatern
# pivot = FALSE
# method = "eigen-I"
# tol = 1e-4

cand.m = c(10, 20, 30)
nsim = 3
n = 15^2
d = 2
p = 2
covmodel = cov_multivariate_flexMatern
pivot = FALSE
method = "eigen-I"
tol = 1e-4

sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2)
nu.mat = matrix(c(1.5, 1, 1, 0.5), 2, 2)
alpha.mat = matrix(c(0.02, 0.02, 0.02, 0.01), 2, 2)

cand.all              <- data.frame(cand.m)
cand.all              <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all)    <- c("index", "m")

sim                   <- list()

no_cores              <- parallel::detectCores() - 2
cl                    <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)

sim <- foreach::foreach(m = cand.all$m, # .export = NULL
                        .packages = c("vecchia", "GPvecchia")) %dopar% simulate_corrvecchia_knownCovparms( nsim = nsim, n = n, d = d, p = p, same.locs = FALSE, m = m, covmodel = covmodel, pivot = pivot, method = method, tol = tol, verbose = FALSE, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat )

parallel::stopCluster(cl)

# KL divergence
n.approx        <- sim[[1]]$n.approx
kls             <- matrix(NA, nrow(cand.all), 1 + n.approx)
kls[, 1]        <- seq(nrow(cand.all))
for(k in 1:length(sim)) {
  kls[k, seq(n.approx) + 1]   <- sim[[k]]$kls.average
}
kls             <- as.data.frame(kls)
colnames(kls)   <- c("index", names(sim[[1]]$kls.average))

out2 <- list(vars = cand.all, kldiv = kls)

vdat1 <- out1$vars %>% left_join(out1$kldiv, by = "index")
vdat2 <- out2$vars %>% left_join(out2$kldiv, by = "index")

vis <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2)
