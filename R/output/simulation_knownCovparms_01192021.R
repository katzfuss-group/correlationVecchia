####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of Vecchia-like approximation in various settings.
###
####################################################################################

# use_build_ignore("R/output/simulation_knownCovparms.R", escape = TRUE)

# set.seed(03222020)
# set.seed(04282020)
# set.seed(01032020)
set.seed(01192021)

library(correlationVecchia)
library(dplyr)
library(foreach)
library(ggplot2)
library(gridExtra)

### anisotropic case #############################################################################################

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.a  <- c(1, 5, 10, 15, 20, 25)

# legend  <- c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN")
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33")
# shape   <- c(18, 15, 17, 16, 8, 13)

out     <- parallel_simulate_anisotropic_knownCovparms(cand.m = cand.m, cand.a = cand.a, nsim = nsim, n = 30^2, d = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_expo_aniso, covparms = c(1, 0.1), method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out, file = "aniso_01192021.RData")

## visualization

vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(a == 10) %>% select(-a)
vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)
vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

ggplot2::ggsave("aniso_01192021.pdf", vis, width = 15.2, height = 5.7)

rm(list = ls())

Sys.time()

### nonstationary case #############################################################################################

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

# legend  <- c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN")
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33")
# shape   <- c(18, 15, 17, 16, 8, 13)

## sim 1: smoothness

kernel <- function(loc) {
  
  d         <- length(loc)
  
  a         <- function(loc) 10 # please use your own a function for the first coordinate
  b         <- function(loc) 10 # please use your own a function for the second coordinate
  angle     <- function(loc) 0 # please use your own spatially varying rotation angle
  
  eta       <- angle(loc)
  rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
  
  range     <- c(a(loc)^(-2), b(loc)^(-2))
  
  return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
}

sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25

smoothness <- function(loc) 0.2 + 1.3 * loc[1]

out1 <- parallel_simulate_nonstationary_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_matern_ns_bruteforce, sigma = sigma, smoothness = smoothness, kernel = kernel, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## sim 2: range

kernel <- function(loc) {
  
  d         <- length(loc)
  
  a         <- function(loc) 10 * (1 + 10 * loc[1]) # please use your own a function for the first coordinate
  b         <- function(loc) 10 * (1 + 10 * loc[1]) # please use your own a function for the second coordinate
  angle     <- function(loc) 0 # please use your own spatially varying rotation angle
  
  eta       <- angle(loc)
  rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
  
  range     <- c(a(loc)^(-2), b(loc)^(-2))
  
  return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
}

sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25

smoothness <- function(loc) 0.5

out2 <- parallel_simulate_nonstationary_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_matern_ns_bruteforce, sigma = sigma, smoothness = smoothness, kernel = kernel, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## sim 3: anisotropy

kernel <- function(loc) {
  
  d         <- length(loc)
  
  a         <- function(loc) 10 * (1 + 10 * loc[1]) # please use your own a function for the first coordinate
  b         <- function(loc) 10 # please use your own a function for the second coordinate
  angle     <- function(loc) 0 # please use your own spatially varying rotation angle
  
  eta       <- angle(loc)
  rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
  
  range     <- c(a(loc)^(-2), b(loc)^(-2))
  
  return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
}

sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25

smoothness <- function(loc) 0.5

out3 <- parallel_simulate_nonstationary_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_matern_ns_bruteforce, sigma = sigma, smoothness = smoothness, kernel = kernel, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## sim 4: rotation

kernel <- function(loc) {
  
  d         <- length(loc)
  
  a         <- function(loc) 10 * 10 # please use your own a function for the first coordinate
  b         <- function(loc) 10 # please use your own a function for the second coordinate
  angle     <- function(loc) pi * loc[1] / 2 # please use your own spatially varying rotation angle
  
  eta       <- angle(loc)
  rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
  
  range     <- c(a(loc)^(-2), b(loc)^(-2))
  
  return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
}

sigma <- function(loc) determinant(kernel(loc), logarithm = FALSE)[[1]][1]^0.25

smoothness <- function(loc) 0.5

out4 <- parallel_simulate_nonstationary_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_matern_ns_bruteforce, sigma = sigma, smoothness = smoothness, kernel = kernel, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out1, out2, out3, out4, file = "nonst_01192021.RData")

## visualization

vdat1   <- out1$vars %>% left_join(out1$kldiv, by = "index")
vdat2   <- out4$vars %>% left_join(out4$kldiv, by = "index")
vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

ggplot2::ggsave("nonst_01192021.pdf", vis, width = 15.2, height = 5.7)

rm(list = ls())

Sys.time()

### Bivariate case #############################################################################################

Sys.time()

nsim    <- 2
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.d  <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

# legend  <- c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN") 
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C")
# shape   <- c(18, 15, 17, 8, 16)

## sim 1: randomly selected locations

out1 <- parallel_simulate_multivariate_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## sim 2: overlapped locations

out2 <- parallel_simulate_multivariate_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = "overlap", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)

save(out1, out2, file = "bivariate_01192021.RData")

## visualization

vdat1     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis_ran   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis_ran   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

vdat1     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis_lap   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis_lap   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("bivar_random_01192021.pdf", vis_ran, width = 15.2, height = 5.7)
ggplot2::ggsave("bivar_overlap_01192021.pdf", vis_lap, width = 15.2, height = 5.7)

rm(list = ls())

Sys.time()

### Trivariate case #############################################################################################

Sys.time()

# nsim    <- 2
# cand.m  <- c(3, 6, 12, 15, 18, 24, 30, 36, 45)
# cand.d  <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

nsim    <- 2
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.d  <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

## sim 1: randomly selected locations

out1 <- parallel_simulate_multivariate_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 3, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## sim 2: overlapped locations

out2 <- parallel_simulate_multivariate_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 3, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = "overlap", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)

save(out1, out2, file = "trivariate_01192021.RData")


vdat1     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis_ran   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

vdat1     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis_lap   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis_lap   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("trivar_random_01192021.pdf", vis_ran, width = 15.2, height = 5.7)
ggplot2::ggsave("trivar_overlap_01192021.pdf", vis_lap, width = 15.2, height = 5.7)

rm(list = ls())

Sys.time()

### spacetime case - 01192021 #############################################################################################

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

# legend  <- c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN") 
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")
# shape   <- c(18, 15, 17, 16)

## sim 1: covmodel = cov_matern_spacetime, covparms = c(1, 0.1, 0.2, 0.5)

gen11 <- generate_gp_spacetime(nsim = 2, n = 30^2, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = c(1, 0.1, 0.2, 0.5))

gen12 <- generate_gp_spacetime(nsim = 2, n = 25, d = 2, t.len = 36, method.locs = "space.random.time.grid", covmodel = cov_matern_spacetime, covparms = c(1, 0.1, 0.2, 0.5))

gen13 <- generate_gp_spacetime(nsim = 2, n = 25, d = 2, t.len = 36, method.locs = "all.grid", covmodel = cov_matern_spacetime, covparms = c(1, 0.1, 0.2, 0.5))

gen14 <- generate_gp_spacetime(nsim = 2, n = 30^2, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_matern_spacetime, covparms = c(1, 0.1, 0.2, 0.5))

# fields::quilt.plot(gen11$sim$sim1$locs[, 1], gen11$sim$sim1$locs[, 2], gen11$sim$sim1$y)
# # plot(gen11$sim$sim1$locs[, 3], gen11$sim$sim1$y)
# 
# fields::quilt.plot(gen12$sim$sim2$locs[, 1], gen12$sim$sim2$locs[, 2], gen12$sim$sim2$y)
# plot(gen12$sim$sim1$locs[, 3], gen12$sim$sim1$y, type = 'o', col = as.factor(gen12$sim$sim1$locs[, 2]))

out11 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_matern_spacetime, covparms = c(1, 0.1, 0.2, 0.5), method.locs = "all.random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out12 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_matern_spacetime, covparms = c(1, 0.1, 0.2, 0.5), method.locs = "space.random.time.grid", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out13 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_matern_spacetime, covparms = c(1, 0.1, 0.2, 0.5), method.locs = "all.grid", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out14 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_matern_spacetime, covparms = c(1, 0.1, 0.2, 0.5), method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out11, out12, out13, out14, file = "spti_01192021.RData")

vdat1         <- out11$vars %>% left_join(out11$kldiv, by = "index")
vdat2         <- out12$vars %>% left_join(out12$kldiv, by = "index")
vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

vdat1         <- out13$vars %>% left_join(out13$kldiv, by = "index") 
vdat2         <- out14$vars %>% left_join(out14$kldiv, by = "index")
vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("spti_sprand_01192021.pdf", vis_sprand, width = 15.2, height = 5.7)
ggplot2::ggsave("spti_spgrid_01192021.pdf", vis_spgrid, width = 15.2, height = 5.7)

