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
set.seed(01032020)

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

save(out, file = "aniso.RData")

## visualization

vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(a == 10) %>% select(-a)
vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)
vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

ggplot2::ggsave("aniso.pdf", vis, width = 15.2, height = 5.7)

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

save(out1, out2, out3, out4, file = "nonst.RData")

## visualization

vdat1   <- out1$vars %>% left_join(out1$kldiv, by = "index")
vdat2   <- out4$vars %>% left_join(out4$kldiv, by = "index")
vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis     <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

ggplot2::ggsave("nonst.pdf", vis, width = 15.2, height = 5.7)

rm(list = ls())

Sys.time()

### Bivariate case #############################################################################################

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.d  <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

# legend  <- c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN") 
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C")
# shape   <- c(18, 15, 17, 8, 16)

## sim 1: randomly selected locations

out1 <- parallel_simulate_multivariate_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
save(out1, file = "temp1.RData")

## sim 2: overlapped locations

out2 <- parallel_simulate_multivariate_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 2, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = "overlap", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)
save(out2, file = "temp2.RData")

save(out1, out2, file = "multi.RData")

## visualization

vdat1     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis_ran   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis_ran   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

vdat1     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis_lap   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis_lap   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("multi_random.pdf", vis_ran, width = 15.2, height = 5.7)
ggplot2::ggsave("multi_overlap.pdf", vis_lap, width = 15.2, height = 5.7)

rm(list = ls())

Sys.time()

### Trivariate case #############################################################################################

Sys.time()

# nsim    <- 2
# cand.m  <- c(3, 6, 12, 15, 18, 24, 30, 36, 45)
# cand.d  <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.d  <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

## sim 1: randomly selected locations

out1 <- parallel_simulate_multivariate_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 3, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
save(out1, file = "temp1.RData")

## sim 2: overlapped locations

out2 <- parallel_simulate_multivariate_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 3, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1), method.locs = "overlap", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)
save(out2, file = "temp2.RData")

save(out1, out2, file = "trivariate.RData")


vdat1     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis_ran   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

vdat1     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis_lap   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis_lap   <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("triva_random.pdf", vis_ran, width = 15.2, height = 5.7)
ggplot2::ggsave("triva_overlap.pdf", vis_lap, width = 15.2, height = 5.7)

rm(list = ls())

Sys.time()

### spacetime case #############################################################################################

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

# legend  <- c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN") 
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")
# shape   <- c(18, 15, 17, 16)

## sim 1: covparms = c(1, 0.75, 50, 25)

out11 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25), method.locs = "all.random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out12 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25), method.locs = "space.random.time.grid", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out13 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25), method.locs = "all.grid", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out14 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25), method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out11, out12, out13, out14, file = "spti1.RData")

## sim 2: covparms = c(1, 0.75, 500, 2.5)

out21 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 500, 2.5), method.locs = "all.random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out22 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 500, 2.5), method.locs = "space.random.time.grid", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out23 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 500, 2.5), method.locs = "all.grid", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out24 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 500, 2.5), method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out21, out22, out23, out24, file = "spti2.RData")

## visualization

vdat1         <- out11$vars %>% left_join(out11$kldiv, by = "index")
vdat2         <- out12$vars %>% left_join(out12$kldiv, by = "index")
vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

vdat1         <- out13$vars %>% left_join(out13$kldiv, by = "index") 
vdat2         <- out14$vars %>% left_join(out14$kldiv, by = "index")
vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("spti1_sprand.pdf", vis_sprand, width = 15.2, height = 5.7)
ggplot2::ggsave("spti1_spgrid.pdf", vis_spgrid, width = 15.2, height = 5.7)

vdat1         <- out21$vars %>% left_join(out21$kldiv, by = "index")
vdat2         <- out22$vars %>% left_join(out22$kldiv, by = "index")
vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

vdat1         <- out23$vars %>% left_join(out23$kldiv, by = "index") 
vdat2         <- out24$vars %>% left_join(out24$kldiv, by = "index")
vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("spti2_sprand.pdf", vis_sprand, width = 15.2, height = 5.7)
ggplot2::ggsave("spti2_spgrid.pdf", vis_spgrid, width = 15.2, height = 5.7)

rm(list = ls())

Sys.time()

### spacetime case - 10052020 #############################################################################################

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

# legend  <- c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN") 
# color   <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")
# shape   <- c(18, 15, 17, 16)

## sim 1: covparms = c(1, 1, 0.1, 0.1)

gen11 <- generate_gp_spacetime(nsim = 3, n = 300, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_spacetime_expo, covparms = c(1, 1, 0.1, 0.1))

gen12 <- generate_gp_spacetime(nsim = 3, n = 25, d = 2, t.len = 10, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 1, 0.1, 0.1), method.modify = "eigen-I", pivot = FALSE, tol = 1e-3)

gen13 <- generate_gp_spacetime(nsim = 3, n = 25, d = 2, t.len = 10, method.locs = "all.grid", covmodel = cov_spacetime_expo, covparms = c(1, 1, 0.1, 0.1), method.modify = "eigen-I", pivot = FALSE, tol = 1e-3)

gen14 <- generate_gp_spacetime(nsim = 3, n = 300, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_spacetime_expo, covparms = c(1, 1, 0.1, 0.1))

# fields::quilt.plot(gen11$sim$sim1$locs[, 1], gen11$sim$sim1$locs[, 2], gen11$sim$sim1$y)
# # plot(gen11$sim$sim1$locs[, 3], gen11$sim$sim1$y)
# 
# fields::quilt.plot(gen12$sim$sim2$locs[, 1], gen12$sim$sim2$locs[, 2], gen12$sim$sim2$y)
# plot(gen12$sim$sim1$locs[, 3], gen12$sim$sim1$y, type = 'o', col = as.factor(gen12$sim$sim1$locs[, 2]))

out11 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 1, 0.1, 0.1), method.locs = "all.random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out12 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 1, 0.1, 0.1), method.locs = "space.random.time.grid", method.modify = "eigen-I", pivot = FALSE, tol = 1e-3, ncores = NULL)

out13 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 1, 0.1, 0.1), method.locs = "all.grid", method.modify = "eigen-I", pivot = FALSE, tol = 1e-3, ncores = NULL)

out14 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 1, 0.1, 0.1), method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out11, out12, out13, out14, file = "spti1_10052020.RData")

## sim 2: covparms = c(1, 1, 15, 1)

gen21 <- generate_gp_spacetime(nsim = 3, n = 300, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_spacetime_expo, covparms = c(1, 1, 15, 1))

gen22 <- generate_gp_spacetime(nsim = 3, n = 25, d = 2, t.len = 10, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 1, 15, 1), method.modify = "eigen-I", pivot = FALSE, tol = 1e-3)

gen23 <- generate_gp_spacetime(nsim = 3, n = 25, d = 2, t.len = 10, method.locs = "all.grid", covmodel = cov_spacetime_expo, covparms = c(1, 1, 15, 1), method.modify = "eigen-I", pivot = FALSE, tol = 1e-3)

gen24 <- generate_gp_spacetime(nsim = 3, n = 300, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_spacetime_expo, covparms = c(1, 1, 15, 1))

fields::quilt.plot(gen21$sim$sim1$locs[, 1], gen21$sim$sim1$locs[, 2], gen21$sim$sim1$y)
# plot(gen21$sim$sim1$locs[, 3], gen21$sim$sim1$y)

fields::quilt.plot(gen22$sim$sim2$locs[, 1], gen22$sim$sim2$locs[, 2], gen22$sim$sim2$y)
plot(gen22$sim$sim1$locs[, 3], gen22$sim$sim1$y, type = 'o', col = as.factor(gen22$sim$sim1$locs[, 2]))

out21 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 1, 15, 1), method.locs = "all.random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out22 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 1, 15, 1), method.locs = "space.random.time.grid", method.modify = "eigen-I", pivot = FALSE, tol = 1e-3, ncores = NULL)

out23 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 1, 15, 1), method.locs = "all.grid", method.modify = "eigen-I", pivot = FALSE, tol = 1e-3, ncores = NULL)

out24 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 1, 15, 1), method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out21, out22, out23, out24, file = "spti2_10052020.RData")

## sim 3: covparms = c(1, 0.1, 0.2, 0.5)

gen31 <- generate_gp_spacetime(nsim = 3, n = 300, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_spacetime_expo, covparms = c(1, 0.1, 0.2, 0.5))

gen32 <- generate_gp_spacetime(nsim = 3, n = 25, d = 2, t.len = 10, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.1, 0.2, 0.5), method.modify = "eigen-I", pivot = FALSE, tol = 1e-3)

gen33 <- generate_gp_spacetime(nsim = 3, n = 25, d = 2, t.len = 10, method.locs = "all.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.1, 0.2, 0.5), method.modify = "eigen-I", pivot = FALSE, tol = 1e-3)

gen34 <- generate_gp_spacetime(nsim = 3, n = 300, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_spacetime_expo, covparms = c(1, 0.1, 0.2, 0.5))

fields::quilt.plot(gen31$sim$sim1$locs[, 1], gen31$sim$sim1$locs[, 2], gen31$sim$sim1$y)
# plot(gen31$sim$sim1$locs[, 3], gen31$sim$sim1$y)

fields::quilt.plot(gen32$sim$sim2$locs[, 1], gen32$sim$sim2$locs[, 2], gen32$sim$sim2$y)
plot(gen32$sim$sim1$locs[, 3], gen32$sim$sim1$y, type = 'o', col = as.factor(gen32$sim$sim1$locs[, 2]))

out31 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.1, 0.2, 0.5), method.locs = "all.random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

out32 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.1, 0.2, 0.5), method.locs = "space.random.time.grid", method.modify = "eigen-I", pivot = FALSE, tol = 1e-3, ncores = NULL)

out33 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t.len = 36, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.1, 0.2, 0.5), method.locs = "all.grid", method.modify = "eigen-I", pivot = FALSE, tol = 1e-3, ncores = NULL)

out34 <- parallel_simulate_spacetime_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t.len = 1, corr.dist = "sqrt(1-abs(rho))", covmodel = cov_spacetime_expo, covparms = c(1, 0.1, 0.2, 0.5), method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

save(out31, out32, out33, out34, file = "spti3_10052020.RData")

## visualization

vdat1         <- out11$vars %>% left_join(out11$kldiv, by = "index")
vdat2         <- out12$vars %>% left_join(out12$kldiv, by = "index")
vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

vdat1         <- out13$vars %>% left_join(out13$kldiv, by = "index") 
vdat2         <- out14$vars %>% left_join(out14$kldiv, by = "index")
vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("spti1_sprand_10052020.pdf", vis_sprand, width = 15.2, height = 5.7)
ggplot2::ggsave("spti1_spgrid_10052020.pdf", vis_spgrid, width = 15.2, height = 5.7)

vdat1         <- out21$vars %>% left_join(out21$kldiv, by = "index")
vdat2         <- out22$vars %>% left_join(out22$kldiv, by = "index")
vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

vdat1         <- out23$vars %>% left_join(out23$kldiv, by = "index") 
vdat2         <- out24$vars %>% left_join(out24$kldiv, by = "index")
vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("spti2_sprand_10052020.pdf", vis_sprand, width = 15.2, height = 5.7)
ggplot2::ggsave("spti2_spgrid_10052020.pdf", vis_spgrid, width = 15.2, height = 5.7)

vdat1         <- out31$vars %>% left_join(out31$kldiv, by = "index")
vdat2         <- out32$vars %>% left_join(out32$kldiv, by = "index")
vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_sprand    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

vdat1         <- out33$vars %>% left_join(out33$kldiv, by = "index") 
vdat2         <- out34$vars %>% left_join(out34$kldiv, by = "index")
vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_spgrid    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("spti3_sprand_10052020.pdf", vis_sprand, width = 15.2, height = 5.7)
ggplot2::ggsave("spti3_spgrid_10052020.pdf", vis_spgrid, width = 15.2, height = 5.7)

rm(list = ls())

Sys.time()

### derivative case #############################################################################################

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.r  <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

#####

n <- 20^2
d <- 1
covmodel <- corr_derivative_matern_2.5_1d

out.d1 <- parallel_simulate_derivative_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = n, d = d, corr.dist = "sqrt(1-abs(rho))", covmodel = covmodel, covparms = c(1, NA), method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6)

vdat1     <- out.d1$vars %>% left_join(out.d1$kldiv, by = "index") %>% filter(r == 0.1) %>% select(-r)
vdat2     <- out.d1$vars %>% left_join(out.d1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.d1    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis.d1    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

save(out.d1, vis.d1, file = "deriv1.RData")

#####

n <- 20^2
d <- 2
covmodel <- corr_derivative_matern_2.5_2d

out.d2 <- parallel_simulate_derivative_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = n, d = d, corr.dist = "sqrt(1-abs(rho))", covmodel = covmodel, covparms = c(1, NA), method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6)

vdat1     <- out.d2$vars %>% left_join(out.d2$kldiv, by = "index") %>% filter(r == 0.1) %>% select(-r)
vdat2     <- out.d2$vars %>% left_join(out.d2$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.d2    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis.d2    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

save(out.d2, vis.d2, file = "deriv2.RData")

#####

n <- 20^2
d <- 1
covmodel <- corr_derivative_matern_4.5_1d

out.d3 <- parallel_simulate_derivative_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = n, d = d, corr.dist = "sqrt(1-abs(rho))", covmodel = covmodel, covparms = c(1, NA), method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6)

vdat1     <- out.d3$vars %>% left_join(out.d3$kldiv, by = "index") %>% filter(r == 0.1) %>% select(-r)
vdat2     <- out.d3$vars %>% left_join(out.d3$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.d3    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis.d3    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

save(out.d3, vis.d3, file = "deriv3.RData")

##### not work

n <- 20^2
d <- 2
covmodel <- corr_derivative_matern_4.5_2d

out.d4 <- parallel_simulate_derivative_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = n, d = d, corr.dist = "sqrt(1-abs(rho))", covmodel = covmodel, covparms = c(1, NA), method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4)

vdat1     <- out.d4$vars %>% left_join(out.d4$kldiv, by = "index") %>% filter(r == 0.1) %>% select(-r)
vdat2     <- out.d4$vars %>% left_join(out.d4$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.d4    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis.d4    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

save(out.d4, vis.d4, file = "deriv4.RData")

#####

n <- 20^2
d <- 1
covmodel <- corr_derivative_squared_expo_1d

out.d5 <- parallel_simulate_derivative_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = n, d = d, corr.dist = "sqrt(1-abs(rho))", covmodel = covmodel, covparms = c(1, NA), method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4)

vdat1     <- out.d5$vars %>% left_join(out.d5$kldiv, by = "index") %>% filter(r == 0.1) %>% select(-r)
vdat2     <- out.d5$vars %>% left_join(out.d5$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.d5    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis.d5    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

save(out.d5, vis.d5, file = "deriv5.RData")

#####

n <- 20^2
d <- 2
covmodel <- corr_derivative_squared_expo_2d

out.d6 <- parallel_simulate_derivative_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = n, d = d, corr.dist = "sqrt(1-abs(rho))", covmodel = covmodel, covparms = c(1, NA), method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4)

vdat1     <- out.d6$vars %>% left_join(out.d6$kldiv, by = "index") %>% filter(r == 0.1) %>% select(-r)
vdat2     <- out.d6$vars %>% left_join(out.d6$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.d6    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))
# vis.d6    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("Baseline 1", "Baseline 2", "Baseline 3", "Baseline 4", "C-Maxmin + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

save(out.d6, vis.d6, file = "deriv6.RData")

#####

save(out.d1, out.d2, out.d3, out.d4, out.d5, out.d6, file = "deriv.RData")

ggplot2::ggsave("deriv1.pdf", vis.d1, width = 15.2, height = 5.7)
ggplot2::ggsave("deriv2.pdf", vis.d2, width = 15.2, height = 5.7)
ggplot2::ggsave("deriv3.pdf", vis.d3, width = 15.2, height = 5.7)
ggplot2::ggsave("deriv4.pdf", vis.d4, width = 15.2, height = 5.7)
ggplot2::ggsave("deriv5.pdf", vis.d5, width = 15.2, height = 5.7)
ggplot2::ggsave("deriv6.pdf", vis.d6, width = 15.2, height = 5.7)

Sys.time()

### periodic case #############################################################################################

Sys.time()

nsim    <- 10
n       <- 30^2

##### d = 2 

cand.m <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.period <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.30, 0.35, 0.4)

out.p1 <- parallel_simulate_wave_knownCovparms(cand.m = cand.m, cand.period = cand.period, nsim = nsim, n = n, d = 2, covmodel = cov_wave, type = "Dampedsine", covparms = c(1), method.locs = 'random', corr.dist = "sqrt(1-abs(rho))", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

vdat1   <- out.p1$vars %>% left_join(out.p1$kldiv, by = "index") %>% filter(period == 0.2) %>% select(-period)
vdat2   <- out.p1$vars %>% left_join(out.p1$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)
vis.p1  <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis.p1  <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

save(out.p1, vis.p1, file = "wave1.RData")

##### d = 2

cand.m <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.period <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

out.p2 <- parallel_simulate_wave_knownCovparms(cand.m = cand.m, cand.period = cand.period, nsim = nsim, n = n, d = 2, covmodel = cov_wave, type = "Dampedcosine", covparms = c(1, 0.1), method.locs = 'random', corr.dist = "sqrt(1-abs(rho))", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

vdat1   <- out.p2$vars %>% left_join(out.p2$kldiv, by = "index") %>% filter(period == 0.2) %>% select(-period)
vdat2   <- out.p2$vars %>% left_join(out.p2$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)
vis.p2  <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis.p2  <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

save(out.p2, vis.p2, file = "wave2.RData")

##### d = 2 ## not work

cand.m <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.period <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

out.p3 <- parallel_simulate_wave_knownCovparms(cand.m = cand.m, cand.period = cand.period, nsim = nsim, n = n, d = 2, covmodel = cov_wave, type = "BesselJ", covparms = c(1, 0), method.locs = 'random', corr.dist = "sqrt(1-abs(rho))", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

vdat1   <- out.p3$vars %>% left_join(out.p3$kldiv, by = "index") %>% filter(period == 0.2) %>% select(-period)
vdat2   <- out.p3$vars %>% left_join(out.p3$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)
vis.p3  <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis.p3  <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

save(out.p3, vis.p3, file = "wave3.RData")

#####

save(out.p1, out.p2, out.p3, vis.p1, vis.p2, vis.p3, file = "wave.RData")

ggplot2::ggsave("wave1.pdf", vis.p1, width = 15.2, height = 5.7)
ggplot2::ggsave("wave2.pdf", vis.p2, width = 15.2, height = 5.7)
ggplot2::ggsave("wave3.pdf", vis.p3, width = 15.2, height = 5.7)

##### d = 1

cand.m <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.period <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.30, 0.35, 0.4)

out.p1d1 <- parallel_simulate_wave_knownCovparms(cand.m = cand.m, cand.period = cand.period, nsim = nsim, n = n, d = 1, covmodel = cov_wave, type = "Dampedsine", covparms = c(1), method.locs = 'random', corr.dist = "sqrt(1-abs(rho))", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

vdat1   <- out.p1d1$vars %>% left_join(out.p1d1$kldiv, by = "index") %>% filter(period == 0.2) %>% select(-period)
vdat2   <- out.p1d1$vars %>% left_join(out.p1d1$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)
vis.p1d1    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis.p1d1    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

save(out.p1d1, vis.p1d1, file = "wave1_d1.RData")

##### d = 1

cand.m <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.period <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

out.p2d1 <- parallel_simulate_wave_knownCovparms(cand.m = cand.m, cand.period = cand.period, nsim = nsim, n = n, d = 1, covmodel = cov_wave, type = "Dampedcosine", covparms = c(1, 0.1), method.locs = 'random', corr.dist = "sqrt(1-abs(rho))", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

vdat1   <- out.p2d1$vars %>% left_join(out.p2d1$kldiv, by = "index") %>% filter(period == 0.2) %>% select(-period)
vdat2   <- out.p2d1$vars %>% left_join(out.p2d1$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)
vis.p2d1    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis.p2d1    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

save(out.p2d1, vis.p2d1, file = "wave2_d1.RData")

##### d = 1

cand.m <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.period <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

out.p3d1 <- parallel_simulate_wave_knownCovparms(cand.m = cand.m, cand.period = cand.period, nsim = nsim, n = n, d = 1, covmodel = cov_wave, type = "BesselJ", covparms = c(1, 0), method.locs = 'random', corr.dist = "sqrt(1-abs(rho))", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

vdat1   <- out.p3d1$vars %>% left_join(out.p3d1$kldiv, by = "index") %>% filter(period == 0.2) %>% select(-period)
vdat2   <- out.p3d1$vars %>% left_join(out.p3d1$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)
vis.p3d1    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))
# vis.p3d1    <- vis_arrange(vdat1 = vdat1, vdat2 = vdat2, legend = c("E-Maxmin + E-NN", "E-Maxmin + C-NN", "C-Maxmin + E-NN", "C-Maxmin + C-NN", "X-Coord + E-NN", "Y-Coord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

save(out.p3d1, vis.p3d1, file = "wave3_d1.RData")

#####

save(out.p1d1, out.p2d1, out.p3d1, vis.p1d1, vis.p2d1, vis.p3d1, file = "wave_d1.RData")

ggplot2::ggsave("wave1_d1.pdf", vis.p1d1, width = 15.2, height = 5.7)
ggplot2::ggsave("wave2_d1.pdf", vis.p2d1, width = 15.2, height = 5.7)
ggplot2::ggsave("wave3_d1.pdf", vis.p3d1, width = 15.2, height = 5.7)

Sys.time()
