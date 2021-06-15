####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of approximations in various settings.
###
####################################################################################

rm(list = ls())
set.seed(05242021)

##

gc()

memory.size() # memory.limit()

##

library(vecchia) ; library(dplyr) ; library(tidyr) ; library(ggplot2) ; library(gridExtra) 

### anisotropic case ############################################################### (~ 45 min)

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.a  <- c(1, 5, 10, 15, 20, 25)

output.aniso <- parSim_aniso_knownCovparms(cand.m = cand.m, cand.a = cand.a, nsim = nsim, n = 30^2, d = 2, covmodel = cov_expo_aniso, covparms = c(1, 0.1), abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## save

save(nsim, cand.m, cand.a, output.aniso, file = "simout_aniso_05242021.RData")

## visualization

out     <- output.aniso
vdat1   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(a == 10) %>% select(-a)
vdat2   <- out$vars %>% left_join(out$kldiv, by = "index") %>% filter(m == 30) %>% select(-m)
vis     <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "a"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

ggplot2::ggsave("visout_aniso_05242021.pdf", vis, width = 15.2, height = 5.7)

### nonstationary case ############################################################# (~ 35 min)

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

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

output.nonst1 <- parSim_nonst_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, covmodel = cov_matern_ns, sigma = sigma, smoothness = smoothness, kernel = kernel, abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

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

output.nonst2 <- parSim_nonst_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, covmodel = cov_matern_ns, sigma = sigma, smoothness = smoothness, kernel = kernel, abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

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

output.nonst3 <- parSim_nonst_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, covmodel = cov_matern_ns, sigma = sigma, smoothness = smoothness, kernel = kernel, abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

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

output.nonst4 <- parSim_nonst_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, covmodel = cov_matern_ns, sigma = sigma, smoothness = smoothness, kernel = kernel, abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## save

save(nsim, cand.m, output.nonst1, output.nonst2, output.nonst3, output.nonst4, file = "simout_nonst_05242021.RData")

## visualization

out1    <- output.nonst1
out2    <- output.nonst4
vdat1   <- out1$vars %>% left_join(out1$kldiv, by = "index")
vdat2   <- out2$vars %>% left_join(out2$kldiv, by = "index")
vis     <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "m"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("E-MM + E-NN", "E-MM + C-NN", "C-MM + E-NN", "C-MM + C-NN", "X-ord + E-NN", "Y-ord + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "#FF7F00", "#FFFF33"), shape = c(18, 15, 17, 16, 8, 13))

ggplot2::ggsave("visout_nonst_05242021.pdf", vis, width = 15.2, height = 5.7)

### multivariate case ############################################################## (~ 6.5 hours)

## sim 1: bivariate

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.d  <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)

# randomly selected locations

output.biv.random <- parSim_mulv_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 2, covmodel = cov_latentDim_biv, covparms = c(1, 0.1), abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

# overlapped locations

output.biv.overlap <- parSim_mulv_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 2, covmodel = cov_latentDim_biv, covparms = c(1, 0.1), abs.corr = FALSE, method.locs = "overlap", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)

## save

save(nsim, cand.m, cand.d, output.biv.random, output.biv.overlap, file = "simout_biv_05242021.RData")

## visualization

out1      <- output.biv.random 
vdat1     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.rand  <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "d"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("visout_biv_rand_05242021.pdf", vis.rand, width = 15.2, height = 5.7)

out2      <- output.biv.overlap 
vdat1     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.over  <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "d"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("visout_biv_over_05242021.pdf", vis.over, width = 15.2, height = 5.7)

## sim 2: trivariate

# randomly selected locations

output.triv.random <- parSim_mulv_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 3, covmodel = cov_latentDim_triv, covparms = c(1, 0.1), abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

# overlapped locations

output.triv.overlap <- parSim_mulv_knownCovparms(cand.m = cand.m, cand.d = cand.d, nsim = nsim, n = 20^2, d = 2, p = 3, covmodel = cov_latentDim_triv, covparms = c(1, 0.1), abs.corr = FALSE, method.locs = "overlap", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)

## save

save(nsim, cand.m, cand.d, output.triv.random, output.triv.overlap, file = "simout_triv_05242021.RData")

## visualization

out1      <- output.triv.random 
vdat1     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out1$vars %>% left_join(out1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.rand  <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "d"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("visout_triv_rand_05242021.pdf", vis.rand, width = 15.2, height = 5.7)

out2      <- output.triv.overlap 
vdat1     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(d == 0.4) %>% select(-d)
vdat2     <- out2$vars %>% left_join(out2$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.over  <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "d"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("visout_triv_over_05242021.pdf", vis.over, width = 15.2, height = 5.7)

### spacetime case ################################################################# (~ 20 min)

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

output.sptm.gen1 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t = 1, covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0), abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
  
output.sptm.gen2 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t = 36, covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0), abs.corr = FALSE, method.locs = "monitoring", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
  
output.sptm.gen3 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t = 36, covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0), abs.corr = FALSE, method.locs = "grid", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
  
output.sptm.gen4 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t = 1, covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0), abs.corr = FALSE, method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## save

save(nsim, cand.m, output.sptm.gen1, output.sptm.gen2, output.sptm.gen3, output.sptm.gen4, file = "simout_sptm_05242021.RData")

## visualization

out11         <- output.sptm.gen1
out12         <- output.sptm.gen2
vdat1         <- out11$vars %>% left_join(out11$kldiv, by = "index")
vdat2         <- out12$vars %>% left_join(out12$kldiv, by = "index")
vis.rand      <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "m"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("visout_sptm_rand_05242021.pdf", vis.rand, width = 15.2, height = 5.7)

out13         <- output.sptm.gen3
out14         <- output.sptm.gen4
vdat1         <- out13$vars %>% left_join(out13$kldiv, by = "index") 
vdat2         <- out14$vars %>% left_join(out14$kldiv, by = "index")
vis.grid      <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "m"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("visout_sptm_grid_05242021.pdf", vis.grid, width = 15.2, height = 5.7)

### derivative case ################################################################

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.r  <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)

output.deriv.matern25.1d <- parSim_deriv_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = 20^2, d = 1, covmodel = cov_sqexpo_deriv, covtype = "matern", covparms = c(1, 2.5), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

output.deriv.matern25.2d <- parSim_deriv_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = 20^2, d = 2, covmodel = cov_sqexpo_deriv, covtype = "matern", covparms = c(1, 2.5), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

output.deriv.matern45.1d <- parSim_deriv_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = 20^2, d = 1, covmodel = cov_sqexpo_deriv, covtype = "matern", covparms = c(1, 4.5), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

output.deriv.matern45.2d <- parSim_deriv_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = 20^2, d = 2, covmodel = cov_sqexpo_deriv, covtype = "matern", covparms = c(1, 4.5), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

output.deriv.sqexpo.1d <- parSim_deriv_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = 20^2, d = 1, covmodel = cov_sqexpo_deriv, covtype = "sqexpo", covparms = c(1), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

output.deriv.sqexpo.2d <- parSim_deriv_knownCovparms(cand.m = cand.m, cand.r = cand.r, nsim = nsim, n = 20^2, d = 2, covmodel = cov_sqexpo_deriv, covtype = "sqexpo", covparms = c(1), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-4, ncores = NULL)

## save

save(nsim, cand.m, cand.r, output.deriv.matern25.1d, output.deriv.matern25.2d, output.deriv.matern45.1d, output.deriv.matern45.2d, output.deriv.sqexpo.1d, output.deriv.sqexpo.2d, file = "simout_deriv_05242021.RData")

## visualization

out.d1    <- output.deriv.sqexpo.1d
vdat1     <- out.d1$vars %>% left_join(out.d1$kldiv, by = "index") %>% filter(r == 0.1) %>% select(-r)
vdat2     <- out.d1$vars %>% left_join(out.d1$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.d1    <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "r"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("visout_deriv_1d_05242021.pdf", vis.d1, width = 15.2, height = 5.7)

out.d2    <- output.deriv.sqexpo.2d
vdat1     <- out.d2$vars %>% left_join(out.d2$kldiv, by = "index") %>% filter(r == 0.1) %>% select(-r)
vdat2     <- out.d2$vars %>% left_join(out.d2$kldiv, by = "index") %>% filter(m == 20) %>% select(-m)
vis.d2    <- vis_arrange_dio(vdat1 = vdat1, vdat2 = vdat2, ftn = log10, xlab = c("m", "r"), ylab = c("log10(KL)", "log10(KL)"), xlim = list(xl1 = NULL, xl2 = NULL), ylim = list(yl1 = NULL, yl2 = NULL), legend = c("S-E-MM + HH-E-NN", "S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "S-E-MM + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#FF7F00", "#E41A1C"), shape = c(18, 15, 17, 8, 16))

ggplot2::ggsave("visout_deriv_2d_05242021.pdf", vis.d2, width = 15.2, height = 5.7)

### wave case ###################################################################### (~ 5.5 hours)

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.l  <- c(0.1, 0.15, 0.2, 0.25, 0.30, 0.35, 0.4)

output.wave.ds.d1 <- parSim_wave_knownCovparms(cand.m = cand.m, cand.l = cand.l, nsim = nsim, n = 30^2, d = 1, covmodel = cov_wave, covtype = "Dampedsine", covparms = c(1), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)

output.wave.dc.d1 <- parSim_wave_knownCovparms(cand.m = cand.m, cand.l = cand.l, nsim = nsim, n = 30^2, d = 1, covmodel = cov_wave, covtype = "Dampedcosine", covparms = c(1, 0.1), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)

output.wave.bj.d1 <- parSim_wave_knownCovparms(cand.m = cand.m, cand.l = cand.l, nsim = nsim, n = 30^2, d = 1, covmodel = cov_wave, covtype = "BesselJ", covparms = c(1, 0), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)

output.wave.ds.d2 <- parSim_wave_knownCovparms(cand.m = cand.m, cand.l = cand.l, nsim = nsim, n = 30^2, d = 2, covmodel = cov_wave, covtype = "Dampedsine", covparms = c(1), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)

output.wave.dc.d2 <- parSim_wave_knownCovparms(cand.m = cand.m, cand.l = cand.l, nsim = nsim, n = 30^2, d = 2, covmodel = cov_wave, covtype = "Dampedcosine", covparms = c(1, 0.1), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-6, ncores = NULL)

output.wave.bj.d2 <- parSim_wave_knownCovparms(cand.m = cand.m, cand.l = cand.l, nsim = nsim, n = 30^2, d = 2, covmodel = cov_wave, covtype = "BesselJ", covparms = c(1, 0), abs.corr = TRUE, method.locs = "random", method.modify = "eigen-I", pivot = FALSE, tol = 1e-5, ncores = NULL)

## save

save(nsim, cand.m, cand.l, output.wave.ds.d1, output.wave.ds.d2, output.wave.dc.d1, output.wave.dc.d2, output.wave.bj.d1, output.wave.bj.d2, file = "simout_wave_05242021.RData")

## visualization

### end ############################################################################

Sys.time()