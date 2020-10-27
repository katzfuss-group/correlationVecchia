####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: .
###
####################################################################################

rm(list = ls())

# set.seed(09292020)
set.seed(10122020)

library(correlationVecchia)
library(GpGp)
library(dplyr)
library(foreach)
library(ggplot2)
library(gridExtra)

### spacetime case - visualization #############################################################################################

# out1  <- generate_gp_spacetime(nsim = 1, n = 100, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
# out2  <- generate_gp_spacetime(nsim = 1, n = 5^2, d = 2, t.len = 5, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
# out3  <- generate_gp_spacetime(nsim = 1, n = 5^2, d = 2, t.len = 5, method.locs = "all.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
# out4  <- generate_gp_spacetime(nsim = 1, n = 100, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
# 
# locs <- list()
# locs[[1]] <- out1$sim$sim1$locs %>% cbind(seq(100)) %>% as.data.frame()
# locs[[2]] <- out2$sim$sim1$locs %>% cbind(seq(5^3)) %>% as.data.frame()
# locs[[3]] <- out3$sim$sim1$locs %>% cbind(seq(5^3)) %>% as.data.frame()
# locs[[4]] <- out4$sim$sim1$locs %>% cbind(seq(100)) %>% as.data.frame()
# 
# for(i in 1:4) colnames(locs[[i]]) <- c("x", "y", "t", "index")
# 
# k <- 1
# p1 <- list()
# for(i in 1:4) {
#   
#   p1[[i]] <- locs[[i]] %>% ggplot() + geom_point(aes(x = x, y = y), size = 0.6) + xlab("x") + ylab("y")
#   k <- k + 1
#   
#   p1[[i+4]] <- locs[[i]] %>% ggplot() + geom_point(aes(x = index, y = t), size = 0.6) + xlab("index") + ylab("time")
#   k <- k + 1
# }
# 
# vis.all   <- arrangeGrob(grobs = p1, nrow = 2, ncol = 4)
# p1        <- grid.arrange(vis.all)

### spacetime case - generation #############################################################################################

# out1  <- generate_gp_spacetime(nsim = 1, n = 20^2, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
# out2  <- generate_gp_spacetime(nsim = 1, n = 20, d = 2, t.len = 4, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
# out3  <- generate_gp_spacetime(nsim = 1, n = 20, d = 2, t.len = 4, method.locs = "all.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
# out4  <- generate_gp_spacetime(nsim = 1, n = 20^2, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
# 
# locs <- list()
# locs[[1]] <- out1$sim$sim1$locs %>% as.matrix()
# locs[[2]] <- out2$sim$sim1$locs %>% as.matrix()
# locs[[3]] <- out3$sim$sim1$locs %>% as.matrix()
# locs[[4]] <- out4$sim$sim1$locs %>% as.matrix()

### spacetime case - process #############################################################################################

# covparms <- c(1, 0.1, 0.01, 3.5, 0) # var, range_space, range_time, nu, nugget
# 
# covmat <- list()
# covmat[[1]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[1]])
# covmat[[2]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[2]])
# covmat[[3]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[3]])
# covmat[[4]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[4]])
# 
# y <- list()
# y[[1]] <- t(chol(covmat[[1]])) %*% rnorm(n = nrow(covmat[[1]]))
# y[[2]] <- t(chol(covmat[[2]])) %*% rnorm(n = nrow(covmat[[2]]))
# y[[3]] <- t(chol(covmat[[3]])) %*% rnorm(n = nrow(covmat[[3]]))
# y[[4]] <- t(factorize(covmat[[1]], method = "eigen-I", tol = 1e-6)$covfactor) %*% rnorm(n = nrow(covmat[[4]])) # y[[4]] <- t(chol(covmat[[4]])) %*% rnorm(n = nrow(covmat[[4]]))

### spacetime case - fit #############################################################################################

# fit <- list()
# fit[[1]] <- fit_corrvecchia(y = y[[1]], inputs = locs[[1]], trend = "intercept", scale = "parms", covfun = "matern_spacetime", max.it = 100)
# fit[[2]] <- fit_corrvecchia(y = y[[2]], inputs = locs[[2]], trend = "intercept", scale = "parms", covfun = "matern_spacetime", max.it = 100)
# fit[[3]] <- fit_corrvecchia(y = y[[3]], inputs = locs[[3]], trend = "intercept", scale = "parms", covfun = "matern_spacetime", max.it = 100)
# fit[[4]] <- fit_corrvecchia(y = y[[4]], inputs = locs[[4]], trend = "intercept", scale = "parms", covfun = "matern_spacetime", max.it = 100)
# 
# fit[[1]]$covparms ; fit[[1]]$conv # randomly selected locations + measured at random times
# fit[[2]]$covparms ; fit[[2]]$conv # randomly selected locations + measured at regular times
# fit[[3]]$covparms ; fit[[3]]$conv # gridded locations + measured at regular times
# fit[[4]]$covparms ; fit[[4]]$conv # satellite-ish 

### spacetime case - baselines #############################################################################################

# fit_b1 <- list()
# fit_b1[[1]] <- fit_spacetime_baseline(y = y[[1]], inputs = locs[[1]], trend = "intercept", scale = "parms", baseline = 1, covfun = "matern_spacetime", max.it = 100)
# fit_b1[[2]] <- fit_spacetime_baseline(y = y[[2]], inputs = locs[[2]], trend = "intercept", scale = "parms", baseline = 1, covfun = "matern_spacetime", max.it = 100)
# fit_b1[[3]] <- fit_spacetime_baseline(y = y[[3]], inputs = locs[[3]], trend = "intercept", scale = "parms", baseline = 1, covfun = "matern_spacetime", max.it = 100)
# fit_b1[[4]] <- fit_spacetime_baseline(y = y[[4]], inputs = locs[[4]], trend = "intercept", scale = "parms", baseline = 1, covfun = "matern_spacetime", max.it = 100)
# 
# fit_b1[[1]]$covparms ; fit_b1[[1]]$conv # randomly selected locations + measured at random times
# fit_b1[[2]]$covparms ; fit_b1[[2]]$conv # randomly selected locations + measured at regular times
# fit_b1[[3]]$covparms ; fit_b1[[3]]$conv # gridded locations + measured at regular times
# fit_b1[[4]]$covparms ; fit_b1[[4]]$conv # satellite-ish 
# 
# fit_b2 <- list()
# fit_b2[[1]] <- fit_spacetime_baseline(y = y[[1]], inputs = locs[[1]], trend = "intercept", scale = "parms", baseline = 2, covfun = "matern_spacetime", max.it = 100)
# fit_b2[[2]] <- fit_spacetime_baseline(y = y[[2]], inputs = locs[[2]], trend = "intercept", scale = "parms", baseline = 2, covfun = "matern_spacetime", max.it = 100)
# fit_b2[[3]] <- fit_spacetime_baseline(y = y[[3]], inputs = locs[[3]], trend = "intercept", scale = "parms", baseline = 2, covfun = "matern_spacetime", max.it = 100)
# fit_b2[[4]] <- fit_spacetime_baseline(y = y[[4]], inputs = locs[[4]], trend = "intercept", scale = "parms", baseline = 2, covfun = "matern_spacetime", max.it = 100)
# 
# fit_b2[[1]]$covparms ; fit_b2[[1]]$conv # randomly selected locations + measured at random times
# fit_b2[[2]]$covparms ; fit_b2[[2]]$conv # randomly selected locations + measured at regular times
# fit_b2[[3]]$covparms ; fit_b2[[3]]$conv # gridded locations + measured at regular times
# fit_b2[[4]]$covparms ; fit_b2[[4]]$conv # satellite-ish 
# 
# fit_b3 <- list()
# fit_b3[[1]] <- fit_spacetime_baseline(y = y[[1]], inputs = locs[[1]], trend = "intercept", scale = "parms", baseline = 3, covfun = "matern_spacetime", max.it = 100)
# fit_b3[[2]] <- fit_spacetime_baseline(y = y[[2]], inputs = locs[[2]], trend = "intercept", scale = "parms", baseline = 3, covfun = "matern_spacetime", max.it = 100)
# fit_b3[[3]] <- fit_spacetime_baseline(y = y[[3]], inputs = locs[[3]], trend = "intercept", scale = "parms", baseline = 3, covfun = "matern_spacetime", max.it = 100)
# fit_b3[[4]] <- fit_spacetime_baseline(y = y[[4]], inputs = locs[[4]], trend = "intercept", scale = "parms", baseline = 3, covfun = "matern_spacetime", max.it = 100)
# 
# fit_b3[[1]]$covparms ; fit_b3[[1]]$conv # randomly selected locations + measured at random times
# fit_b3[[2]]$covparms ; fit_b3[[2]]$conv # randomly selected locations + measured at regular times
# fit_b3[[3]]$covparms ; fit_b3[[3]]$conv # gridded locations + measured at regular times
# fit_b3[[4]]$covparms ; fit_b3[[4]]$conv # satellite-ish 



### simulation with different ms #############################################################################################

nsim      <- 100
covparms  <- c(1, 0.1, 0.1, 3.5, 0) # var, range_space, range_time, nu, nugget
m         <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)

# covparms: not good
# c(1, 0.1, 0.3, 1.5, 0)
# c(1, 0.1, 0.3, 0.5, 0)
# c(1, 0.1, 0.2, 0.5, 0)
# c(1, 0.1, 0.1, 0.5, 0)
# c(1, 0.1, 0.5, 0.5, 0)
# c(1, 0.1, 1.0, 0.5, 0)
# c(1, 1.0, 1.0, 0.5, 0)
# c(1, 1.0, 0.1, 0.5, 0)
# c(1, 0.5, 0.1, 0.5, 0)
# c(1, 0.5, 0.1, 0.5, 0)

out2      <- generate_gp_spacetime(nsim = nsim, n = 25, d = 2, t.len = 30, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))

locs      <- list()
for(i in 1:nsim) locs[[i]] <- out2$sim[[i]]$locs %>% as.matrix()

covmat    <- list()
for(i in 1:nsim) covmat[[i]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[i]])

y         <- list()
for(i in 1:nsim) y[[i]] <- t(chol(covmat[[i]])) %*% rnorm(n = nrow(covmat[[i]])) # t(factorize(covmat[[i]], method = "eigen-I", tol = 1e-6)$covfactor)

# plot
fields::quilt.plot(as.numeric(locs[[1]][, 1]), as.numeric(locs[[1]][, 2]), as.numeric(y[[1]]))
plot(as.numeric(locs[[1]][, 3]), as.numeric(y[[1]]), type = 'o', col = as.factor(locs[[1]][, 1]))

# baseline_1
fit_b1 <- list()
for(j in 1:length(m)){
  fit_b1[[j]] <- list()
  for(i in 1:nsim) {
    fit_b1[[j]][[i]] <- fit_spacetime_baseline(y = y[[i]], inputs = locs[[i]], ms = m[j], trend = "intercept", scale = "parms", baseline = 1, covfun = "matern_spacetime", max.it = 100)
  }
}

perf_b1 <- rep(0, length(m))
for(j in 1:length(m)){
  for(i in 1:nsim){
    perf_b1[[j]] <- perf_b1[[j]] + sum((covparms - fit_b1[[j]][[i]]$covparms)^2) 
  }
} ; perf_b1 <- perf_b1 / nsim


# baseline_2
fit_b2 <- list()
for(j in 1:length(m)){
  fit_b2[[j]] <- list()
  for(i in 1:nsim) {
    fit_b2[[j]][[i]] <- fit_spacetime_baseline(y = y[[i]], inputs = locs[[i]], ms = m[j], trend = "intercept", scale = "parms", baseline = 2, covfun = "matern_spacetime", max.it = 100)
  }
}

perf_b2 <- rep(0, length(m))
for(j in 1:length(m)){
  for(i in 1:nsim){
    perf_b2[[j]] <- perf_b2[[j]] + sum((covparms - fit_b2[[j]][[i]]$covparms)^2) 
  }
} ; perf_b2 <- perf_b2 / nsim


# baseline_3
fit_b3 <- list()
for(j in 1:length(m)){
  fit_b3[[j]] <- list()
  for(i in 1:nsim) {
    fit_b3[[j]][[i]] <- fit_spacetime_baseline(y = y[[i]], inputs = locs[[i]], ms = m[j], trend = "intercept", scale = "parms", baseline = 3, covfun = "matern_spacetime", max.it = 100)
  }
}

perf_b3 <- rep(0, length(m))
for(j in 1:length(m)){
  for(i in 1:nsim){
    perf_b3[[j]] <- perf_b3[[j]] + sum((covparms - fit_b3[[j]][[i]]$covparms)^2) 
  }
} ; perf_b3 <- perf_b3 / nsim


# corrvecchia
fit_cc <- list()
for(j in 1:length(m)){
  fit_cc[[j]] <- list()
  for(i in 1:nsim) {
    fit_cc[[j]][[i]] <- fit_corrvecchia(y = y[[i]], inputs = locs[[i]], ms = m[j], trend = "intercept", scale = "parms", covfun = "matern_spacetime", max.it = 100)
  }
}

perf_cc <- rep(0, length(m))
for(j in 1:length(m)){
  for(i in 1:nsim){
    perf_cc[[j]] <- perf_cc[[j]] + sum((covparms - fit_cc[[j]][[i]]$covparms)^2) 
  }
} ; perf_cc <- perf_cc / nsim

save(nsim, covparms, m, out2, locs, covmat, y, fit_b1, fit_b2, fit_b3, fit_cc, perf_b1, perf_b2, perf_b3, perf_cc, file = "estims.RData")

plot(m, perf_b1, col = 1, type = 'o', ylim = range(c(perf_b1, perf_b2, perf_b3, perf_cc)), lwd = 2, ylab = "MSE(theta)", main = paste0("(variance, range_space, range_time, nu, nugget)", " = ", "(", paste0(covparms, collapse = ", "), ")"))
points(m, perf_b2, col = 3, type = 'o', lwd = 2)
points(m, perf_b3, col = 4, type = 'o', lwd = 2)
points(m, perf_cc, col = 2, type = 'o', lwd = 2)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c(1, 3, 4, 2), lty = 1, lwd = 2)

### KLdiv ###

kldiv_b1 <- list()
for(j in 1:length(m)) {
  
  kldiv_b1[[j]] <- list()
  for(i in 1:nsim) {
    
    out_specify <- baseline_1_spacetime_specify(locs = locs[[i]], m = m[j])
    
    indmat <- matrix(0, nrow = nrow(locs[[i]]), ncol = nrow(locs[[i]]))
    for(k in 1:nrow(locs[[i]])) indmat[k, out_specify$U.prep$revNNarray[k, ]] <- 1
    indmat <- indmat + t(indmat) ; diag(indmat) <- 1
    # fields::image.plot(indmat)
    
    # ?GpGp::matern_spacetime
    covmat_est  <- GpGp::matern_spacetime(covparms = fit_b1[[j]][[i]]$covparms, locs = locs[[i]])
    mu_est      <- rep(fit_b1[[j]][[i]]$betahat, nrow(locs[[i]]))
    
    revord      <- order(out_specify$ord)
    covmat_ord  <- covmat_est[out_specify$ord, out_specify$ord]
    
    U           <- GPvecchia::createU(out_specify, c(1), 0, covmodel = covmat_ord)$U
    covmat_est  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    # ?correlationVecchia::kldiv
    kldiv_b1[[j]][[i]] <- kldiv(covmat0 = covmat[[i]], covmat1 = covmat_est, mu0 = rep(0, length(mu_est)), mu1 = mu_est)
  }
}

kldiv_b2 <- list()
for(j in 1:length(m)) {
  
  kldiv_b2[[j]] <- list()
  for(i in 1:nsim) {
    
    out_specify <- baseline_2_spacetime_specify(locs = locs[[i]], m = m[j])
    
    indmat <- matrix(0, nrow = nrow(locs[[i]]), ncol = nrow(locs[[i]]))
    for(k in 1:nrow(locs[[i]])) indmat[k, out_specify$U.prep$revNNarray[k, ]] <- 1
    indmat <- indmat + t(indmat) ; diag(indmat) <- 1
    # fields::image.plot(indmat)
    
    # ?GpGp::matern_spacetime
    covmat_est  <- GpGp::matern_spacetime(covparms = fit_b2[[j]][[i]]$covparms, locs = locs[[i]])
    mu_est      <- rep(fit_b2[[j]][[i]]$betahat, nrow(locs[[i]]))
    
    revord      <- order(out_specify$ord)
    covmat_ord  <- covmat_est[out_specify$ord, out_specify$ord]
    
    U           <- GPvecchia::createU(out_specify, c(1), 0, covmodel = covmat_ord)$U
    covmat_est  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    # ?correlationVecchia::kldiv
    kldiv_b2[[j]][[i]] <- kldiv(covmat0 = covmat[[i]], covmat1 = covmat_est, mu0 = rep(0, length(mu_est)), mu1 = mu_est)
  }
}
  
kldiv_b3 <- list()
for(j in 1:length(m)) {
  
  kldiv_b3[[j]] <- list()
  for(i in 1:nsim) {
    
    out_specify <- baseline_3_spacetime_specify(locs = locs[[i]], m = m[j], covmodel = GpGp::matern_spacetime, covparms = fit_b3[[j]][[i]]$covparms)
    
    indmat <- matrix(0, nrow = nrow(locs[[i]]), ncol = nrow(locs[[i]]))
    for(k in 1:nrow(locs[[i]])) indmat[k, out_specify$U.prep$revNNarray[k, ]] <- 1
    indmat <- indmat + t(indmat) ; diag(indmat) <- 1
    # fields::image.plot(indmat)
    
    # ?GpGp::matern_spacetime
    covmat_est  <- GpGp::matern_spacetime(covparms = fit_b3[[j]][[i]]$covparms, locs = locs[[i]])
    mu_est      <- rep(fit_b3[[j]][[i]]$betahat, nrow(locs[[i]]))
    
    revord      <- order(out_specify$ord)
    covmat_ord  <- covmat_est[out_specify$ord, out_specify$ord]
    
    U           <- GPvecchia::createU(out_specify, c(1), 0, covmodel = covmat_ord)$U
    covmat_est  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    # ?correlationVecchia::kldiv
    kldiv_b3[[j]][[i]] <- kldiv(covmat0 = covmat[[i]], covmat1 = covmat_est, mu0 = rep(0, length(mu_est)), mu1 = mu_est)
  }
}

covmodel    <- function(locs, covparms) GpGp::matern_spacetime(covparms = covparms, locs = locs)

kldiv_cc <- list()
for(j in 1:length(m)) {
  
  kldiv_cc[[j]] <- list()
  for(i in 1:nsim) {
    
    out_specify <- corrvecchia_specify_knownCovparms_2(locs = locs[[i]], m = m[j], covmodel = covmodel, covparms = fit_b3[[j]][[i]]$covparms)
    
    indmat <- matrix(0, nrow = nrow(locs[[i]]), ncol = nrow(locs[[i]]))
    for(k in 1:nrow(locs[[i]])) indmat[k, out_specify$U.prep$revNNarray[k, ]] <- 1
    indmat <- indmat + t(indmat) ; diag(indmat) <- 1
    # fields::image.plot(indmat)
    
    # ?GpGp::matern_spacetime
    covmat_est  <- GpGp::matern_spacetime(covparms = fit_cc[[j]][[i]]$covparms, locs = locs[[i]])
    mu_est      <- rep(fit_cc[[j]][[i]]$betahat, nrow(locs[[i]]))
    
    revord      <- order(out_specify$ord)
    covmat_ord  <- covmat_est[out_specify$ord, out_specify$ord]
    
    U           <- GPvecchia::createU(out_specify, c(1), 0, covmodel = covmat_ord)$U
    covmat_est  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    cholfactor  <- t(chol(covmat_est/2 + t(covmat_est)/2, tol = 1e-8))
    
    # ?correlationVecchia::kldiv
    kldiv_cc[[j]][[i]]  <- kldiv(covmat0 = covmat[[i]], chol1 = cholfactor, mu0 = rep(0, length(mu_est)), mu1 = mu_est)
  }
}

meankld_b1 <- rep(0, length(m))
for(j in 1:length(m)) {
  for(i in 1:nsim) {
    meankld_b1[j] <- meankld_b1[j] + kldiv_b1[[j]][[i]]
  }
} ; meankld_b1 <- meankld_b1 / length(m)

meankld_b2 <- rep(0, length(m))
for(j in 1:length(m)) {
  for(i in 1:nsim) {
    meankld_b2[j] <- meankld_b2[j] + kldiv_b2[[j]][[i]]
  }
} ; meankld_b2 <- meankld_b2 / length(m)

meankld_b3 <- rep(0, length(m))
for(j in 1:length(m)) {
  for(i in 1:nsim) {
    meankld_b3[j] <- meankld_b3[j] + kldiv_b3[[j]][[i]]
  }
} ; meankld_b3 <- meankld_b3 / length(m)

meankld_cc <- rep(0, length(m))
for(j in 1:length(m)) {
  for(i in 1:nsim) {
    meankld_cc[j] <- meankld_cc[j] + kldiv_cc[[j]][[i]]
  }
} ; meankld_cc <- meankld_cc / length(m)

save(nsim, covparms, m, out2, locs, covmat, y, fit_b1, fit_b2, fit_b3, fit_cc, perf_b1, perf_b2, perf_b3, perf_cc, kldiv_b1, kldiv_b2, kldiv_b3, kldiv_cc, meankld_b1, meankld_b2, meankld_b3, meankld_cc, file = "estims_with_kld.RData")

mvis <- m[1:10] ; viskld_b1 <- meankld_b1[1:10] ; viskld_b2 <- meankld_b2[1:10] ; viskld_b3 <- meankld_b3[1:10] ; viskld_cc <- meankld_cc[1:10]

plot(mvis, log10(viskld_b1), col = 1, type = 'o', ylim = log10(range(c(viskld_b1, viskld_b2, viskld_b3, viskld_cc))), lwd = 2, ylab = "log10(KLdiv)", main = paste0("(variance, range_space, range_time, nu, nugget)", " = ", "(", paste0(covparms, collapse = ", "), ")"))
points(mvis, log10(viskld_b2), col = 3, type = 'o', lwd = 2)
points(mvis, log10(viskld_b3), col = 4, type = 'o', lwd = 2)
points(mvis, log10(viskld_cc), col = 2, type = 'o', lwd = 2)
legend("bottomleft", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c(1, 3, 4, 2), lty = 1, lwd = 2)


