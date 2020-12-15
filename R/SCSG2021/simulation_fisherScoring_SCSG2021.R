####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: 
###
####################################################################################

rm(list = ls())
set.seed(12122020)

library(correlationVecchia) ; library(GpGp)

library(dplyr) ; library(ggplot2) ; library(gridExtra)

library(foreach)

### setting #############################################################################################

nsim      <- 100
m         <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

covparms  <- c(1, 0.1, 0.2, 0.5, 0) # var, range_space, range_time, nu, nugget

### ncores #############################################################################################

no_cores  <- parallel::detectCores() - 2

# sim                   <- list()
# cl                    <- parallel::makeCluster(no_cores)
# 
# doParallel::registerDoParallel(cl)
# sim <- foreach::foreach(m = cand.all$m, a = cand.all$a, .packages = c("correlationVecchia", "GPvecchia")) %dopar% simulate_univariate_knownCovparms(nsim = nsim, n = n, d = d, m = m, method.locs = method.locs, corr.dist = corr.dist, method.modify = method.modify, pivot = pivot, tol = tol, verbose = FALSE, covmodel = covmodel, covparms = covparms, a = c(a, rep(1, d - 1)) )
# parallel::stopCluster(cl)

### simulation with different ms #############################################################################################

Sys.time()

out2      <- generate_gp_spacetime(nsim = nsim, n = 25, d = 2, t.len = 36, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
# out2      <- generate_gp_spacetime(nsim = nsim, n = 900, d = 2, t.len = 1, method.locs = "satellite",covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))

locs      <- list()
for(i in 1:nsim) locs[[i]] <- out2$sim[[i]]$locs %>% as.matrix()

covmat    <- list()
# for(i in 1:nsim) covmat[[i]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[i]])
cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
covmat    <- foreach::foreach(i = seq(nsim), .packages = c("GPvecchia", "GpGp")) %dopar% matern_spacetime(covparms = covparms, locs = locs[[i]])
parallel::stopCluster(cl)

y         <- list()
# for(i in 1:nsim) y[[i]] <- t(chol(covmat[[i]])) %*% rnorm(n = nrow(covmat[[i]])) # t(factorize(covmat[[i]], method = "eigen-I", tol = 1e-6)$covfactor)
cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
y         <- foreach::foreach(i = seq(nsim)) %dopar% ( t(chol(covmat[[i]])) %*% rnorm(n = nrow(covmat[[i]])) )
parallel::stopCluster(cl)

## plot
fields::quilt.plot(as.numeric(locs[[1]][, 1]), as.numeric(locs[[1]][, 2]), as.numeric(y[[1]]))
plot(as.numeric(locs[[1]][, 3]), as.numeric(y[[1]]), type = 'o', col = as.factor(locs[[1]][, 1]))

Sys.time()

### baseline_1 #############################################################################################

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

fit_b1 <- list()
for(j in 1:length(m)){
  fit_b1[[j]] <- list()
  fit_b1[[j]] <- foreach::foreach(i = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% fit_spacetime_baseline(y = y[[i]], inputs = locs[[i]], ms = m[j], trend = "intercept", scale = "parms", baseline = 1, nu = covparms[4], covfun = "matern_spacetime", max.it = 100)
}

parallel::stopCluster(cl)

Sys.time()

### baseline_2 #############################################################################################

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

fit_b2 <- list()
for(j in 1:length(m)){
  fit_b2[[j]] <- list()
  fit_b2[[j]] <- foreach::foreach(i = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% fit_spacetime_baseline(y = y[[i]], inputs = locs[[i]], ms = m[j], trend = "intercept", scale = "parms", baseline = 2, nu = covparms[4], covfun = "matern_spacetime", max.it = 100)
}

parallel::stopCluster(cl)

Sys.time()

### baseline_3 #############################################################################################

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

fit_b3 <- list()
for(j in 1:length(m)){
  fit_b3[[j]] <- list()
  fit_b3[[j]] <- foreach::foreach(i = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% fit_spacetime_baseline(y = y[[i]], inputs = locs[[i]], ms = m[j], trend = "intercept", scale = "parms", baseline = 3, nu = covparms[4], covfun = "matern_spacetime", max.it = 100)
}

parallel::stopCluster(cl)

Sys.time()

### corrvecchia #############################################################################################

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

fit_cc <- list()
for(j in 1:length(m)){
  fit_cc[[j]] <- list()
  fit_cc[[j]] <- foreach::foreach(i = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% fit_corrvecchia(y = y[[i]], inputs = locs[[i]], ms = m[j], trend = "intercept", scale = "parms", nu = covparms[4], covfun = "matern_spacetime", max.it = 100)
}

parallel::stopCluster(cl)

Sys.time()

### save #############################################################################################

save(locs, covmat, fit_b1, fit_b2, fit_b3, fit_cc, file = "R/SCSG2021/fisherscoring.RData")

### KLdiv #############################################################################################

kldiv_loop <- function(nsim, m, locs, fit, type) {
  
  result <- list()
  for(j in 1:length(m)) {
    
    print(j)
    
    result[[j]] <- list()
    for(i in 1:nsim) {
      
      if(type == "b1") {
        
        out_specify <- baseline_1_spacetime_specify(locs = locs[[i]], m = m[j])
        
      } else if(type == "b2") {
        
        out_specify <- baseline_2_spacetime_specify(locs = locs[[i]], m = m[j])
        
      } else if(type == "b3") {
        
        out_specify <- baseline_3_spacetime_specify(locs = locs[[i]], m = m[j], covmodel = GpGp::matern_spacetime, covparms = fit[[j]][[i]]$covparms)
        
      } else if(type == "cc") {
        
        covmodel    <- function(locs, covparms) GpGp::matern_spacetime(covparms = covparms, locs = locs)
        out_specify <- corrvecchia_specify_knownCovparms_2(locs = locs[[i]], m = m[j], covmodel = covmodel, covparms = fit[[j]][[i]]$covparms)
        
      } else {
        
        stop("The argument type is invalid!")
        
      }
      
      indmat <- matrix(0, nrow = nrow(locs[[i]]), ncol = nrow(locs[[i]]))
      for(k in 1:nrow(locs[[i]])) indmat[k, out_specify$U.prep$revNNarray[k, ]] <- 1
      indmat <- indmat + t(indmat) ; diag(indmat) <- 1
      # fields::image.plot(indmat)
      
      # ?GpGp::matern_spacetime
      covmat_est  <- GpGp::matern_spacetime(covparms = fit[[j]][[i]]$covparms, locs = locs[[i]])
      mu_est      <- rep(fit[[j]][[i]]$betahat, nrow(locs[[i]]))
      
      revord      <- order(out_specify$ord)
      covmat_ord  <- covmat_est[out_specify$ord, out_specify$ord]
      
      U           <- GPvecchia::createU(out_specify, c(1), 0, covmodel = covmat_ord)$U
      covmat_est  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
      
      # ?correlationVecchia::kldiv
      result[[j]][[i]] <- kldiv(covmat0 = covmat[[i]], covmat1 = covmat_est, mu0 = rep(0, length(mu_est)), mu1 = mu_est)
    }
  }
  
  return(result)
}

meankld_loop <- function(nsim, m, kldiv_obj) {
  
  result <- rep(0, length(m))
  
  for(j in 1:length(m)) {
    for(i in 1:nsim) {
      result[j] <- result[j] + kldiv_obj[[j]][[i]]
    }
  }
  
  return(result / length(m))
}

Sys.time()

kldiv_b1 <- kldiv_loop(nsim = nsim, m = m, locs = locs, fit = fit_b1, type = "b1") ; Sys.time()
kldiv_b2 <- kldiv_loop(nsim = nsim, m = m, locs = locs, fit = fit_b2, type = "b2") ; Sys.time()
kldiv_b3 <- kldiv_loop(nsim = nsim, m = m, locs = locs, fit = fit_b3, type = "b3") ; Sys.time()
kldiv_cc <- kldiv_loop(nsim = nsim, m = m, locs = locs, fit = fit_cc, type = "cc") ; Sys.time()

meankld_b1 <- meankld_loop(nsim = nsim, m = m, kldiv_obj = kldiv_b1)
meankld_b2 <- meankld_loop(nsim = nsim, m = m, kldiv_obj = kldiv_b2)
meankld_b3 <- meankld_loop(nsim = nsim, m = m, kldiv_obj = kldiv_b3)
meankld_cc <- meankld_loop(nsim = nsim, m = m, kldiv_obj = kldiv_cc)

Sys.time()

### save and output #############################################################################################

save(locs, covmat, y, fit_b1, fit_b2, fit_b3, fit_cc, kldiv_b1, kldiv_b2, kldiv_b3, kldiv_cc, meankld_b1, meankld_b2, meankld_b3, meankld_cc, kldiv_loop, meankld_loop, file = "R/SCSG2021/fisherscoring_with_kld.RData")

mvis <- m[1:10] ; viskld_b1 <- meankld_b1[1:10] ; viskld_b2 <- meankld_b2[1:10] ; viskld_b3 <- meankld_b3[1:10] ; viskld_cc <- meankld_cc[1:10]

plot(mvis, log10(viskld_b1), col = 1, type = 'o', ylim = log10(range(c(viskld_b1, viskld_b2, viskld_b3, viskld_cc))), lwd = 2, ylab = "log10(KLdiv)", main = paste0("(variance, range_space, range_time, nu, nugget)", " = ", "(", paste0(covparms, collapse = ", "), ")"))
points(mvis, log10(viskld_b2), col = 3, type = 'o', lwd = 2)
points(mvis, log10(viskld_b3), col = 4, type = 'o', lwd = 2)
points(mvis, log10(viskld_cc), col = 2, type = 'o', lwd = 2)
legend("bottomleft", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c(1, 3, 4, 2), lty = 1, lwd = 2)

vdat1 <- data.frame(index = 1:10, m = mvis, approx_1 = viskld_b1, approx_2 = viskld_b2, approx_3 = viskld_b3, approx_4 = viskld_cc)
vis_est <- vis_SCSG(vdat1 = vdat1, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("R/SCSG2021/spest.pdf", vis_est, width = 9.0, height = 5.7)

Sys.time()

