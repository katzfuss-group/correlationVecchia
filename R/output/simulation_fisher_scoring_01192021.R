####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: .
###
####################################################################################

rm(list = ls())
set.seed(01192021)

library(correlationVecchia) ; library(GpGp)

library(dplyr) ; library(ggplot2) ; library(gridExtra)

library(foreach)

### ncores #############################################################################################

no_cores  <- parallel::detectCores() - 2

### setting #############################################################################################

fit_spacetime_baseline_range <- function(y, inputs, ms = c(30), baseline, sig2 = 1, srange.ini, trange.ini, nu = 0.5, nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
{
  ### dimensions
  n = nrow(inputs)
  d = ncol(inputs)
  
  # X = as.matrix(sample(c(-1, 1), n, replace = TRUE))
  X = as.matrix(rep(1, n))
  
  fix.nug     <- TRUE
  fix.nu      <- TRUE
  
  cur.var     <- sig2
  cur.ranges  <- c(srange.ini, trange.ini) ; cur.srange <- srange.ini ; cur.trange <- trange.ini
  cur.oth     <- c(nu, nug)
  
  ### decrease or remove m values larger than n
  ms = unique(ifelse(ms < length(y), ms, length(y)-1))
  
  ### for increasing m
  for(i.m in 1:length(ms)) {
    
    m = ms[i.m]
    if(i.m < length(ms)){ tol = 10^(-tol.dec-2) } else {tol = 10^(-tol.dec)}
    
    ### increase maxit until convergence
    conv = FALSE
    maxit = 2
    while(conv == FALSE & maxit <= max.it){
      
      if(print.level > 0) {
        print(paste0('m = ', m, ', maxit = ', maxit)); print(cur.ranges)
      }
      
      ## define scale parameters using range parameters
      scales = 1/c(rep(cur.srange, d-1), cur.trange)
      
      ## order and condition based on current params
      ord = order_time(locs = inputs, coordinate = NULL) # timely ordering
      inputs.ord = inputs[ord, , drop = FALSE]
      y.ord = y[ord]
      X.ord = X[ord, , drop = FALSE]
      
      if(baseline == 1) { # T-ord + T-NN
        NNarray = matrix(NA, nrow = n, ncol = m + 1)
        for(i in 1:n) {
          ind               <- seq(from = 1, to = min(i, m + 1), by = 1)
          NNarray[i, ind]   <- seq(from = i, by = -1, length.out = length(ind)) 
        }
      } else if(baseline == 2) { # T-ord + E-NN
        NNarray = GpGp::find_ordered_nn(inputs.ord, m)
      } else if(baseline == 3) { # T-ord + C-NN
        NNarray = GpGp::find_ordered_nn(t(t(inputs.ord)*scales), m)
      } else {
        stop("The argument baseline must be 1, 2, or 3.")
      }
      
      ## starting and fixed parameters
      cur.parms     <- c(sig2, cur.ranges, cur.oth)
      ind.parms     <- c(rep("var", length(sig2)), rep("range", length(c(cur.ranges))), rep("oth", length(cur.oth)))
      
      fixed         <- c(1, 4, 5)
      
      ## fisher scoring
      fit = GpGp::fit_model(y.ord, inputs.ord, X.ord, NNarray = NNarray, m_seq = m, convtol = tol, start_parms = cur.parms, max_iter = maxit, covfun_name = "matern_spacetime", silent = (print.level < 2), reorder = FALSE, fixed_parms = fixed)
      cur.var = fit$covparms[ind.parms == "var"]
      cur.ranges = fit$covparms[ind.parms == "range"] ; cur.srange <- cur.ranges[1] ; cur.trange <- cur.ranges[2]
      cur.oth = fit$covparms[ind.parms == "oth"]
      conv = fit$conv
      maxit = maxit * 2
    }
  }
  
  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var, cur.ranges, cur.oth)
  fit$locs = inputs.ord
  return(fit)
}

fit_spacetime_corrvecchia_range <- function(y, inputs, ms = c(30), sig2 = 1, srange.ini, trange.ini, nu = 0.5, nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
{
  ### dimensions
  n = nrow(inputs)
  d = ncol(inputs)
  
  # X = as.matrix(sample(c(-1, 1), n, replace = TRUE))
  X = as.matrix(rep(1, n))
  
  fix.nug     <- TRUE
  fix.nu      <- TRUE
  
  cur.var     <- sig2
  cur.ranges  <- c(srange.ini, trange.ini) ; cur.srange <- srange.ini ; cur.trange <- trange.ini
  cur.oth     <- c(nu, nug)
  
  ### decrease or remove m values larger than n
  ms = unique(ifelse(ms < length(y), ms, length(y)-1))
  
  ### for increasing m
  for(i.m in 1:length(ms)) {
    
    m = ms[i.m]
    if(i.m < length(ms)){ tol = 10^(-tol.dec-2) } else {tol = 10^(-tol.dec)}
    
    ### increase maxit until convergence
    conv = FALSE
    maxit = 2
    while(conv == FALSE & maxit <= max.it){
      
      if(print.level > 0) {
        print(paste0('m = ', m, ', maxit = ', maxit)); print(cur.ranges)
      }
      
      ## define scale parameters using range parameters
      scales = 1/c(rep(cur.srange, d-1), cur.trange)
      
      ## order and condition based on current params
      ord = GPvecchia::order_maxmin_exact(t(t(inputs)*scales))
      inputs.ord = inputs[ord, , drop = FALSE]
      y.ord = y[ord]
      X.ord = X[ord, , drop = FALSE]
      NNarray = GpGp::find_ordered_nn(t(t(inputs.ord)*scales), m)
      
      ## starting and fixed parameters
      cur.parms     <- c(sig2, cur.ranges, cur.oth)
      ind.parms     <- c(rep("var", length(sig2)), rep("range", length(c(cur.ranges))), rep("oth", length(cur.oth)))
      
      fixed         <- c(1, 4, 5)
      
      ## fisher scoring
      fit = GpGp::fit_model(y.ord, inputs.ord, X.ord, NNarray = NNarray, m_seq = m, convtol = tol, start_parms = cur.parms, max_iter = maxit, covfun_name = "matern_spacetime", silent = (print.level < 2), reorder = FALSE, fixed_parms = fixed)
      cur.var = fit$covparms[ind.parms == "var"]
      cur.ranges = fit$covparms[ind.parms == "range"] ; cur.srange <- cur.ranges[1] ; cur.trange <- cur.ranges[2]
      cur.oth = fit$covparms[ind.parms == "oth"]
      conv = fit$conv
      maxit = maxit * 2
    }
  }
  
  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var, cur.ranges, cur.oth)
  fit$locs = inputs.ord
  return(fit)
}

fit_spacetime_exactGP_range <- function(y, inputs, sig2 = 1, srange.ini, trange.ini, nu = 0.5, nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
{
  ### dimensions
  n = nrow(inputs)
  d = ncol(inputs)
  ms = n - 1
  
  # X = as.matrix(sample(c(-1, 1), n, replace = TRUE))
  X = as.matrix(rep(1, n))
  
  fix.nug     <- TRUE
  fix.nu      <- TRUE
  
  cur.var     <- sig2
  cur.ranges  <- c(srange.ini, trange.ini) ; cur.srange <- srange.ini ; cur.trange <- trange.ini
  cur.oth     <- c(nu, nug)
  
  ### decrease or remove m values larger than n
  ms = unique(ifelse(ms < length(y), ms, length(y)-1))
  
  ### for increasing m
  for(i.m in 1:length(ms)) {
    
    m = ms[i.m]
    if(i.m < length(ms)){ tol = 10^(-tol.dec-2) } else {tol = 10^(-tol.dec)}
    
    ### increase maxit until convergence
    conv = FALSE
    maxit = 2
    while(conv == FALSE & maxit <= max.it){
      
      if(print.level > 0) {
        print(paste0('m = ', m, ', maxit = ', maxit)); print(cur.ranges)
      }
      
      ## define scale parameters using range parameters
      scales = 1/c(rep(cur.srange, d-1), cur.trange)
      
      ## order and condition based on current params
      ord = seq(n)
      inputs.ord = inputs[ord, , drop = FALSE]
      y.ord = y[ord]
      X.ord = X[ord, , drop = FALSE]
      NNarray = GpGp::find_ordered_nn(t(t(inputs.ord)*scales), n)
      
      ## starting and fixed parameters
      cur.parms     <- c(sig2, cur.ranges, cur.oth)
      ind.parms     <- c(rep("var", length(sig2)), rep("range", length(c(cur.ranges))), rep("oth", length(cur.oth)))
      
      fixed         <- c(1, 4, 5)
      
      ## fisher scoring
      fit = GpGp::fit_model(y.ord, inputs.ord, X.ord, NNarray = NNarray, m_seq = m, convtol = tol, start_parms = cur.parms, max_iter = maxit, covfun_name = "matern_spacetime", silent = (print.level < 2), reorder = FALSE, fixed_parms = fixed)
      cur.var = fit$covparms[ind.parms == "var"]
      cur.ranges = fit$covparms[ind.parms == "range"] ; cur.srange <- cur.ranges[1] ; cur.trange <- cur.ranges[2]
      cur.oth = fit$covparms[ind.parms == "oth"]
      conv = fit$conv
      maxit = maxit * 2
    }
  }
  
  ### prepare fit object for subsequent prediction
  fit$covparms=c(cur.var, cur.ranges, cur.oth)
  fit$locs = inputs.ord
  return(fit)
}

### generation #############################################################################################

nsim      <- 200
m         <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

covparms  <- c(1, 0.1, 0.2, 0.5, 0) # var, range_space, range_time, nu, nugget

# out2      <- generate_gp_spacetime(nsim = nsim, n = 25, d = 2, t.len = 36, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out2      <- generate_gp_spacetime(nsim = nsim, n = 30^2, d = 2, t.len = 1, method.locs = "satellite",covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))

locs      <- list()
for(i in 1:nsim) locs[[i]] <- out2$sim[[i]]$locs %>% as.matrix()

covmat    <- list()
for(i in 1:nsim) covmat[[i]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[i]])

y         <- list()
for(i in 1:nsim) y[[i]] <- t(chol(covmat[[i]])) %*% rnorm(n = nrow(covmat[[i]])) 

### experiment #############################################################################################

srange.ini <- 0.15
trange.ini <- 0.15
verbose <- 0

Sys.time()

#####

# fit_b1 <- list()
# for(i in 1:length(m)) {
#   
#   fit_b1[[i]] <- list()
#   for(j in 1:nsim) {
#     fit_b1[[i]][[j]] <- fit_spacetime_baseline_range(y[[j]], locs[[j]], ms = m[i], baseline = 1, sig2 = 1, srange.ini = srange.ini, trange.ini = trange.ini, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)
#   }
#   print(i)
# }

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

fit_b1 <- list()
for(i in 1:length(m)) {
  fit_b1[[i]] <- list()
  fit_b1[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% fit_spacetime_baseline_range(y[[j]], locs[[j]], ms = m[i], baseline = 1, sig2 = 1, srange.ini = srange.ini, trange.ini = trange.ini, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)
  
  print(i)
}

parallel::stopCluster(cl)

Sys.time()

#####

# fit_b2 <- list()
# for(i in 1:length(m)) {
#   
#   fit_b2[[i]] <- list()
#   for(j in 1:nsim) {
#     fit_b2[[i]][[j]] <- fit_spacetime_baseline_range(y[[j]], locs[[j]], ms = m[i], baseline = 2, sig2 = 1, srange.ini = srange.ini, trange.ini = trange.ini, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)
#   }
#   print(i)
# }

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

fit_b2 <- list()
for(i in 1:length(m)) {
  fit_b2[[i]] <- list()
  fit_b2[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% fit_spacetime_baseline_range(y[[j]], locs[[j]], ms = m[i], baseline = 2, sig2 = 1, srange.ini = srange.ini, trange.ini = trange.ini, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)
  
  print(i)
}

parallel::stopCluster(cl)

Sys.time()

#####

# fit_b3 <- list()
# for(i in 1:length(m)) {
#   
#   fit_b3[[i]] <- list()
#   for(j in 1:nsim) {
#     fit_b3[[i]][[j]] <- fit_spacetime_baseline_range(y[[j]], locs[[j]], ms = m[i], baseline = 3, sig2 = 1, srange.ini = srange.ini, trange.ini = trange.ini, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)
#   }
#   print(i)
# }

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

fit_b3 <- list()
for(i in 1:length(m)) {
  fit_b3[[i]] <- list()
  fit_b3[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% fit_spacetime_baseline_range(y[[j]], locs[[j]], ms = m[i], baseline = 3, sig2 = 1, srange.ini = srange.ini, trange.ini = trange.ini, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)
  
  print(i)
}

parallel::stopCluster(cl)

Sys.time()

#####

# fit_cc <- list()
# for(i in 1:length(m)) {
#   
#   fit_cc[[i]] <- list()
#   for(j in 1:nsim) {
#     fit_cc[[i]][[j]] <- fit_spacetime_corrvecchia_range(y[[j]], locs[[j]], ms = m[i], sig2 = 1, srange.ini = 0.2, trange.ini = 0.2, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)
#   }
#   print(i)
# }

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

fit_cc <- list()
for(i in 1:length(m)) {
  fit_cc[[i]] <- list()
  fit_cc[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% fit_spacetime_corrvecchia_range(y[[j]], locs[[j]], ms = m[i], sig2 = 1, srange.ini = 0.2, trange.ini = 0.2, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)
  
  print(i)
}

parallel::stopCluster(cl)

Sys.time()

#####

# fit_gp <- list()
# for(i in 1:length(m)) {
#   
#   fit_gp[[i]] <- list()
#   for(j in 1:nsim) {
#     fit_gp[[i]][[j]] <- fit_spacetime_exactGP_range(y[[j]], locs[[j]], sig2 = 1, srange.ini = 0.2, trange.ini = 0.2, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)
#   }
#   print(i)
# }

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

fit_gp <- list()
fit_gp <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% fit_spacetime_exactGP_range(y[[j]], locs[[j]], sig2 = 1, srange.ini = 0.2, trange.ini = 0.2, nu = 0.5, nug = 0, print.level = verbose, max.it = 100, tol.dec = 4)

fit_temp <- list()
for(i in 1:length(m)) {
  fit_temp[[i]] <- list()
  for(j in 1:nsim) fit_temp[[i]][[j]] <- fit_gp[[j]]
}

parallel::stopCluster(cl)

fit_gp <- fit_temp ; rm(fit_temp)

Sys.time()

### save #############################################################################################

save(nsim, m, covparms, locs, covmat, y, fit_b1, fit_b2, fit_b3, fit_cc, fit_gp, fit_spacetime_baseline_range, fit_spacetime_corrvecchia_range, fit_spacetime_exactGP_range, file = "pilotstudy_fisherscoring.RData")

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

# kldiv_b1 <- kldiv_loop(nsim = nsim, m = m, locs = locs, fit = fit_b1, type = "b1") ; Sys.time()
# kldiv_b2 <- kldiv_loop(nsim = nsim, m = m, locs = locs, fit = fit_b2, type = "b2") ; Sys.time()
# kldiv_b3 <- kldiv_loop(nsim = nsim, m = m, locs = locs, fit = fit_b3, type = "b3") ; Sys.time()
# kldiv_cc <- kldiv_loop(nsim = nsim, m = m, locs = locs, fit = fit_cc, type = "cc") ; Sys.time()

approxs <- c("b1", "b2", "b3", "cc")
fit_list <- list() ; fit_list[[1]] <- fit_b1 ; fit_list[[2]] <- fit_b2 ; fit_list[[3]] <- fit_b3 ; fit_list[[4]] <- fit_cc

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

kld_list <- foreach::foreach(j = seq(length(approxs)), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% kldiv_loop(nsim = nsim, m = m, locs = locs, fit = fit_list[[j]], type = approxs[j])

parallel::stopCluster(cl)

#####

# kldiv_gp <- list()
# for(i in 1:length(m)) {
#   
#   kldiv_gp[[i]] <- list()
#   for(j in 1:nsim) {
#     kldiv_gp[[i]][[j]] <- correlationVecchia::kldiv(covmat0 = covmat[[i]], covmat1 = GpGp::matern_spacetime(covparms = fit_gp[[i]][[j]]$covparms, locs = locs[[j]]), mu0 = rep(0, length(y[[j]])), mu1 = rep(fit_gp[[i]][[j]]$betahat, length(y[[j]])))
#   }
#   print(i)
# }

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

kldiv_gp <- list()
for(i in 1:length(m)) {
  
  kldiv_gp[[i]] <- list()
  kldiv_gp[[i]] <- foreach::foreach(j = seq(nsim), .packages = c("correlationVecchia", "GPvecchia", "GpGp")) %dopar% correlationVecchia::kldiv(covmat0 = covmat[[j]], covmat1 = GpGp::matern_spacetime(covparms = fit_gp[[i]][[j]]$covparms, locs = locs[[j]]), mu0 = rep(0, length(y[[j]])), mu1 = rep(fit_gp[[i]][[j]]$betahat, length(y[[j]])))
  
}

parallel::stopCluster(cl)

#####

kldiv_b1 <- kld_list[[1]] ; kldiv_b2 <- kld_list[[2]] ; kldiv_b3 <- kld_list[[3]] ; kldiv_cc <- kld_list[[4]]

meankld_b1 <- meankld_loop(nsim = nsim, m = m, kldiv_obj = kldiv_b1)
meankld_b2 <- meankld_loop(nsim = nsim, m = m, kldiv_obj = kldiv_b2)
meankld_b3 <- meankld_loop(nsim = nsim, m = m, kldiv_obj = kldiv_b3)
meankld_cc <- meankld_loop(nsim = nsim, m = m, kldiv_obj = kldiv_cc)
meankld_gp <- meankld_loop(nsim = nsim, m = m, kldiv_obj = kldiv_gp)

Sys.time()

### save and output #############################################################################################

save(nsim, m, covparms, kldiv_b1, kldiv_b2, kldiv_b3, kldiv_cc, kldiv_gp, meankld_b1, meankld_b2, meankld_b3, meankld_cc, meankld_gp, kldiv_loop, meankld_loop, file = "pilotstudy_fisherscoring_with_kld.RData")

plot(m, log10(meankld_b1), col = "#984EA3", type = 'o', ylim = log10(range(c(meankld_b1, meankld_b2, meankld_b3, meankld_cc, meankld_gp))), lwd = 2, ylab = "log10(KLdiv)", main = paste0("(variance, range_space, range_time, nu, nugget)", " = ", "(", paste0(covparms, collapse = ", "), ")"))
points(m, log10(meankld_b2), col = "#4DAF4A", type = 'o', lwd = 2)
points(m, log10(meankld_b3), col = "#377EB8", type = 'o', lwd = 2)
points(m, log10(meankld_cc), col = "#E41A1C", type = 'o', lwd = 2)
points(m, log10(meankld_gp), col = "gray70", type = 'o', lwd = 2)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

### visualization - Likelihood #############################################################################################

compute_loglik <- function(y, locs, fit, approx)
{
  covmat_est  <- GpGp::matern_spacetime(covparms = fit$covparms, locs = locs)
  mu_est      <- rep(fit$betahat, nrow(locs))
  
  revord      <- order(approx$ord)
  covmat_ord  <- covmat_est[approx$ord, approx$ord]
  
  U           <- GPvecchia::createU(approx, c(1), 0, covmodel = covmat_ord)$U
  covmat_est  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
  
  return(mvtnorm::dmvnorm(as.numeric(y), mu_est, covmat_est, log = TRUE))
}

lik_b1 <- list()
for(i in 1:length(m)) {
  
  lik_b1[[i]] <- list()
  for(j in 1:nsim) {
    Sigma <- GpGp::matern_spacetime(fit_b1[[i]][[j]]$covparms, locs[[j]])
    approx <- correlationVecchia::baseline_1_spacetime_specify(locs = locs[[j]], m = m[i])
    
    # lik_b1[[i]][[j]] <- GPvecchia::vecchia_likelihood(y[[j]] - fit_b1[[i]][[j]]$betahat, approx, fit_b1[[i]][[j]]$covparms, nuggets = 0, covmodel = Sigma)
    lik_b1[[i]][[j]] <- compute_loglik(y[[j]], locs[[j]], fit_b1[[i]][[j]], approx)
  }
  print(i)
}

# lik_b1 <- list()
# for(i in 1:length(m)) {
# 
#   lik_b1[[i]] <- list()
#   for(j in 1:nsim) {
#     lik_b1[[i]][[j]] <- fit_b1[[i]][[j]]$loglik
#   }
#   print(i)
# }

Sys.time()

lik_b2 <- list()
for(i in 1:length(m)) {
  
  lik_b2[[i]] <- list()
  for(j in 1:nsim) {
    Sigma <- GpGp::matern_spacetime(fit_b2[[i]][[j]]$covparms, locs[[j]])
    approx <- correlationVecchia::baseline_2_spacetime_specify(locs = locs[[j]], m = m[i])
    
    # lik_b2[[i]][[j]] <- GPvecchia::vecchia_likelihood(y[[j]] - fit_b2[[i]][[j]]$betahat, approx, fit_b2[[i]][[j]]$covparms, nuggets = 0, covmodel = Sigma)
    lik_b2[[i]][[j]] <- compute_loglik(y[[j]], locs[[j]], fit_b2[[i]][[j]], approx)
  }
  print(i)
}

# lik_b2 <- list()
# for(i in 1:length(m)) {
# 
#   lik_b2[[i]] <- list()
#   for(j in 1:nsim) {
#     lik_b2[[i]][[j]] <- fit_b2[[i]][[j]]$loglik
#   }
#   print(i)
# }

Sys.time()

lik_b3 <- list()
for(i in 1:length(m)) {
  
  lik_b3[[i]] <- list()
  for(j in 1:nsim) {
    Sigma <- GpGp::matern_spacetime(fit_b3[[i]][[j]]$covparms, locs[[j]])
    approx <- correlationVecchia::baseline_3_spacetime_specify(locs = locs[[j]], m = m[i], covmodel = GpGp::matern_spacetime, covparms = fit_b3[[i]][[j]]$covparms)
    
    # lik_b3[[i]][[j]] <- GPvecchia::vecchia_likelihood(y[[j]] - fit_b3[[i]][[j]]$betahat, approx, fit_b3[[i]][[j]]$covparms, nuggets = 0, covmodel = Sigma)
    lik_b3[[i]][[j]] <- compute_loglik(y[[j]], locs[[j]], fit_b3[[i]][[j]], approx)
  }
  print(i)
}

# lik_b3 <- list()
# for(i in 1:length(m)) {
# 
#   lik_b3[[i]] <- list()
#   for(j in 1:nsim) {
#     lik_b3[[i]][[j]] <- fit_b3[[i]][[j]]$loglik
#   }
#   print(i)
# }

Sys.time()

lik_cc <- list()
for(i in 1:length(m)) {
  
  lik_cc[[i]] <- list()
  for(j in 1:nsim) {
    Sigma <- GpGp::matern_spacetime(fit_cc[[i]][[j]]$covparms, locs[[j]])
    
    approx <- correlationVecchia::corrvecchia_specify_knownCovparms_2(locs = locs[[j]], m = m[i], covmodel = Sigma, covparms = fit_cc[[i]][[j]]$covparms)
    # covmodel <- function(locs, covparms) GpGp::matern_spacetime(covparms = covparms, locs = locs)
    # approx <- correlationVecchia::corrvecchia_specify_knownCovparms_2(locs = locs[[j]], m = m[i], covmodel = covmodel, covparms = fit_cc[[i]][[j]]$covparms)
    
    # lik_cc[[i]][[j]] <- GPvecchia::vecchia_likelihood(y[[j]] - fit_cc[[i]][[j]]$betahat, approx, fit_cc[[i]][[j]]$covparms, nuggets = 0, covmodel = Sigma)
    lik_cc[[i]][[j]] <- compute_loglik(y[[j]], locs[[j]], fit_cc[[i]][[j]], approx)
  }
  print(i)
}

# lik_cc <- list()
# for(i in 1:length(m)) {
# 
#   lik_cc[[i]] <- list()
#   for(j in 1:nsim) {
#     lik_cc[[i]][[j]] <- fit_cc[[i]][[j]]$loglik
#   }
#   print(i)
# }

Sys.time()

lik_gp <- list()
for(i in 1:length(m)) {
  
  lik_gp[[i]] <- list()
  for(j in 1:nsim) {
    Sigma <- GpGp::matern_spacetime(fit_gp[[i]][[j]]$covparms, locs[[j]])
    covmodel <- function(locs, covparms) GpGp::matern_spacetime(covparms = covparms, locs = locs)
    approx <- correlationVecchia::baseline_1_spacetime_specify(locs = locs[[j]], m = length(y[[j]]))
    
    # lik_gp[[i]][[j]] <- GPvecchia::vecchia_likelihood(y[[j]] - fit_gp[[i]][[j]]$betahat, approx, fit_gp[[i]][[j]]$covparms, nuggets = 0, covmodel = Sigma)
    lik_gp[[i]][[j]] <- compute_loglik(y[[j]], locs[[j]], fit_gp[[i]][[j]], approx)
  }
  print(i)
}

# lik_gp <- list()
# for(i in 1:length(m)) {
# 
#   lik_gp[[i]] <- list()
#   for(j in 1:nsim) {
#     lik_gp[[i]][[j]] <- fit_gp[[i]][[j]]$loglik
#   }
#   print(i)
# }

Sys.time()

# vislik_b1 <- rep(0, length(m))
# for(i in 1:length(m)) {
#   vislik_b1[i] <- lik_b1[[i]][[1]] 
# }

meanlik_b1 <- rep(0, length(m))
for(i in 1:length(m)) {
  for(j in 1:nsim) {
    meanlik_b1[i] <- meanlik_b1[i] + lik_b1[[i]][[j]] 
  }
  
  meanlik_b1[i] <- meanlik_b1[i] / nsim
}

# vislik_b2 <- rep(0, length(m))
# for(i in 1:length(m)) {
#   vislik_b2[i] <- lik_b2[[i]][[1]] 
# }

meanlik_b2 <- rep(0, length(m))
for(i in 1:length(m)) {
  for(j in 1:nsim) {
    meanlik_b2[i] <- meanlik_b2[i] + lik_b2[[i]][[j]] 
  }
  
  meanlik_b2[i] <- meanlik_b2[i] / nsim
}

# vislik_b3 <- rep(0, length(m))
# for(i in 1:length(m)) {
#   vislik_b3[i] <- lik_b3[[i]][[1]] 
# }

meanlik_b3 <- rep(0, length(m))
for(i in 1:length(m)) {
  for(j in 1:nsim) {
    meanlik_b3[i] <- meanlik_b3[i] + lik_b3[[i]][[j]] 
  }
  
  meanlik_b3[i] <- meanlik_b3[i] / nsim
}

# vislik_cc <- rep(0, length(m))
# for(i in 1:length(m)) {
#   vislik_cc[i] <- lik_cc[[i]][[1]] 
# }

meanlik_cc <- rep(0, length(m))
for(i in 1:length(m)) {
  for(j in 1:nsim) {
    meanlik_cc[i] <- meanlik_cc[i] + lik_cc[[i]][[j]] 
  }
  
  meanlik_cc[i] <- meanlik_cc[i] / nsim
}

# vislik_gp <- rep(0, length(m))
# for(i in 1:length(m)) {
#   vislik_gp[i] <- lik_gp[[i]][[1]] 
# }

meanlik_gp <- rep(0, length(m))
for(i in 1:length(m)) {
  for(j in 1:nsim) {
    meanlik_gp[i] <- meanlik_gp[i] + lik_gp[[i]][[j]] 
  }
  
  meanlik_gp[i] <- meanlik_gp[i] / nsim
}

### save and output #############################################################################################

save(nsim, m, covparms, lik_b1, lik_b2, lik_b3, lik_cc, lik_gp, meanlik_b1, meanlik_b2, meanlik_b3, meanlik_cc, meanlik_gp, compute_loglik, file = "pilotstudy_fisherscoring_with_lik.RData")

plot(m, -meanlik_b1, col = "#984EA3", type = 'o', ylim = range(c(-meanlik_b1, -meanlik_b2, -meanlik_b3, -meanlik_cc, -meanlik_gp)), lwd = 2, ylab = "negative log-liklihood", main = paste0("(variance, range_space, range_time, nu, nugget)", " = ", "(", paste0(covparms, collapse = ", "), ")"), xlab = "m")
points(m, -meanlik_b2, col = "#4DAF4A", type = 'o', lwd = 2)
points(m, -meanlik_b3, col = "#377EB8", type = 'o', lwd = 2)
points(m, -meanlik_cc, col = "#E41A1C", type = 'o', lwd = 2)
points(m, -meanlik_gp, col = "gray70", type = 'o', lwd = 2)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

### visualization - MSE #############################################################################################

mse_srange_b1 <- c()
mse_trange_b1 <- c()
for(i in 1:length(m)) {
  
  mse_srange_b1[i] <- 0
  mse_trange_b1[i] <- 0
  for(j in 1:nsim) {
    mse_srange_b1[i] <- mse_srange_b1[i] + (fit_b1[[i]][[j]]$covparms[2] - covparms[2])^2
    mse_trange_b1[i] <- mse_trange_b1[i] + (fit_b1[[i]][[j]]$covparms[3] - covparms[3])^2
  }
  
  mse_srange_b1[i] <- mse_srange_b1[i] / nsim
  mse_trange_b1[i] <- mse_trange_b1[i] / nsim
}

mse_srange_b2 <- c()
mse_trange_b2 <- c()
for(i in 1:length(m)) {
  
  mse_srange_b2[i] <- 0
  mse_trange_b2[i] <- 0
  for(j in 1:nsim) {
    mse_srange_b2[i] <- mse_srange_b2[i] + (fit_b2[[i]][[j]]$covparms[2] - covparms[2])^2
    mse_trange_b2[i] <- mse_trange_b2[i] + (fit_b2[[i]][[j]]$covparms[3] - covparms[3])^2
  }
  
  mse_srange_b2[i] <- mse_srange_b2[i] / nsim
  mse_trange_b2[i] <- mse_trange_b2[i] / nsim
}

mse_srange_b3 <- c()
mse_trange_b3 <- c()
for(i in 1:length(m)) {
  
  mse_srange_b3[i] <- 0
  mse_trange_b3[i] <- 0
  for(j in 1:nsim) {
    mse_srange_b3[i] <- mse_srange_b3[i] + (fit_b3[[i]][[j]]$covparms[2] - covparms[2])^2
    mse_trange_b3[i] <- mse_trange_b3[i] + (fit_b3[[i]][[j]]$covparms[3] - covparms[3])^2
  }
  
  mse_srange_b3[i] <- mse_srange_b3[i] / nsim
  mse_trange_b3[i] <- mse_trange_b3[i] / nsim
}

mse_srange_cc <- c()
mse_trange_cc <- c()
for(i in 1:length(m)) {
  
  mse_srange_cc[i] <- 0
  mse_trange_cc[i] <- 0
  for(j in 1:nsim) {
    mse_srange_cc[i] <- mse_srange_cc[i] + (fit_cc[[i]][[j]]$covparms[2] - covparms[2])^2
    mse_trange_cc[i] <- mse_trange_cc[i] + (fit_cc[[i]][[j]]$covparms[3] - covparms[3])^2
  }
  
  mse_srange_cc[i] <- mse_srange_cc[i] / nsim
  mse_trange_cc[i] <- mse_trange_cc[i] / nsim
}

mse_srange_gp <- c()
mse_trange_gp <- c()
for(i in 1:length(m)) {
  
  mse_srange_gp[i] <- 0
  mse_trange_gp[i] <- 0
  for(j in 1:nsim) {
    mse_srange_gp[i] <- mse_srange_gp[i] + (fit_gp[[i]][[j]]$covparms[2] - covparms[2])^2
    mse_trange_gp[i] <- mse_trange_gp[i] + (fit_gp[[i]][[j]]$covparms[3] - covparms[3])^2
  }
  
  mse_srange_gp[i] <- mse_srange_gp[i] / nsim
  mse_trange_gp[i] <- mse_trange_gp[i] / nsim
}

mse_srange <- list()
mse_srange[[1]] <- mse_srange_b1 ; rm(mse_srange_b1)
mse_srange[[2]] <- mse_srange_b2 ; rm(mse_srange_b2)
mse_srange[[3]] <- mse_srange_b3 ; rm(mse_srange_b3)
mse_srange[[4]] <- mse_srange_cc ; rm(mse_srange_cc)
mse_srange[[5]] <- mse_srange_gp ; rm(mse_srange_gp)

mse_trange <- list()
mse_trange[[1]] <- mse_trange_b1 ; rm(mse_trange_b1)
mse_trange[[2]] <- mse_trange_b2 ; rm(mse_trange_b2)
mse_trange[[3]] <- mse_trange_b3 ; rm(mse_trange_b3)
mse_trange[[4]] <- mse_trange_cc ; rm(mse_trange_cc)
mse_trange[[5]] <- mse_trange_gp ; rm(mse_trange_gp)

### save and output #############################################################################################

save(nsim, m, covparms, mse_srange, mse_trange, file = "pilotstudy_fisherscoring_with_mse.RData")

par(mfrow = c(1, 3))
plot(m, mse_srange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "MSE for spatial range parameter estimation", xlab = "m", ylab = "MSE", ylim = range(c(mse_srange[[1]], mse_srange[[2]], mse_srange[[3]], mse_srange[[4]], mse_srange[[5]])))
lines(m, mse_srange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, mse_srange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mse_srange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mse_srange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, mse_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "MSE for temporal range parameter estimation", xlab = "m", ylab = "MSE", ylim = range(c(mse_trange[[1]], mse_trange[[2]], mse_trange[[3]], mse_trange[[4]], mse_trange[[5]])))
lines(m, mse_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, mse_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mse_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mse_trange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, mse_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "MSE for temporal range parameter estimation", xlab = "m", ylab = "MSE", ylim = range(c(mse_trange[[2]], mse_trange[[3]], mse_trange[[4]], mse_trange[[5]])))
lines(m, mse_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, mse_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mse_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mse_trange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)
par(mfrow = c(1, 1))

### visualization - MAE #############################################################################################

mae_srange_b1 <- c()
mae_trange_b1 <- c()
for(i in 1:length(m)) {
  
  mae_srange_b1[i] <- 0
  mae_trange_b1[i] <- 0
  for(j in 1:nsim) {
    mae_srange_b1[i] <- mae_srange_b1[i] + abs(fit_b1[[i]][[j]]$covparms[2] - covparms[2])
    mae_trange_b1[i] <- mae_trange_b1[i] + abs(fit_b1[[i]][[j]]$covparms[3] - covparms[3])
  }
  
  mae_srange_b1[i] <- mae_srange_b1[i] / nsim
  mae_trange_b1[i] <- mae_trange_b1[i] / nsim
}

mae_srange_b2 <- c()
mae_trange_b2 <- c()
for(i in 1:length(m)) {
  
  mae_srange_b2[i] <- 0
  mae_trange_b2[i] <- 0
  for(j in 1:nsim) {
    mae_srange_b2[i] <- mae_srange_b2[i] + abs(fit_b2[[i]][[j]]$covparms[2] - covparms[2])
    mae_trange_b2[i] <- mae_trange_b2[i] + abs(fit_b2[[i]][[j]]$covparms[3] - covparms[3])
  }
  
  mae_srange_b2[i] <- mae_srange_b2[i] / nsim
  mae_trange_b2[i] <- mae_trange_b2[i] / nsim
}

mae_srange_b3 <- c()
mae_trange_b3 <- c()
for(i in 1:length(m)) {
  
  mae_srange_b3[i] <- 0
  mae_trange_b3[i] <- 0
  for(j in 1:nsim) {
    mae_srange_b3[i] <- mae_srange_b3[i] + abs(fit_b3[[i]][[j]]$covparms[2] - covparms[2])
    mae_trange_b3[i] <- mae_trange_b3[i] + abs(fit_b3[[i]][[j]]$covparms[3] - covparms[3])
  }
  
  mae_srange_b3[i] <- mae_srange_b3[i] / nsim
  mae_trange_b3[i] <- mae_trange_b3[i] / nsim
}

mae_srange_cc <- c()
mae_trange_cc <- c()
for(i in 1:length(m)) {
  
  mae_srange_cc[i] <- 0
  mae_trange_cc[i] <- 0
  for(j in 1:nsim) {
    mae_srange_cc[i] <- mae_srange_cc[i] + abs(fit_cc[[i]][[j]]$covparms[2] - covparms[2])
    mae_trange_cc[i] <- mae_trange_cc[i] + abs(fit_cc[[i]][[j]]$covparms[3] - covparms[3])
  }
  
  mae_srange_cc[i] <- mae_srange_cc[i] / nsim
  mae_trange_cc[i] <- mae_trange_cc[i] / nsim
}

mae_srange_gp <- c()
mae_trange_gp <- c()
for(i in 1:length(m)) {
  
  mae_srange_gp[i] <- 0
  mae_trange_gp[i] <- 0
  for(j in 1:nsim) {
    mae_srange_gp[i] <- mae_srange_gp[i] + abs(fit_gp[[i]][[j]]$covparms[2] - covparms[2])
    mae_trange_gp[i] <- mae_trange_gp[i] + abs(fit_gp[[i]][[j]]$covparms[3] - covparms[3])
  }
  
  mae_srange_gp[i] <- mae_srange_gp[i] / nsim
  mae_trange_gp[i] <- mae_trange_gp[i] / nsim
}

mae_srange <- list()
mae_srange[[1]] <- mae_srange_b1 ; rm(mae_srange_b1)
mae_srange[[2]] <- mae_srange_b2 ; rm(mae_srange_b2)
mae_srange[[3]] <- mae_srange_b3 ; rm(mae_srange_b3)
mae_srange[[4]] <- mae_srange_cc ; rm(mae_srange_cc)
mae_srange[[5]] <- mae_srange_gp ; rm(mae_srange_gp)

mae_trange <- list()
mae_trange[[1]] <- mae_trange_b1 ; rm(mae_trange_b1)
mae_trange[[2]] <- mae_trange_b2 ; rm(mae_trange_b2)
mae_trange[[3]] <- mae_trange_b3 ; rm(mae_trange_b3)
mae_trange[[4]] <- mae_trange_cc ; rm(mae_trange_cc)
mae_trange[[5]] <- mae_trange_gp ; rm(mae_trange_gp)

### save and output #############################################################################################

save(nsim, m, covparms, mae_srange, mae_trange, file = "pilotstudy_fisherscoring_with_mae.RData")

par(mfrow = c(1, 3))
plot(m, mae_srange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "MAE for spatial range parameter estimation", xlab = "m", ylab = "mae", ylim = range(c(mae_srange[[1]], mae_srange[[2]], mae_srange[[3]], mae_srange[[4]], mae_srange[[5]])))
lines(m, mae_srange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, mae_srange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mae_srange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mae_srange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, mae_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "MAE for temporal range parameter estimation", xlab = "m", ylab = "mae", ylim = range(c(mae_trange[[1]], mae_trange[[2]], mae_trange[[3]], mae_trange[[4]], mae_trange[[5]])))
lines(m, mae_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, mae_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mae_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mae_trange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, mae_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "MAE for temporal range parameter estimation", xlab = "m", ylab = "mae", ylim = range(c(mae_trange[[2]], mae_trange[[3]], mae_trange[[4]], mae_trange[[5]])))
lines(m, mae_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, mae_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mae_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mae_trange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)
par(mfrow = c(1, 1))

### visualization - MAD, MSD #############################################################################################

vislist_srange <- list()
vislist_trange <- list()
for(j in 1:nsim) {
  vislist_srange[[j]] <- data.frame(nsim = j, m = m, est_b1 = NA, est_b2 = NA, est_b3 = NA, est_cc = NA, est_gp = NA)
  vislist_trange[[j]] <- data.frame(nsim = j, m = m, est_b1 = NA, est_b2 = NA, est_b3 = NA, est_cc = NA, est_gp = NA)
  for(i in 1:length(m)) {
    vislist_srange[[j]]$est_b1[i] <- fit_b1[[i]][[j]]$covparms[2]
    vislist_srange[[j]]$est_b2[i] <- fit_b2[[i]][[j]]$covparms[2]
    vislist_srange[[j]]$est_b3[i] <- fit_b3[[i]][[j]]$covparms[2]
    vislist_srange[[j]]$est_cc[i] <- fit_cc[[i]][[j]]$covparms[2]
    vislist_srange[[j]]$est_gp[i] <- fit_gp[[i]][[j]]$covparms[2]
    
    vislist_trange[[j]]$est_b1[i] <- fit_b1[[i]][[j]]$covparms[3]
    vislist_trange[[j]]$est_b2[i] <- fit_b2[[i]][[j]]$covparms[3]
    vislist_trange[[j]]$est_b3[i] <- fit_b3[[i]][[j]]$covparms[3]
    vislist_trange[[j]]$est_cc[i] <- fit_cc[[i]][[j]]$covparms[3]
    vislist_trange[[j]]$est_gp[i] <- fit_gp[[i]][[j]]$covparms[3]
  }
}

visdat_srange <- vislist_srange[[1]]
visdat_trange <- vislist_trange[[1]]
for(j in 2:nsim) {
  visdat_srange <- rbind(visdat_srange, vislist_srange[[j]])
  visdat_trange <- rbind(visdat_trange, vislist_trange[[j]])
}

library(dplyr)

visdat_srange <- visdat_srange %>% mutate(ad_b1 = abs(est_b1 - est_gp), ad_b2 = abs(est_b2 - est_gp), ad_b3 = abs(est_b3 - est_gp), ad_cc = abs(est_cc - est_gp))
visdat_trange <- visdat_trange %>% mutate(ad_b1 = abs(est_b1 - est_gp), ad_b2 = abs(est_b2 - est_gp), ad_b3 = abs(est_b3 - est_gp), ad_cc = abs(est_cc - est_gp))

visdat_srange <- visdat_srange %>% mutate(sd_b1 = (est_b1 - est_gp)^2, sd_b2 = (est_b2 - est_gp)^2, sd_b3 = (est_b3 - est_gp)^2, sd_cc = (est_cc - est_gp)^2)
visdat_trange <- visdat_trange %>% mutate(sd_b1 = (est_b1 - est_gp)^2, sd_b2 = (est_b2 - est_gp)^2, sd_b3 = (est_b3 - est_gp)^2, sd_cc = (est_cc - est_gp)^2)

mad_srange <- list()

mad_srange[[1]] <- rep(NA, length(m))
for(i in 1:length(m)) mad_srange[[1]][i] <- visdat_srange %>% filter(m == m[i]) %>% pull(ad_b1) %>% mean()

mad_srange[[2]] <- rep(NA, length(m))
for(i in 1:length(m)) mad_srange[[2]][i] <- visdat_srange %>% filter(m == m[i]) %>% pull(ad_b2) %>% mean()

mad_srange[[3]] <- rep(NA, length(m))
for(i in 1:length(m)) mad_srange[[3]][i] <- visdat_srange %>% filter(m == m[i]) %>% pull(ad_b3) %>% mean()

mad_srange[[4]] <- rep(NA, length(m))
for(i in 1:length(m)) mad_srange[[4]][i] <- visdat_srange %>% filter(m == m[i]) %>% pull(ad_cc) %>% mean()

msd_srange <- list()

msd_srange[[1]] <- rep(NA, length(m))
for(i in 1:length(m)) msd_srange[[1]][i] <- visdat_srange %>% filter(m == m[i]) %>% pull(sd_b1) %>% mean()

msd_srange[[2]] <- rep(NA, length(m))
for(i in 1:length(m)) msd_srange[[2]][i] <- visdat_srange %>% filter(m == m[i]) %>% pull(sd_b2) %>% mean()

msd_srange[[3]] <- rep(NA, length(m))
for(i in 1:length(m)) msd_srange[[3]][i] <- visdat_srange %>% filter(m == m[i]) %>% pull(sd_b3) %>% mean()

msd_srange[[4]] <- rep(NA, length(m))
for(i in 1:length(m)) msd_srange[[4]][i] <- visdat_srange %>% filter(m == m[i]) %>% pull(sd_cc) %>% mean()

mad_trange <- list()

mad_trange[[1]] <- rep(NA, length(m))
for(i in 1:length(m)) mad_trange[[1]][i] <- visdat_trange %>% filter(m == m[i]) %>% pull(ad_b1) %>% mean()

mad_trange[[2]] <- rep(NA, length(m))
for(i in 1:length(m)) mad_trange[[2]][i] <- visdat_trange %>% filter(m == m[i]) %>% pull(ad_b2) %>% mean()

mad_trange[[3]] <- rep(NA, length(m))
for(i in 1:length(m)) mad_trange[[3]][i] <- visdat_trange %>% filter(m == m[i]) %>% pull(ad_b3) %>% mean()

mad_trange[[4]] <- rep(NA, length(m))
for(i in 1:length(m)) mad_trange[[4]][i] <- visdat_trange %>% filter(m == m[i]) %>% pull(ad_cc) %>% mean()

msd_trange <- list()

msd_trange[[1]] <- rep(NA, length(m))
for(i in 1:length(m)) msd_trange[[1]][i] <- visdat_trange %>% filter(m == m[i]) %>% pull(sd_b1) %>% mean()

msd_trange[[2]] <- rep(NA, length(m))
for(i in 1:length(m)) msd_trange[[2]][i] <- visdat_trange %>% filter(m == m[i]) %>% pull(sd_b2) %>% mean()

msd_trange[[3]] <- rep(NA, length(m))
for(i in 1:length(m)) msd_trange[[3]][i] <- visdat_trange %>% filter(m == m[i]) %>% pull(sd_b3) %>% mean()

msd_trange[[4]] <- rep(NA, length(m))
for(i in 1:length(m)) msd_trange[[4]][i] <- visdat_trange %>% filter(m == m[i]) %>% pull(sd_cc) %>% mean()

### save and output #############################################################################################

save(nsim, m, covparms, mad_srange, mad_trange, file = "pilotstudy_fisherscoring_with_mad.RData")

save(nsim, m, covparms, msd_srange, msd_trange, file = "pilotstudy_fisherscoring_with_msd.RData")

par(mfrow = c(1, 3))
plot(m, mad_srange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "mad for spatial range parameter estimation", xlab = "m", ylab = "mad", ylim = range(c(mad_srange[[1]], mad_srange[[2]], mad_srange[[3]], mad_srange[[4]])))
lines(m, mad_srange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, mad_srange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mad_srange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

plot(m, mad_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "mad for temporal range parameter estimation", xlab = "m", ylab = "mad", ylim = range(c(mad_trange[[1]], mad_trange[[2]], mad_trange[[3]], mad_trange[[4]])))
lines(m, mad_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, mad_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mad_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

plot(m, mad_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "mad for temporal range parameter estimation", xlab = "m", ylab = "mad", ylim = range(c(mad_trange[[2]], mad_trange[[3]], mad_trange[[4]])))
lines(m, mad_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, mad_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mad_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)
par(mfrow = c(1, 1))

par(mfrow = c(1, 3))
plot(m, msd_srange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "msd for spatial range parameter estimation", xlab = "m", ylab = "msd", ylim = range(c(msd_srange[[1]], msd_srange[[2]], msd_srange[[3]], msd_srange[[4]])))
lines(m, msd_srange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, msd_srange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, msd_srange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

plot(m, msd_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "msd for temporal range parameter estimation", xlab = "m", ylab = "msd", ylim = range(c(msd_trange[[1]], msd_trange[[2]], msd_trange[[3]], msd_trange[[4]])))
lines(m, msd_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, msd_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, msd_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

plot(m, msd_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "msd for temporal range parameter estimation", xlab = "m", ylab = "msd", ylim = range(c(msd_trange[[2]], msd_trange[[3]], msd_trange[[4]])))
lines(m, msd_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, msd_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, msd_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)
par(mfrow = c(1, 1))

### comparison #############################################################################################

par(mfrow = c(2, 3))

plot(m, mae_srange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "MAE for spatial range parameter estimation", xlab = "m", ylab = "mae", ylim = range(c(mae_srange[[1]], mae_srange[[2]], mae_srange[[3]], mae_srange[[4]], mae_srange[[5]])))
lines(m, mae_srange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, mae_srange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mae_srange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mae_srange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, mse_srange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "MSE for spatial range parameter estimation", xlab = "m", ylab = "MSE", ylim = range(c(mse_srange[[1]], mse_srange[[2]], mse_srange[[3]], mse_srange[[4]], mse_srange[[5]])))
lines(m, mse_srange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, mse_srange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mse_srange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mse_srange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, sqrt(mse_srange[[1]]), col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "RMSE for spatial range parameter estimation", xlab = "m", ylab = "RMSE", ylim = sqrt(range(c(mse_srange[[1]], mse_srange[[2]], mse_srange[[3]], mse_srange[[4]], mse_srange[[5]]))))
lines(m, sqrt(mse_srange[[2]]), col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(mse_srange[[3]]), col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(mse_srange[[4]]), col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(mse_srange[[5]]), col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, mad_srange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "MAD for spatial range parameter estimation", xlab = "m", ylab = "mad", ylim = range(c(mad_srange[[1]], mad_srange[[2]], mad_srange[[3]], mad_srange[[4]])))
lines(m, mad_srange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, mad_srange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mad_srange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

plot(m, msd_srange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "MSD for spatial range parameter estimation", xlab = "m", ylab = "MSD", ylim = range(c(msd_srange[[1]], msd_srange[[2]], msd_srange[[3]], msd_srange[[4]])))
lines(m, msd_srange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, msd_srange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, msd_srange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

plot(m, sqrt(msd_srange[[1]]), col = "#984EA3", type = 'o', lwd = 2, pch = 16, main = "RMSD for spatial range parameter estimation", xlab = "m", ylab = "RMSD", ylim = sqrt(range(c(msd_srange[[1]], msd_srange[[2]], msd_srange[[3]], msd_srange[[4]]))))
lines(m, sqrt(msd_srange[[2]]), col = "#4DAF4A", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(msd_srange[[3]]), col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(msd_srange[[4]]), col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

par(mfrow = c(1, 1))


par(mfrow = c(2, 3))

plot(m, mae_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "MAE for temporal range parameter estimation", xlab = "m", ylab = "mae", ylim = range(c(mae_trange[[2]], mae_trange[[3]], mae_trange[[4]], mae_trange[[5]])))
lines(m, mae_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, mae_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mae_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mae_trange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, mse_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "MSE for temporal range parameter estimation", xlab = "m", ylab = "MSE", ylim = range(c(mse_trange[[2]], mse_trange[[3]], mse_trange[[4]], mse_trange[[5]])))
lines(m, mse_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, mse_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mse_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, mse_trange[[5]], col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, sqrt(mse_trange[[2]]), col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "RMSE for temporal range parameter estimation", xlab = "m", ylab = "RMSE", ylim = sqrt(range(c(mse_trange[[2]], mse_trange[[3]], mse_trange[[4]], mse_trange[[5]]))))
lines(m, sqrt(mse_trange[[1]]), col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(mse_trange[[3]]), col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(mse_trange[[4]]), col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(mse_trange[[5]]), col = "gray70", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "exact GP"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lty = 1, lwd = 2)

plot(m, mad_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "MAD for temporal range parameter estimation", xlab = "m", ylab = "MAD", ylim = range(c(mad_trange[[2]], mad_trange[[3]], mad_trange[[4]])))
lines(m, mad_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, mad_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, mad_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

plot(m, msd_trange[[2]], col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "MSD for temporal range parameter estimation", xlab = "m", ylab = "MSD", ylim = range(c(msd_trange[[2]], msd_trange[[3]], msd_trange[[4]])))
lines(m, msd_trange[[1]], col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, msd_trange[[3]], col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, msd_trange[[4]], col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

plot(m, sqrt(msd_trange[[2]]), col = "#4DAF4A", type = 'o', lwd = 2, pch = 16, main = "RMSD for temporal range parameter estimation", xlab = "m", ylab = "RMSD", ylim = sqrt(range(c(msd_trange[[2]], msd_trange[[3]], msd_trange[[4]]))))
lines(m, sqrt(msd_trange[[1]]), col = "#984EA3", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(msd_trange[[3]]), col = "#377EB8", type = 'o', lwd = 2, pch = 16)
lines(m, sqrt(msd_trange[[4]]), col = "#E41A1C", type = 'o', lwd = 2, pch = 16)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lty = 1, lwd = 2)

par(mfrow = c(1, 1))

### figures - KLD, MSE, MSD #############################################################################################

vis_SCSG <- function(vdat1, legend, color, shape, alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
  vars1     <- vdat1 %>% select(-index, -starts_with("approx")) %>% colnames()
  if(length(vars1) > 1) stop("Please check the argument vdat1.")
  
  xlabel1   <- vdat1 %>% pull(vars1) %>% unique() %>% sort()
  vis1      <- vdat1 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -vars1)
  
  plot1 <- vis1 %>% ggplot(aes(x = get(vars1), y = log10(kldiv), col = approx, shape = approx)) + 
    geom_point(size = size.point) + 
    geom_line(size = size.line, alpha = alpha) + 
    scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
    scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 2, byrow = TRUE)) + 
    scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 2, byrow = TRUE)) +
    xlab(vars1) + ylab('log10(KL)') + 
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text), 
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text), 
          legend.title = element_blank(), 
          legend.position="top",
          legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
          legend.direction = 'horizontal', 
          legend.spacing.x = unit(15, 'pt'),
          plot.margin = unit(size.margin, "pt"))
  
  return(plot1)
}

vis_mxx <- function(vdat1, vdat2 = NULL, ylab_vdat1 = "MSE for spatial range", ylab_vdat2 = "MSE for temporal range", ylim_vdat2, legend, color, shape, alpha = 0.7, size.point = 4, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))
{
  # library(ggplot2) ; library(dplyr) ; library(RColorBrewer) ; library(gridExtra)
  
  # n.space   <- 10
  # legend    <- paste0( legend, paste0(rep(" ", n.space), collapse = "") )
  
  if(is.null(vdat2)) {
    
    stop("not yet!")
    
  } else {
    
    ### Visualize the first data table ###
    
    # Gathering vdat1
    vars1     <- vdat1 %>% select(-index, -starts_with("approx")) %>% colnames()
    if(length(vars1) > 1) stop("Please check the argument vdat1.")
    
    xlabel1   <- vdat1 %>% pull(vars1) %>% unique() %>% sort()
    vis1      <- vdat1 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -vars1)
    
    # plot vis1
    plot1     <- vis1 %>% ggplot(aes(x = get(vars1), y = kldiv, col = approx, shape = approx)) + 
      geom_point(size = size.point) + 
      geom_line(size = size.line, alpha = alpha) + 
      scale_x_continuous(name = vars1, limits = range(xlabel1), breaks = xlabel1) +
      scale_color_manual(values = color, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) + 
      scale_shape_manual(values = shape, labels = legend, guide = guide_legend(nrow = 1, byrow = TRUE)) +
      xlab(vars1) + ylab(ylab_vdat1) + 
      theme(axis.title.x = element_text(size = size.lab), 
            axis.text.x = element_text(size = size.text), 
            axis.title.y = element_text(size = size.lab), 
            axis.text.y = element_text(size = size.text), 
            legend.title = element_blank(), 
            legend.text = element_text(size = size.legend, margin = margin(r = 25, unit = 'pt')), 
            legend.direction = 'horizontal', 
            legend.spacing.x = unit(15, 'pt'),
            plot.margin = unit(size.margin, "pt")) # t, r, b, l
    
    ### Visualize the second data table ###
    
    # Gathering vdat2
    vars2     <- vdat2 %>% select(-index, -starts_with("approx")) %>% colnames()
    if(length(vars2) > 1) stop("Please check the argument vdat2.")
    
    xlabel2   <- vdat2 %>% pull(vars2) %>% unique() %>% sort()
    vis2      <- vdat2 %>% tidyr::gather(key = "approx", value = "kldiv", -index, -vars2)
    
    # plot vis2
    plot2     <- vis2 %>% ggplot(aes(x = get(vars2), y = kldiv, col = approx, shape = approx)) + 
      geom_point(size = size.point) + 
      geom_line(size = size.line, alpha = alpha) + 
      scale_x_continuous(name = vars2, limits = range(xlabel2), breaks = xlabel2) +
      scale_color_manual(values = color, labels = legend) + 
      scale_shape_manual(values = shape, labels = legend) +
      xlab(vars2) + ylab(ylab_vdat2) + coord_cartesian(ylim = ylim_vdat2) +
      theme(axis.title.x = element_text(size = size.lab), 
            axis.text.x = element_text(size = size.text), 
            axis.title.y = element_text(size = size.lab), 
            axis.text.y = element_text(size = size.text), 
            legend.title = element_blank(), 
            legend.text = element_text(size = size.legend), 
            legend.direction = 'horizontal', 
            plot.margin = unit(size.margin, "pt")) # t, r, b, l
    
    ### Creating the legend ###
    
    tmp       <- ggplot_gtable(ggplot_build(plot1))
    leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    mylegend  <- tmp$grobs[[leg]]
    
    ### Merge the two plots ###
    
    result    <- grid.arrange(mylegend, arrangeGrob(plot1 + theme(legend.position="none"), plot2 + theme(legend.position="none"), nrow=1), nrow=2, heights=c(1, 10))
    
  }
  
  return(result)
}

vdat1 <- data.frame(index = 1:10, m = m, approx_1 = meankld_b1, approx_2 = meankld_b2, approx_3 = meankld_b3, approx_4 = meankld_cc, approx_5 = meankld_gp)
vis_est <- vis_SCSG(vdat1 = vdat1, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), shape = c(18, 15, 17, 16, NA))

ggplot2::ggsave("spest_kldiv.pdf", vis_est, width = 9.0, height = 5.7)

vdat1     <- data.frame(index = seq(length(m)), m = m, approx_1 = mse_srange[[1]], approx_2 = mse_srange[[2]], approx_3 = mse_srange[[3]], approx_4 = mse_srange[[4]], approx_5 = mse_srange[[5]])
vdat2     <- data.frame(index = seq(length(m)), m = m, approx_1 = mse_trange[[1]], approx_2 = mse_trange[[2]], approx_3 = mse_trange[[3]], approx_4 = mse_trange[[4]], approx_5 = mse_trange[[5]])
vis_est   <- vis_mxx(vdat1 = vdat1, vdat2 = vdat2, ylab_vdat1 = "MSE for spatial range", ylab_vdat2 = "MSE for temporal range", ylim_vdat2 = range(c(mse_trange[[1]], mse_trange[[2]], mse_trange[[3]], mse_trange[[4]], mse_trange[[5]])), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP"), color =  c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), shape = c(18, 15, 17, 16, NA))
# vis_est   <- vis_mxx(vdat1 = vdat1, vdat2 = vdat2, ylab_vdat1 = "MSE for spatial range", ylab_vdat2 = "MSE for temporal range", ylim_vdat2 = range(c(mse_trange[[2]], mse_trange[[3]], mse_trange[[4]], mse_trange[[5]])), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP"), color =  c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), shape = c(18, 15, 17, 16, NA))

ggplot2::ggsave("spest_mse.pdf", vis_est, width = 15.2, height = 5.7)

vdat1     <- data.frame(index = seq(length(m)), m = m, approx_1 = msd_srange[[1]], approx_2 = msd_srange[[2]], approx_3 = msd_srange[[3]], approx_4 = msd_srange[[4]])
vdat2     <- data.frame(index = seq(length(m)), m = m, approx_1 = msd_trange[[1]], approx_2 = msd_trange[[2]], approx_3 = msd_trange[[3]], approx_4 = msd_trange[[4]])
vis_est   <- vis_mxx(vdat1 = vdat1, vdat2 = vdat2, ylab_vdat1 = "MSD for spatial range", ylab_vdat2 = "MSD for temporal range", ylim_vdat2 = range(c(msd_trange[[1]], msd_trange[[2]], msd_trange[[3]], msd_trange[[4]])), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color =  c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
# vis_est   <- vis_mxx(vdat1 = vdat1, vdat2 = vdat2, ylab_vdat1 = "MSD for spatial range", ylab_vdat2 = "MSD for temporal range", ylim_vdat2 = range(c(msd_trange[[2]], msd_trange[[3]], msd_trange[[4]])), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color =  c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))

ggplot2::ggsave("spest_msd.pdf", vis_est, width = 15.2, height = 5.7)


