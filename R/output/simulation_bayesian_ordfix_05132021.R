####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: 
###
####################################################################################

rm(list = ls())
# set.seed(05042021)
set.seed(05072021)

### load packages
library(correlationVecchia)
library(foreach)

library(RandomFields)
library(mvtnorm)
library(LaplacesDemon)

##### ncores #######################################################################

no_cores  <- parallel::detectCores() - 2

##### cov_matern_spacetime #########################################################

ms              <- c(5, 10, 20, 30, 40)
approxs         <- c("b1", "b2", "b3", "cc", "ex", "tr")
candid          <- expand.grid(ms, approxs) ; colnames(candid) <- c("m", "approx")

n               <- 30^2
covparms        <- c(1, 0.1, 1.0, 0.5)
nugget          <- 0.00

### fisher scoring ####################################################

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

### bayesian ##########################################################

.posterior_fisherord_for_spatialRange <- function(z, locs, covparms, nugget, m, approx, N, xlim, sdlog)
{
  ### prior information
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[2])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  ### ordering / conditioning with fisher-scoring approach
  if(approx == "b1") {
    
    out.fisher      <- fit_spacetime_baseline_range(z, locs, ms = c(m), baseline = 1, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- correlationVecchia::baseline_1_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b2") {
    
    out.fisher      <- fit_spacetime_baseline_range(z, locs, ms = c(m), baseline = 2, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- correlationVecchia::baseline_2_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b3") {
    
    out.fisher      <- fit_spacetime_baseline_range(z, locs, ms = c(m), baseline = 3, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- correlationVecchia::baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = covparms.est)
    
  } else if(approx == "cc") {
    
    out.fisher      <- fit_spacetime_corrvecchia_range(z, locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = covparms.est) + nugget * diag(n)
    vecchia.approx  <- correlationVecchia::corrvecchia_specify_knownCovparms_2(locs = locs, m = m, covmodel = Sigma.true, covparms = covparms.est)
    
  } else if(approx == "ex") {
    
    out.fisher      <- fit_spacetime_exactGP_range(z, locs, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- correlationVecchia::baseline_1_spacetime_specify(locs = locs, m = n - 1)
    
  } else if(approx == "tr") {
    
    covparms.est    <- covparms
    
    vecchia.approx  <- correlationVecchia::baseline_1_spacetime_specify(locs = locs, m = n - 1)
    
  } else {
    
    stop("Something wrong with approx!")
    
  }
  
  ### posterior
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    print(i)
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], alpha, covparms[3], covparms[4])) + nugget * diag(n)
    
    ## exact
    loglik          <- dmvnorm(z, mean = rep(0, n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = c(1), nuggets = nugget, covmodel = Sigma.ord)$U
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.vecchia, log = TRUE)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  ### return
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}



.posterior_fisherord_for_temporalRange <- function(z, locs, covparms, nugget, m, approx, N, xlim, sdlog)
{
  ### prior information
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[3])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  ### ordering / conditioning with fisher-scoring approach
  if(approx == "b1") {
    
    out.fisher      <- fit_spacetime_baseline_range(z, locs, ms = c(m), baseline = 1, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- correlationVecchia::baseline_1_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b2") {
    
    out.fisher      <- fit_spacetime_baseline_range(z, locs, ms = c(m), baseline = 2, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- correlationVecchia::baseline_2_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b3") {
    
    out.fisher      <- fit_spacetime_baseline_range(z, locs, ms = c(m), baseline = 3, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- correlationVecchia::baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = covparms.est)
    
  } else if(approx == "cc") {
    
    out.fisher      <- fit_spacetime_corrvecchia_range(z, locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = covparms.est) + nugget * diag(n)
    vecchia.approx  <- correlationVecchia::corrvecchia_specify_knownCovparms_2(locs = locs, m = m, covmodel = Sigma.true, covparms = covparms.est)
    
  } else if(approx == "ex") {
    
    out.fisher      <- fit_spacetime_exactGP_range(z, locs, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = 0, print.level = 0, max.it = 100, tol.dec = 4)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- correlationVecchia::baseline_1_spacetime_specify(locs = locs, m = nrow(locs) - 1)
    
  } else if(approx == "tr") {
    
    covparms.est    <- covparms
    
    vecchia.approx  <- correlationVecchia::baseline_1_spacetime_specify(locs = locs, m = n - 1)
    
  } else {
    
    stop("Something wrong with approx!")
    
  }
  
  ### posterior
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    print(i)
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], covparms[2], alpha, covparms[4])) + nugget * diag(n)
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = c(1), nuggets = nugget, covmodel = Sigma.ord)$U
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.vecchia, log = TRUE)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  ### return
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}

posterior_fisherord_for_cov_matern_spacetime <- function(z, locs, covparms, nugget = 0, m = 30, approx = "cc", target = "spatialRange", N = 100, xlim = c(0.05, 0.15), sdlog = 0.6)
{
  if(target == "spatialRange" | target == "spatial") {
    
    output <- .posterior_fisherord_for_spatialRange(z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, approx = approx, N = N, xlim = xlim, sdlog = sdlog)
    
  } else if(target == "temporalRange" | target == "temporal") {
    
    output <- .posterior_fisherord_for_temporalRange(z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, approx = approx, N = N, xlim = xlim, sdlog = sdlog)
    
  } else {
    
    stop("Check the target!")
    
  }
  
  return(output)
}

### simulation ########################################################

# n               <- 30^2
# m               <- 30
# covparms        <- c(1, 0.1, 1.0, 0.5)
# process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)
# locs            <- process$sim$sim1$locs
# 
# nugget          <- 0.0001
# z               <- process$sim$sim1$y + nugget * rnorm(n)
# 
# approx          <- "b1"
# 
# N = 100 ; xlim = c(0.05, 0.15) ; sdlog = 0.6

### n = 900, monitoring station #######################################

process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "space.random.time.grid", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs

z               <- process$sim$sim1$y + nugget * rnorm(n) ; fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.spatial.n900.monitoring <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_fisherord_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.05, 0.12), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.temporal.n900.monitoring <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_fisherord_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.05, 3.05), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

### n = 900, satellite ################################################

process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs

z               <- process$sim$sim1$y + nugget * rnorm(n) # fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.spatial.n900.satellite <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_fisherord_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.08, 0.12), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.temporal.n900.satellite <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_fisherord_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.05, 3.05), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

### n = 900, all random locs ##########################################

process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs

z               <- process$sim$sim1$y + nugget * rnorm(n) ; fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.spatial.n900.allrandom <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_fisherord_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.08, 0.16), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.temporal.n900.allrandom <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_fisherord_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.05, 2.25), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

##### Save #########################################################################

save(candid, n, covparms, nugget, out.spatial.n900.allrandom, out.temporal.n900.allrandom, out.spatial.n900.monitoring, out.temporal.n900.monitoring, out.spatial.n900.satellite, out.temporal.n900.satellite, file = "bayesian_fisherord_05132021.RData")

##### Visualization  1 #############################################################

par(mfrow = c(6, 5))
for(i in 1:length(out.spatial.n900.allrandom)) {
  temp <- out.spatial.n900.allrandom[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(6, 5))
for(i in 1:length(out.temporal.n900.allrandom)) {
  temp <- out.temporal.n900.allrandom[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(6, 5))
for(i in 1:length(out.spatial.n900.monitoring)) {
  temp <- out.spatial.n900.monitoring[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(6, 5))
for(i in 1:length(out.temporal.n900.monitoring)) {
  temp <- out.temporal.n900.monitoring[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(6, 5))
for(i in 1:length(out.spatial.n900.satellite)) {
  temp <- out.spatial.n900.satellite[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(6, 5))
for(i in 1:length(out.temporal.n900.satellite)) {
  temp <- out.temporal.n900.satellite[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

##### Visualization  3 #############################################################

m <- unique(candid$m)
cols <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70", "gray30")
legends <- c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Optimal GP", "True GP")

idx.list <- list()
for(i in 1:length(m)) idx.list[[i]] <- which(candid$m == m[i])


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.spatial.n900.allrandom[[idx.list[[i]][1]]]$alpha, out.spatial.n900.allrandom[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][5]
  k <- 5
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")
  
  j <- idx.list[[i]][6]
  k <- 6
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")
  
}
par(mfrow = c(1, 1))


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.temporal.n900.allrandom[[idx.list[[i]][1]]]$alpha, out.temporal.n900.allrandom[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][5]
  k <- 5
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")
  
  j <- idx.list[[i]][6]
  k <- 6
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")
  
}
par(mfrow = c(1, 1))


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.spatial.n900.monitoring[[idx.list[[i]][1]]]$alpha, out.spatial.n900.monitoring[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topleft", legend = legends, col = cols, lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "twodash")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "longdash")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotdash")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dashed")
  
  j <- idx.list[[i]][5]
  k <- 5
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotted")
  
  j <- idx.list[[i]][6]
  k <- 6
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotted")
  
}
par(mfrow = c(1, 1))


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.temporal.n900.monitoring[[idx.list[[i]][1]]]$alpha, out.temporal.n900.monitoring[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "twodash")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "longdash")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotdash")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dashed")
  
  j <- idx.list[[i]][5]
  k <- 5
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotted")
  
  j <- idx.list[[i]][6]
  k <- 6
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotted")
  
}
par(mfrow = c(1, 1))


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.spatial.n900.satellite[[idx.list[[i]][1]]]$alpha, out.spatial.n900.satellite[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "twodash")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "longdash")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotdash")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 1, lty = "solid")
  
  j <- idx.list[[i]][5]
  k <- 5
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotted")
  
  j <- idx.list[[i]][6]
  k <- 6
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotted")
  
}
par(mfrow = c(1, 1))


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.temporal.n900.satellite[[idx.list[[i]][1]]]$alpha, out.temporal.n900.satellite[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "twodash")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "longdash")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotdash")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 1, lty = "solid")
  
  j <- idx.list[[i]][5]
  k <- 5
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 1, lty = "dotted")
  
  j <- idx.list[[i]][6]
  k <- 6
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 1, lty = "dotted")
  
}
par(mfrow = c(1, 1))
