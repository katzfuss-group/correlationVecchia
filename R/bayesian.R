####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to obtain posterior densities.
###
###   Contents: 
###     posterior_mulv / .posterior_mulv_range / .posterior_mulv_latentDim
###     posterior_spacetime / .posterior_spacetime_srange / .posterior_spacetime_trange
###     posterior_ordfix_mulv / .posterior_ordfix_mulv_range / .posterior_ordfix_mulv_latentDim
###     posterior_ordfix_spacetime / .posterior_ordfix_spacetime_srange / .posterior_ordfix_spacetime_trange
###
####################################################################################

# posterior_mulv <- function()
# {}
# 
# .posterior_mulv_range <- function()
# {}
# 
# .posterior_mulv_latentDim <- function()
# {}


#' @title Computing the approximate posterior of the spatial/temporal range parameter for cov_matern_spacetime
#'
#' @param target srange or trange
#' @param approx b1, b2, b3, or cc
#' @param method.U "simple" or "ic0"
#' @param z process
#' @param locs locs
#' @param covparms covariance parameters
#' @param nugget nugget (noise)
#' @param m the size of conditioning sets 
#' @param N the number of grid points
#' @param xlim extent of the search
#' @param sdlog prior standard deviation of log(range parameter) 
#' @param verbose logical
#'
#' @import RandomFields
#' @import mvtnorm
#' @import LaplacesDemon
#' 
#' @return list
#' 
#' @export
#'
#' @examples
#' 1 + 1
posterior_spacetime <- function(target = "srange", approx = "cc", method.U = "simple", z, locs, covparms, nugget = 0, m = 30, N = 100, xlim = c(0.05, 0.15), sdlog = 0.6, verbose = TRUE)
{
  if( target %in% c("srange", "spatial", "space", "s") ) {
    
    if( method.U %in% c("simple", "straightforward") ) {
      
      output <- .posterior_spacetime_srange(approx = approx, z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, N = N, xlim = xlim, sdlog = sdlog, verbose = verbose)
      
    } else if( method.U %in% c("ic0", "involved", "alternative", "two-step") ) {
      
      output <- .posterior_ic0_spacetime_srange(approx = approx, z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, N = N, xlim = xlim, sdlog = sdlog, verbose = verbose)
      
    } else {
      
      stop("Check the argument method. It must be simple or ic0.")
    }
    
  } else if( target %in% c("trange", "temporal", "time", "t") ) {
    
    if( method.U %in% c("simple", "straightforward") ) {
      
      output <- .posterior_spacetime_trange(approx = approx, z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, N = N, xlim = xlim, sdlog = sdlog, verbose = verbose)
      
    } else if( method.U %in% c("ic0", "involved", "alternative", "two-step") ) {
      
      output <- .posterior_ic0_spacetime_trange(approx = approx, z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, N = N, xlim = xlim, sdlog = sdlog, verbose = verbose)
      
    } else {
      
      stop("Check the argument method. It must be simple or ic0.")
    }
    
  } else {
    
    stop("Check the argument target!")
  }
  
  return(output)
}

.posterior_spacetime_srange <- function(approx, z, locs, covparms, nugget, m, N, xlim, sdlog, verbose)
{
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[2])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    if(verbose == TRUE) { print(i) ; Sys.time() }
    
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], alpha, covparms[3], covparms[4])) + nugget * diag(n)
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    ## Vecchia
    if(approx == "b1") {
      
      vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = m)
      
    } else if(approx == "b2") {
      
      vecchia.approx  <- baseline_2_spacetime_specify(locs = locs, m = m)
      
    } else if(approx == "b3") {
      
      vecchia.approx  <- baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = c(covparms[1], alpha, covparms[3], covparms[4]))
    
    } else {
      
      vecchia.approx  <- cvecchia_m_specify(locs = locs, m = m, covmodel = Sigma.true, covparms = c(covparms[1], alpha, covparms[3], covparms[4]))
    }
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = covparms[1], nuggets = c(0), covmodel = Sigma.ord)$U
    
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.vecchia, log = TRUE) ; rm(Sigma.vecchia)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  # post            <- exp(logpost - mean(logpost))
  # post.norm       <- post
  # for(j in 1:2) post.norm[, j] <- post[, j] / sum(post[, j])
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}

.posterior_spacetime_trange <- function(approx, z, locs, covparms, nugget, m, N, xlim, sdlog, verbose)
{
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[3])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    if(verbose == TRUE) { print(i) ; Sys.time() }
    
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], covparms[2], alpha, covparms[4])) + nugget * diag(n)
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    ## Vecchia
    if(approx == "b1") {
      
      vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = m)
      
    } else if(approx == "b2") {
      
      vecchia.approx  <- baseline_2_spacetime_specify(locs = locs, m = m)
      
    } else if(approx == "b3") {
      
      vecchia.approx  <- baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = c(covparms[1], covparms[2], alpha, covparms[4]))
      
    } else {
      
      vecchia.approx  <- cvecchia_m_specify(locs = locs, m = m, covmodel = Sigma.true, covparms = c(covparms[1], covparms[2], alpha, covparms[4]))
    }
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = covparms[1], nuggets = c(0), covmodel = Sigma.ord)$U
    
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.vecchia, log = TRUE) ; rm(Sigma.vecchia)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  # post            <- exp(logpost - mean(logpost))
  # post.norm       <- post
  # for(j in 1:2) post.norm[, j] <- post[, j] / sum(post[, j])
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}

.posterior_ic0_spacetime_srange <- function(approx, z, locs, covparms, nugget, m, N, xlim, sdlog, verbose)
{
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[2])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    if(verbose == TRUE) { print(i) ; Sys.time() }
    
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], alpha, covparms[3], covparms[4])) 
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true + nugget * diag(n), log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    ## Vecchia
    if(approx == "b1") {
      
      vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = m)
      
    } else if(approx == "b2") {
      
      vecchia.approx  <- baseline_2_spacetime_specify(locs = locs, m = m)
      
    } else if(approx == "b3") {
      
      vecchia.approx  <- baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = c(covparms[1], alpha, covparms[3], covparms[4]))
      
    } else {
      
      vecchia.approx  <- cvecchia_m_specify(locs = locs, m = m, covmodel = Sigma.true, covparms = c(covparms[1], alpha, covparms[3], covparms[4]))
    }
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    
    vecchia.approx$ic0 <- TRUE
    U.obj           <- GPvecchia::createU(vecchia.approx, covparms = covparms[1], nuggets = nugget, covmodel = Sigma.ord)
    loglik          <- GPvecchia:::vecchia_likelihood_U(z, U.obj)

    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}

.posterior_ic0_spacetime_trange <- function(approx, z, locs, covparms, nugget, m, N, xlim, sdlog, verbose)
{
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[3])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    if(verbose == TRUE) { print(i) ; Sys.time() }
    
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], covparms[2], alpha, covparms[4])) + nugget * diag(n)
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    ## Vecchia
    if(approx == "b1") {
      
      vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = m)
      
    } else if(approx == "b2") {
      
      vecchia.approx  <- baseline_2_spacetime_specify(locs = locs, m = m)
      
    } else if(approx == "b3") {
      
      vecchia.approx  <- baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = c(covparms[1], covparms[2], alpha, covparms[4]))
      
    } else {
      
      vecchia.approx  <- cvecchia_m_specify(locs = locs, m = m, covmodel = Sigma.true, covparms = c(covparms[1], covparms[2], alpha, covparms[4]))
    }
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]

    vecchia.approx$ic0 <- TRUE
    U.obj           <- GPvecchia::createU(vecchia.approx, covparms = covparms[1], nuggets = nugget, covmodel = Sigma.ord)
    loglik          <- GPvecchia:::vecchia_likelihood_U(z, U.obj)
    
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}

####################################################################################

# posterior_ordfix_mulv <- function()
# {}
# 
# .posterior_ordfix_mulv_range <- function()
# {}
# 
# .posterior_ordfix_mulv_latentDim <- function()
# {}


#' @title Computing the approximate posterior of the spatial/temporal range parameter for cov_matern_spacetime
#'
#' @param target srange or trange
#' @param approx b1, b2, b3, cc, ex, or tr
#' @param method.U "simple" or "ic0"
#' @param z process
#' @param locs locs
#' @param covparms covariance parameters
#' @param nugget nugget (noise)
#' @param m the size of conditioning sets 
#' @param N the number of grid points
#' @param xlim extent of the search
#' @param sdlog prior standard deviation of log(range parameter) 
#' @param verbose logical
#' @param tol.dec tol.dec
#'
#' @import RandomFields
#' @import mvtnorm
#' @import LaplacesDemon
#' 
#' @return list
#' 
#' @export
#'
#' @examples
#' 1 + 1
posterior_ordfix_spacetime <- function(target = "srange", approx = "cc", method.U = "simple", z, locs, covparms, nugget = 0, m = 30, N = 100, xlim = c(0.05, 0.15), sdlog = 0.6, verbose = TRUE, tol.dec = 4)
{
  if( target %in% c("srange", "spatial", "space", "s") ) {
    
    if( method.U %in% c("simple", "straightforward") ) {
      
      output <- .posterior_ordfix_spacetime_srange(approx = approx, z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, N = N, xlim = xlim, sdlog = sdlog, verbose = verbose, tol.dec = tol.dec)
      
    } else if( method.U %in% c("ic0", "involved", "alternative", "two-step") ) {
      
      output <- .posterior_ordfix_ic0_spacetime_srange(approx = approx, z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, N = N, xlim = xlim, sdlog = sdlog, verbose = verbose, tol.dec = tol.dec)
      
    } else {
      
      stop("Check the argument method. It must be simple or ic0.")
    }
    
  } else if( target %in% c("trange", "temporal", "time", "t") ) {
    
    if( method.U %in% c("simple", "straightforward") ) {
      
      output <- .posterior_ordfix_spacetime_trange(approx = approx, z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, N = N, xlim = xlim, sdlog = sdlog, verbose = verbose, tol.dec = tol.dec)
      
    } else if( method.U %in% c("ic0", "involved", "alternative", "two-step") ) {
      
      output <- .posterior_ordfix_ic0_spacetime_trange(approx = approx, z = z, locs = locs, covparms = covparms, nugget = nugget, m = m, N = N, xlim = xlim, sdlog = sdlog, verbose = verbose, tol.dec = tol.dec)
      
    } else {
      
      stop("Check the argument method. It must be simple or ic0.")
    }
    
  } else {
    
    stop("Check the argument target!")
  }
  
  return(output)
}

.posterior_ordfix_spacetime_srange <- function(approx, z, locs, covparms, nugget, m, N, xlim, sdlog, verbose, tol.dec)
{
  ### prior information
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[2])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  ### ordering / conditioning with fisher-scoring approach
  if(approx == "b1") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 1, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b2") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 2, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_2_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b3") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 3, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = covparms.est)

  } else if(approx == "cc") {
    
    out.fisher      <- .fit_scaled_spacetime_cvecchia_range(y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]

    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = covparms.est) + nugget * diag(n)
    vecchia.approx  <- cvecchia_m_specify(locs = locs, m = m, covmodel = Sigma.true, covparms = covparms.est, abs.corr = FALSE, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor")

  } else if(approx == "ex") {
    
    out.fisher      <- .fit_scaled_spacetime_exactgp_range(y = z, inputs = locs, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = n-1)
    
  } else if(approx == "tr") {
    
    covparms.est    <- covparms
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = n-1)
    
  } else {
    
    stop("Something wrong with approx!")
    
  }
  
  ### posterior
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    if(verbose == TRUE) { print(i) ; Sys.time() }
    
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], alpha, covparms[3], covparms[4])) + nugget * diag(n)
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = c(1), nuggets = c(0), covmodel = Sigma.ord)$U
    
    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.vecchia, log = TRUE) ; rm(Sigma.vecchia)
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

.posterior_ordfix_spacetime_trange <- function(approx, z, locs, covparms, nugget, m, N, xlim, sdlog, verbose, tol.dec)
{
  ### prior information
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[3])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  ### ordering / conditioning with fisher-scoring approach
  if(approx == "b1") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 1, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b2") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 2, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_2_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b3") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 3, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = covparms.est)
    
  } else if(approx == "cc") {
    
    out.fisher      <- .fit_scaled_spacetime_cvecchia_range(y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = nugget, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = covparms.est) + nugget * diag(n)
    vecchia.approx  <- cvecchia_m_specify(locs = locs, m = m, covmodel = Sigma.true, covparms = covparms.est, abs.corr = FALSE, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor")
    
  } else if(approx == "ex") {
    
    out.fisher      <- .fit_scaled_spacetime_exactgp_range(y = z, inputs = locs, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = n-1)
    
  } else if(approx == "tr") {
    
    covparms.est    <- covparms
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = n-1)
    
  } else {
    
    stop("Something wrong with approx!")
    
  }
  
  ### posterior
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    if(verbose == TRUE) { print(i) ; Sys.time() }
    
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], covparms[2], alpha, covparms[4])) + nugget * diag(n)
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true, log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]
    U               <- GPvecchia::createU(vecchia.approx, covparms = c(1), nuggets = c(0), covmodel = Sigma.ord)$U

    revord          <- order(vecchia.approx$ord)
    Sigma.vecchia   <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.vecchia, log = TRUE) ; rm(Sigma.vecchia)
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  # post            <- exp(logpost - mean(logpost))
  # post.norm       <- post
  # for(j in 1:2) post.norm[, j] <- post[, j] / sum(post[, j])
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}

.posterior_ordfix_ic0_spacetime_srange <- function(approx, z, locs, covparms, nugget, m, N, xlim, sdlog, verbose, tol.dec)
{
  ### prior information
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[2])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  ### ordering / conditioning with fisher-scoring approach
  if(approx == "b1") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 1, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b2") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 2, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_2_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b3") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 3, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = covparms.est)
    
  } else if(approx == "cc") {
    
    out.fisher      <- .fit_scaled_spacetime_cvecchia_range(y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = covparms.est) + nugget * diag(n)
    vecchia.approx  <- cvecchia_m_specify(locs = locs, m = m, covmodel = Sigma.true, covparms = covparms.est, abs.corr = FALSE, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor")
    
  } else if(approx == "ex") {
    
    out.fisher      <- .fit_scaled_spacetime_exactgp_range(y = z, inputs = locs, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = n-1)
    
  } else if(approx == "tr") {
    
    covparms.est    <- covparms
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = n-1)
    
  } else {
    
    stop("Something wrong with approx!")
    
  }
  
  ### posterior
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    if(verbose == TRUE) { print(i) ; Sys.time() }
    
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], alpha, covparms[3], covparms[4]))
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true + nugget * diag(n), log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]

    vecchia.approx$ic0 <- TRUE
    U.obj           <- GPvecchia::createU(vecchia.approx, covparms = covparms[1], nuggets = nugget, covmodel = Sigma.ord)
    loglik          <- GPvecchia:::vecchia_likelihood_U(z, U.obj)
    
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

.posterior_ordfix_ic0_spacetime_trange <- function(approx, z, locs, covparms, nugget, m, N, xlim, sdlog, verbose, tol.dec)
{
  ### prior information
  n               <- nrow(locs)
  
  meanlog         <- log(covparms[3])
  logprior        <- function(alpha) dlnorm(alpha, meanlog, sdlog, log = TRUE)
  
  ### ordering / conditioning with fisher-scoring approach
  if(approx == "b1") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 1, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b2") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 2, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_2_spacetime_specify(locs = locs, m = m)
    
  } else if(approx == "b3") {
    
    out.fisher      <- .fit_scaled_spacetime_baseline_range(approx = 3, y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_3_spacetime_specify(locs = locs, m = m, covmodel = cov_matern_spacetime, covparms = covparms.est)
    
  } else if(approx == "cc") {
    
    out.fisher      <- .fit_scaled_spacetime_cvecchia_range(y = z, inputs = locs, ms = c(m), sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = covparms.est) + nugget * diag(n)
    vecchia.approx  <- cvecchia_m_specify(locs = locs, m = m, covmodel = Sigma.true, covparms = covparms.est, abs.corr = FALSE, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor")
    
  } else if(approx == "ex") {
    
    out.fisher      <- .fit_scaled_spacetime_exactgp_range(y = z, inputs = locs, sig2 = covparms[1], srange.ini = covparms[2], trange.ini = covparms[3], nu = covparms[4], nug = nugget, print.level = 0, max.it = 100, tol.dec = tol.dec)
    covparms.est    <- out.fisher$covparms[1:4]
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = n-1)
    
  } else if(approx == "tr") {
    
    covparms.est    <- covparms
    
    vecchia.approx  <- baseline_1_spacetime_specify(locs = locs, m = n-1)
    
  } else {
    
    stop("Something wrong with approx!")
    
  }
  
  ### posterior
  alpha.grid      <- seq(xlim[1], xlim[2], length = N)
  logpost         <- array(dim = c(N, 2))
  for(i in 1:N) {
    
    if(verbose == TRUE) { print(i) ; Sys.time() }
    
    alpha           <- alpha.grid[i]
    
    Sigma.true      <- cov_matern_spacetime(locs = locs, covparms = c(covparms[1], covparms[2], alpha, covparms[4]))
    
    ## exact
    loglik          <- mvtnorm::dmvnorm(z, mean = rep(0, n), sigma = Sigma.true + nugget * diag(n), log = TRUE)
    logpost[i, 1]   <- loglik + logprior(alpha)
    
    Sigma.ord       <- Sigma.true[vecchia.approx$ord, vecchia.approx$ord]

    vecchia.approx$ic0 <- TRUE
    U.obj           <- GPvecchia::createU(vecchia.approx, covparms = covparms[1], nuggets = nugget, covmodel = Sigma.ord)
    loglik          <- GPvecchia:::vecchia_likelihood_U(z, U.obj)
    
    logpost[i, 2]   <- loglik + logprior(alpha)
  }
  
  post            <- logpost - mean(logpost)
  A               <- mean(post[, 1])
  B               <- mean(post[, 2])
  
  post.norm       <- post
  post.norm[, 1]  <- exp(post[, 1] - A) / sum(exp(post[, 1] - A))
  post.norm[, 2]  <- exp(post[, 2] - B) / sum(exp(post[, 2] - B))
  
  post            <- exp(post)
  
  return(list(alpha = alpha.grid, logpost = logpost, post = post, post.norm = post.norm))
}

