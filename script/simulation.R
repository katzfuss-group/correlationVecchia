####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to conduct simulation.
###
####################################################################################



generate_gp <- function(nsim = 1, n = 15^2, d = 2, p = 1, same.locs = FALSE, covmodel, pivot = FALSE, method = NULL, tol = .Machine$double.eps, return.err = TRUE, verbose = TRUE, ...)
{
  ### random generation
  if(p == 1) {
    
    locs.full     <- matrix(runif(nsim * n * d, 0, 1), nsim * n, d)
    noise.full    <- rnorm(nsim * n)
    
  } else if(p > 1 & same.locs == FALSE) {
    
    locs.full     <- matrix(runif(nsim * n * p * d, 0, 1), nsim * n * p, d)
    noise.full    <- rnorm(nsim * n * p)
    
  } else if(p > 1 & same.locs == TRUE) {
    
    locs.full     <- matrix(runif(nsim * n * d, 0, 1), nsim * n, d)
    
    # ind <- unlist(lapply(lapply(seq(from = 1, to = nsim * n, by = n), function(x) seq(from = x, length.out = n)), function(x) rep(x = x, times = p)))
    # ind <- unlist(rep(lapply(seq(from = 1, to = nsim * n, by = n), function(x) seq(from = x, length.out = n)), each = p))
    ind           <- unlist(rep(split(seq(nsim * n), rep(seq(nsim), each = n)), each = p))
    
    locs.full     <- locs.full[ind, ]
    noise.full    <- rnorm(nsim * n * p)
    
  } else {
    
    stop("Either p or same.locs is invalid.")
    
  }
  
  ### assignment
  sim <- list()
  if(p == 1) {
    
    for(k in 1:nsim) {
      ind           <- seq(from = 1 + (k - 1) * n, to = k * n, by = 1)
      locs          <- locs.full[ind, ]
      
      covmat        <- covmodel(locs, ...)
      covfac        <- factorize(covmat = covmat, pivot = pivot, method = method, tol = tol, return.err = return.err, verbose = verbose)
      
      y             <- as.numeric(t(covfac$covfactor) %*% noise.full[ind])
      
      if(return.err == TRUE) decomp.err <- covfac$decomp.err
      
      # sim[[k]] <- ind
      sim[[k]]      <- list(locs = locs, y = y, covmat = covmat, covfac = covfac)
    }
    
  } else if(p > 1) {
    
    for(k in 1:nsim) {
      
      locs  <- list()
      for(i in 1:p) {
        ind           <- seq(from = 1 + (k - 1) * n * p + (i - 1) * n, to = (k - 1) * n * p + i * n, by = 1)
        
        # locs[[i]] <- ind
        locs[[i]]     <- locs.full[ind, ]
      }
      names(locs)   <- paste0("locs", seq(p))
      
      covmat        <- covmodel(locs, ...)
      covfac        <- factorize(covmat = covmat, pivot = pivot, method = method, tol = tol, return.err = return.err, verbose = verbose)
      
      ind           <- seq(from = 1 + (k - 1) * n * p, to = k * n * p, by = 1)
      # y <- ind
      y             <- as.numeric(t(covfac$covfactor) %*% noise.full[ind])
      
      sim[[k]]      <- list(locs = locs, y = y, covmat = covmat, covfac = covfac)
    }
    
  } else {
    
    stop("p is invalid.")
    
  }
  
  # return
  names(sim)  <- paste0("sim", seq(nsim))
  return(list(nsim = nsim, n = n, d = d, p = p, sim = sim))
}



simulate_corrvecchia_knownCovparms <- function(nsim = 1, n = 15^2, d = 2, p = 1, same.locs = FALSE, m = 10, covmodel, pivot = FALSE, method = NULL, tol = .Machine$double.eps, verbose = TRUE, ...)
{
  time.tot  <- proc.time()
  time.sim  <- list()
  time.gen  <- proc.time()
  
  ### generation
  if(verbose == TRUE) {
    message("")
    message("------------------------------------------------------------")
    message(paste0("-----                    ", "generation", "                    -----"))
    message("------------------------------------------------------------")
  }
  
  if(verbose == TRUE) message(paste0("System: Simulation starts. [", Sys.time(), "]"))
  
  realization   <- generate_gp(nsim = nsim, n = n, d = d, p = p, same.locs = same.locs, covmodel = covmodel, pivot = pivot, method = method, tol = tol, return.err = TRUE, verbose = verbose, ...)
  
  if(verbose == TRUE) message(paste0("System: Realizations of GPs are generated. [", Sys.time(), "]"))
  time.gen  <- proc.time() - time.gen
  
  ### simulation
  kls.full      <- matrix(NA, nsim, 20)
  decomp.err    <- rep(NA, nsim)
  if(p == 1) {
    
    for(k in 1:nsim) {
      
      time.sim[[k]] <- proc.time()
      
      if(verbose == TRUE) {
        message("")
        message("------------------------------------------------------------")
        message(paste0("-----                   ", "simulation ", k, "                   -----"))
        message("------------------------------------------------------------")
      }
      
      ### basic information
      locs              <- realization$sim[[k]]$locs
      covmat.modified   <- base::crossprod( realization$sim[[k]]$covfac$covfactor ) 
      decomp.err[k]     <- realization$sim[[k]]$covfac$decomp.err
      
      ### specify vecchia approximations
      approx        <- list()
      
      approx[[1]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: E-Maxmin + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      approx[[2]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: E-Maxmin + C-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      approx[[3]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: C-Maxmin + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      approx[[4]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: C-Maxmin + C-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      
      approx[[5]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "coord", ordering.method = "euclidean", coordinate = c(1), abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: X-coord + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      approx[[6]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "coord", ordering.method = "euclidean", coordinate = c(2), abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: Y-coord + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      
      ### compute approximate covariance matrices
      n.approx    <- length(approx)
      covmat.hat  <- list()
      kls         <- rep(NA, n.approx)
      for(i in 1:n.approx) {
        
        covmat.ord.modified   <- covmat.modified[approx[[i]]$ord, approx[[i]]$ord]
        U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified)$U
        
        # # previous version (wrong but noticeable)
        # covmat.ord            <- covmodel(approx[[i]]$locsord, ...)
        # covmat.ord.modified   <- modify(covmat.ord, pivot = pivot, method = method, tol = tol, return.err = TRUE, verbose = verbose)
        # U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified$covmat.modified)$U
        
        revord                <- order(approx[[i]]$ord)
        covmat.hat[[i]]       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
        
        kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]])
      }
      
      kls.full[k, 1:n.approx] <- kls
      if(verbose == TRUE) message(paste0("System: All vecchia approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
      
      time.sim[[k]] <- proc.time() - time.sim[[k]]
    }
    
  } else if(p > 1) {
    
    for(k in 1:nsim) {
      
      time.sim[[k]] <- proc.time()
      
      if(verbose == TRUE) {
        message("")
        message("------------------------------------------------------------")
        message(paste0("-----                   ", "simulation ", k, "                   -----"))
        message("------------------------------------------------------------")
      }
      
      ### basic information
      locs              <- do.call(rbind, realization$sim[[k]]$locs)
      covmat.modified   <- base::crossprod( realization$sim[[k]]$covfac$covfactor)
      decomp.err[k]     <- realization$sim[[k]]$covfac$decomp.err
      
      ### specify vecchia approximations
      approx        <- list()
      
      approx[[1]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: E-Maxmin + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      approx[[2]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: E-Maxmin + C-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      approx[[3]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: C-Maxmin + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      approx[[4]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: C-Maxmin + C-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      
      approx[[5]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "coord", ordering.method = "euclidean", coordinate = c(1), abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: X-coord + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      approx[[6]]   <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "coord", ordering.method = "euclidean", coordinate = c(2), abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = covmat.modified, covparms = c(1))
      if(verbose == TRUE) message(paste0("System: Y-coord + E-NN Vecchia approximation is accomplished. [", Sys.time(), "]"))
      
      ### compute approximate covariance matrices
      n.approx    <- length(approx)
      covmat.hat  <- list()
      kls         <- rep(NA, n.approx)
      for(i in 1:n.approx) {
        
        locsord <- list()
        for(j in 1:p) {
          locsord[[j]] <- approx[[i]]$locsord[seq(from = 1 + (j-1) * n, to = j * n, by = 1), , drop = FALSE]
        }
        
        covmat.ord.modified   <- covmat.modified[approx[[i]]$ord, approx[[i]]$ord]
        U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified)$U
        
        # # previous version (wrong but noticeable)
        # covmat.ord            <- covmodel(locsord, ...)
        # covmat.ord.modified   <- modify(covmat.ord, pivot = pivot, method = method, tol = tol, return.err = TRUE, verbose = verbose)
        # U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord.modified$covmat.modified)$U
        
        revord                <- order(approx[[i]]$ord)
        covmat.hat[[i]]       <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
        
        kls[i]                <- kldiv(covmat.modified, covmat.hat[[i]])
      }
      
      kls.full[k, 1:n.approx] <- kls
      if(verbose == TRUE) message(paste0("System: All vecchia approximations are evaluated in terms of KL divergence. [", Sys.time(), "]"))
      
      time.sim[[k]] <- proc.time() - time.sim[[k]]
    }
    
  } else {
    
    stop("p is invalid.")
    
  }
  
  kls.full <- data.frame(kls.full[, 1:n.approx, drop = FALSE])
  colnames(kls.full) <- paste0("approx_", 1:n.approx)
  
  kls.average <- colMeans(kls.full)
  names(kls.average) <- paste0("approx_", 1:n.approx)
  
  time.tot  <- proc.time() - time.tot
  names(time.sim) <- paste0("simulation_", 1:nsim)
  
  # return
  result              <- list()
  result$nsim         <- nsim
  result$n            <- n
  result$d            <- d
  result$p            <- p
  result$process      <- realization
  result$n.approx     <- n.approx
  result$kls.average  <- kls.average
  result$kls.full     <- kls.full
  result$decomp.error <- decomp.err
  result$time.tot     <- time.tot
  result$time.gen     <- time.gen
  result$time.sim     <- time.sim
  
  if(verbose == TRUE) {
    message("")
    message(paste0("System: Simulation is finished. [", Sys.time(), "]"))
  }
  
  return(result)
}