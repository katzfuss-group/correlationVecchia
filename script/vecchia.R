corrvecchia_specify_knownCovparms <- function(locs, m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel, covparms) 
{
  locs  <- as.matrix(locs)
  p     <- ncol(locs)
  n     <- nrow(locs)
  
  if(ordering.method == "correlation" | conditioning.method == "correlation") {
    
    rho <- .correlation(locs = locs, covmodel = covmodel, covparms = covparms, abs.corr = abs.corr)
    
  }
  
  # ordering
  if(ordering == "coord") {
    
    if(is.null(coordinate)) coordinate <- seq(p)
    ord         <- order_coordinate(locs = locs, coordinate = coordinate)
    
  } else if(ordering == "maxmin" & ordering.method == "euclidean") {
    
    ord         <- order_maxmin_euclidean(locs = locs)
    
  } else if(ordering == "maxmin" & ordering.method == "correlation") {
    
    ord         <- order_maxmin_correlation_inverseDist(locs = locs, d.inv = rho, initial.pt = initial.pt)
    
  } else {
    
    stop("Please check the ordering method.")
    
  }
  
  # order locations
  locsord       <- locs[ord, , drop = FALSE]
  
  if(is.matrix(covmodel)) {
    covmodel    <- covmodel[ord, ord]
  }
  
  # conditioning
  if(conditioning.method == "euclidean") {
    
    cond.sets   <- GpGp::find_ordered_nn(locs = locsord, m = m)
    
  } else if(conditioning.method == "correlation") {
    
    rho         <- rho[ord, ord]
    cond.sets   <- conditioning_nn(m = m, d = 1 - rho)
    
  } else {
    
    stop("Please check the conditioning method.")
    
  }
  
  # return 
  Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = ordering, ordering.method = ordering.method, conditioning = 'NN', conditioning.method = conditioning.method)
  return(vecchia.approx)
}

.distance_correlation <- function(locs, covmodel, covparms, abs.corr) 
{
  ### CAUTION: This function can cause numerical issue. Please use the 'correlation()' function, instead.
  
  if(is.function(covmodel)) {
    
    covparms[1] <- 1 # only use correlation
    if(abs.corr == FALSE) {
      dist.matrix <- 1 - covmodel(locs, covparms) # 1-rho
    } else if(abs.corr == TRUE) {
      dist.matrix <- 1 - abs(covmodel(locs, covparms)) # 1-|rho|
    } else {
      stop("abs.corr must be logical.")
    }
    
  } else if(is.matrix(covmodel)) {
    
    if(abs.corr == FALSE) {
      dist.matrix <- 1 - covmodel / covparms[1] # 1-rho
    } else if(abs.corr == TRUE) {
      dist.matrix <- 1 - abs(covmodel) / covparms[1] # 1-|rho|
    } else {
      stop("abs.corr must be logical.")
    }
    
  } else {
    
    stop("covmodel must be either function or matrix.")
    
  }
  
  return(dist.matrix)
}

.correlation <- function(locs, covmodel, covparms, abs.corr) 
{
  if(is.function(covmodel)) {
    
    covparms[1]   <- 1 # correlation function
    if(abs.corr == FALSE) {
      corr.matrix <- covmodel(locs, covparms) # rho
    }  else if(abs.corr == TRUE) {
      corr.matrix <- abs(covmodel(locs, covparms)) # |rho|
    } else {
      stop("abs.corr must be logical.")
    }
    
  } else if(is.matrix(covmodel)) {
    
    if(abs.corr == FALSE) {
      corr.matrix <- covmodel / covparms[1] # rho
    }  else if(abs.corr == TRUE) {
      corr.matrix <- abs(covmodel) / covparms[1] # |rho|
    } else {
      stop("abs.corr must be logical.")
    }
    
  } else {
    
    stop("covmodel must be either function or matrix.")
    
  }
  
  return(corr.matrix)
}

order_coordinate <- function(locs, coordinate) 
{
  if(ncol(locs) == 1){
    ord <- order(rowSums(locs[, drop = FALSE]))
  } else {
    ord <- order(rowSums(locs[, coordinate, drop = FALSE]))
  }
  
  return(ord)
}

order_maxmin_euclidean <- function(locs)
{
  n             <- nrow(locs)
  p             <- ncol(locs)
  ord           <- rep(NA, n)
  cen           <- t(as.matrix(colMeans(locs)))
  
  ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  cand.argmax   <- seq(n)[seq(n) != ord[1]]
  
  cdist         <- fields::rdist(locs[cand.argmax, ], matrix(locs[ord[1], ], nrow = 1, ncol = 2))
  ord[2]        <- cand.argmax[which.max(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]
  
  for(j in 3:(n-1)){
    cdist             <- fields::rdist(locs[cand.argmax, ], locs[ord[seq(j-1)], ])
    cdist             <- Rfast::rowMins(cdist, value = T)
    ord[j]            <- cand.argmax[which.max(cdist)]
    cand.argmax       <- cand.argmax[cand.argmax != ord[j]]
  } 
  
  ord[n]        <- cand.argmax
  return(ord)
}

order_maxmin_correlation_straightforward <- function(locs, d, initial.pt)
{
  n             <- nrow(locs)
  p             <- ncol(locs)
  ord           <- rep(NA, n)
  
  # first step
  if( is.null(initial.pt) ){
    ord[1]          <-  which.min(rowSums(d))
  } else if( is.numeric(initial.pt) ) {
    ord[1]          <- initial.pt
  } else if( initial.pt == 'center' ) {
    cen             <- t(as.matrix(colMeans(locs)))
    ord[1]          <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  } else if( initial.pt == 'random') { # at random
    ord[1]          <- sample(1:n, 1)
  } else {
    ord[1]          <- which.min(rowSums(d))
  }
  cand.argmax   <- seq(n)[seq(n) != ord[1]]
  
  ind           <- as.matrix(expand.grid(cand.argmax, ord[1]))
  cdist         <- d[ind]
  ord[2]        <- cand.argmax[which.max(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]
  
  for(j in 3:(n-1)){
    ind               <- as.matrix(expand.grid(cand.argmax, ord[seq(j-1)]))
    cdist             <- matrix(d[ind], nrow = length(cand.argmax), ncol = j-1, byrow = F)
    cdist             <- Rfast::rowMins(cdist, value = T)
    ord[j]            <- cand.argmax[which.max(cdist)]
    cand.argmax       <- cand.argmax[cand.argmax != ord[j]]
  } 
  
  ord[n]        <- cand.argmax
  return(ord)
}

order_maxmin_correlation_inverseDist <- function(locs, d.inv, initial.pt)
{
  n             <- nrow(locs)
  p             <- ncol(locs)
  ord           <- rep(NA, n)
  
  # first step
  if( is.null(initial.pt) ){
    ord[1]        <-  which.max(rowSums(d.inv))
  } else if( is.numeric(initial.pt) ) {
    ord[1]        <- initial.pt
  } else if( initial.pt == 'center' ) {
    cen           <- t(as.matrix(colMeans(locs)))
    ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  } else if( initial.pt == 'random') { # at random
    ord[1]        <- sample(1:n, 1)
  } else {
    ord[1]        <-  which.max(rowSums(d.inv))
  }
  cand.argmax   <- seq(n)[seq(n) != ord[1]]
  
  ind           <- as.matrix(expand.grid(cand.argmax, ord[1]))
  cdist         <- d.inv[ind]
  ord[2]        <- cand.argmax[which.min(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]
  
  for(j in 3:(n-1)){
    ind               <- as.matrix(expand.grid(cand.argmax, ord[seq(j-1)]))
    cdist             <- matrix(d.inv[ind], nrow = length(cand.argmax), ncol = j-1, byrow = F)
    cdist             <- Rfast::rowMaxs(cdist, value = T)
    ord[j]            <- cand.argmax[which.min(cdist)]
    cand.argmax       <- cand.argmax[cand.argmax != ord[j]]
  } 
  
  ord[n]        <- cand.argmax
  return(ord)
}

conditioning_nn <- function(m, d)
{
  # initialize an output object NN which is a n*n matrix
  n     <- nrow(d)
  NN    <- matrix(rep(NA, n * (m + 1)), nrow = n, ncol = m + 1) ; NN[1, 1] <- 1
  
  # Find the nearest neighbor conditioning set for each i-th location using the 'dist_to_knn()' function 
  for(i in 2:n) {
    k               <- min(i, m + 1) # the number of neighbors of the i-th observation
    NN[i, seq(k)]   <- scanstatistics::dist_to_knn(d[seq(i), seq(i)], k = k)[i, seq(k)]
  }
  
  return(NN)
}