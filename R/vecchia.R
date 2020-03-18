####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to implement vecchia approximations
###
###   Contents:
###       FTN - corrvecchia_specify_knownCovparms() consisting of .distance_correlation() and correlation()
###       FTN - order_coordinate() 
###       FTN - order_maxmin_euclidean() 
###       FTN - order_maxmin_correlation_straightforward() 
###       FTN - order_maxmin_correlation_inverseDist() 
###       FTN - conditioning_nn()
###
####################################################################################



#' @title specify a general vecchia approximation with known parameters
#'
#' @param locs A matrix of locations
#' @param m Number of nearby points to condition on
#' @param ordering 'maxmin' or 'coord'
#' @param ordering.method either 'euclidean' or 'correlation'
#' @param coordinate integer or vector of integers in \code{1,...,d}
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|rho|. If \code{FALSE} then distane = 1-rho
#' @param initial.pt NULL = which.min(rowMeans(d)), center = euclidean-based center, integer = specify the first obs, 'random' = at random, and else = which.min(rowMeans(d))
#' @param conditioning 'NN' (nearest neighbor)
#' @param conditioning.method either 'euclidean' or 'correlation'
#' @param covmodel covariance function (or matrix)
#' @param covparms covariance parameters as a vector. The first entry must be its overall variance
#'
#' @return An object that specifies the vecchia approximation for later use in likelihood evaluation or prediction.
#' @export
#'
#' @examples
#' n             <- 15^2
#' m             <- 15
#' locs          <- matrix(runif(n * 2, 0, 1), n, 2)
#' covparms      <- c(1, 0.1, 10)
#' 
#' # true cov matrix
#' Sigma <- cov_expo_aniso(locs, covparms)
#' 
#' # Visualize the process
#' y <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
#' fields::quilt.plot(locs[,1], locs[,2], y)
#' 
#' out1 <- corrvecchia_specify_knownCovparms(locs = locs, m = m, 
#'                                           ordering = "maxmin", ordering.method = "euclidean", 
#'                                           coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, 
#'                                           conditioning = "NN", conditioning.method = "euclidean", 
#'                                           covmodel = cov_expo_aniso, covparms = covparms)
#' out2 <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", 
#'                                           ordering.method = "euclidean", 
#'                                           coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, 
#'                                           conditioning = "NN", conditioning.method = "correlation", 
#'                                           covmodel = cov_expo_aniso, covparms = covparms)
#' out3 <- corrvecchia_specify_knownCovparms(locs = locs, m = m, ordering = "maxmin", 
#'                                           ordering.method = "correlation", 
#'                                           coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, 
#'                                           conditioning = "NN", conditioning.method = "euclidean", 
#'                                           covmodel = cov_expo_aniso, covparms = covparms)
#' out4 <- corrvecchia_specify_knownCovparms(locs = locs, m = m, 
#'                                           ordering = "maxmin", ordering.method = "correlation", 
#'                                           coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, 
#'                                           conditioning = "NN", conditioning.method = "correlation", 
#'                                           covmodel = cov_expo_aniso, covparms = covparms)
#' out5 <- corrvecchia_specify_knownCovparms(locs = locs, m = m, 
#'                                           ordering = "coord", ordering.method = "euclidean", 
#'                                           coordinate = NULL, abs.corr = FALSE, initial.pt = NULL, 
#'                                           conditioning = "NN", conditioning.method = "euclidean", 
#'                                           covmodel = cov_expo_aniso, covparms = covparms)
#' out6 <- corrvecchia_specify_knownCovparms(locs = locs, m = m, 
#'                                           ordering = "coord", ordering.method = "euclidean", 
#'                                           coordinate = c(1), abs.corr = FALSE, initial.pt = NULL, 
#'                                           conditioning = "NN", conditioning.method = "euclidean", 
#'                                           covmodel = cov_expo_aniso, covparms = covparms)
#' out7 <- corrvecchia_specify_knownCovparms(locs = locs, m = m, 
#'                                           ordering = "coord", ordering.method = "euclidean", 
#'                                           coordinate = c(2), abs.corr = FALSE, initial.pt = NULL, 
#'                                           conditioning = "NN", conditioning.method = "euclidean", 
#'                                           covmodel = cov_expo_aniso, covparms = covparms)
#' 
#' kls.euc.euc    <- performance(out = out1, locs = locs, 
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.cor.euc    <- performance(out = out2, locs = locs, 
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.euc.cor    <- performance(out = out3, locs = locs, 
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.cor.cor    <- performance(out = out4, locs = locs, 
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.coord.xy   <- performance(out = out5, locs = locs, 
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.coord.x    <- performance(out = out6, locs = locs, 
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.coord.y    <- performance(out = out7, locs = locs, 
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' 
#' c(kls.euc.euc, kls.euc.cor, kls.cor.euc, kls.cor.cor, kls.coord.xy, kls.coord.x, kls.coord.y)
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



#' @title Coordinate-based Ordering of locations
#'
#' @param locs A matrix of locations
#' @param coordinate integer or vector of integers in \code{1,...,d}
#'
#' @return A vector of indices giving the coordinate-based ordering
#' @export
#'
#' @examples
#' locs <- matrix(runif(100 * 2), 100, 2)
#' 
#' identical(order_coordinate(locs, NULL), GpGp::order_coordinate(locs, NULL))
#' identical(order_coordinate(locs, 1), GpGp::order_coordinate(locs, 1))
#' identical(order_coordinate(locs, 2), GpGp::order_coordinate(locs, 2))
#' identical(order_coordinate(locs, c(1, 2)), GpGp::order_coordinate(locs, c(1, 2)))
order_coordinate <- function(locs, coordinate) 
{
  if(ncol(locs) == 1){
    ord <- order(rowSums(locs[, drop = FALSE]))
  } else {
    ord <- order(rowSums(locs[, coordinate, drop = FALSE]))
  }
  
  return(ord)
}



#' @title Maximum-minimum (maxmin) ordering of locations in terms of the Euclidean distance
#'
#' @param locs A matrix of locations
#'
#' @return A vector of indices giving the euclidean-based maxmin ordering 
#' @export
#'
#' @examples
#' locs <- matrix(runif(100 * 2), 100, 2)
#' 
#' identical(order_maxmin_euclidean(locs), GPvecchia::order_maxmin_exact(locs))
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



#' @title Straightforward maximum-minimum (maxmin) ordering of locations with respect to a user-defined distance matrix
#'
#' @param locs A matrix of locations
#' @param d A matrix of distances between locations (for instance, 1-rho)
#' @param initial.pt NULL = which.min(rowMeans(d)), center = euclidean-based center, integer = specify the first obs, 'random' = at random, and else = which.min(rowMeans(d))
#'
#' @return A vector of indices giving the maxmin ordering with respect to the user-defined distance matrix
#' @export
#'
#' @examples
#' covparms <- c(1, 0.1)
#' locs <- matrix(runif(100 * 2), 100, 2)
#'
#' rho <- vecchia:::.correlation(locs, cov_expo_iso, covparms, FALSE)
#' identical(order_maxmin_euclidean(locs), 
#'           order_maxmin_correlation_straightforward(locs, 1 - rho, 
#'                                                    initial.pt = order_maxmin_euclidean(locs)[1]))
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



#' @title Reversely implemented Maximum-minimum (maxmin) ordering of locations with respect to a user-defined distance matrix
#'
#' @param locs A matrix of locations
#' @param d.inv A matrix of inverse of distances between locations (for instance, rho) 
#' @param initial.pt NULL = which.min(rowMeans(d)), center = euclidean-based center, integer = specify the first obs, 'random' = at random, and else = which.min(rowMeans(d))
#'
#' @return A vector of indices giving the correlation-based maxmin ordering 
#' @export
#'
#' @examples
#' covparms <- c(1, 0.1)
#' locs <- matrix(runif(100 * 2), 100, 2)
#' 
#' rho <- vecchia:::.correlation(locs, cov_expo_iso, covparms, FALSE)
#' identical(order_maxmin_correlation_inverseDist(locs, rho, 
#'                        initial.pt = GPvecchia::order_maxmin_exact(locs)[1]), 
#'           order_maxmin_correlation_straightforward(locs, 1 - rho, 
#'                        initial.pt = GPvecchia::order_maxmin_exact(locs)[1]))
#' 
#' covparms <- c(1, 0.1, 10)
#' locs.trans <- cbind(locs[ ,1] * covparms[3], locs[,2]) / covparms[2]
#' 
#' rho.new <- vecchia:::.correlation(locs, cov_expo_aniso, covparms, FALSE)
#' identical(order_maxmin_euclidean(locs.trans),
#'           order_maxmin_correlation_inverseDist(locs, rho.new, 
#'                        initial.pt = order_maxmin_euclidean(locs.trans)[1]))
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



#' @title Nearest-Neighbor (NN) conditioning with respect to user-specified distance matrix
#'
#' @param m Number of nearby points to condition on
#' @param d A matrix of distances between locations (for instance, 1-rho)
#'
#' @return A matrix of indices giving NN conditioning sets
#' @export
#'
#' @examples
#' locs <- matrix(runif(100 * 2), 100, 2)
#' d <- as.matrix(dist(locs))
#' which(conditioning_nn(5, d) != GpGp::find_ordered_nn_brute(locs, 5), arr.ind = TRUE)
#' which(conditioning_nn(5, d) != GpGp::find_ordered_nn(locs, 5), arr.ind = TRUE)
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
