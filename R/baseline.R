####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to implement baseline models
###
###   Contents:
###
####################################################################################


#' @title The first baseline approximation for multivaraite cases
#'
#' @param locs A matrix of locations
#' @param m A size of conditioning sets
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' locs1  <- matrix(runif(20), 10, 2)
#' locs2  <- matrix(runif(30), 15, 2)
#' locs   <- list(locs1 = locs1, locs2 = locs2)
#' m      <- 5
#' 
#' # baseline 1 = separately E-maxmin ordering + half-and-half E-NN conditioning
#' # baseline 2 = separately E-maxmin ordering + unified E-NN conditioning 
#' # baseline 3 = separately E-maxmin ordering + separately E-NN conditioning
#' # baseline 4 = separately E-maxmin ordering + unified C-NN conditioning
#' 
#' # Note 1 = For all baseline approximations, the first process is not relevant to the second one
#' # Note 2 = For the third baseline approximations, two processes are totally irrelevant
#' 
#' ### ordering
#' ord1 <- baseline_1_multivariate_specify(locs = locs, m = m)$ord
#' ord2 <- baseline_2_multivariate_specify(locs = locs, m = m)$ord
#' ord3 <- baseline_3_multivariate_specify(locs = locs, m = m)$ord
#' ord4 <- baseline_4_multivariate_specify(locs = locs, m = m, 
#'                                         covmodel = cov_bivariate_flexMatern, 
#'                                         sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                         nu.mat = matrix(0.5, 2, 2), 
#'                                         alpha.mat = matrix(1, 2, 2))$ord
#' 
#' all.equal(ord1, ord2) ; all.equal(ord1, ord3) ; all.equal(ord1, ord4)
#' 
#' ### conditioning
#' ind <- seq(from = m, to = 1, by = -1)
#' cond1 <- baseline_1_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond2 <- baseline_2_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond3 <- baseline_3_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond4 <- baseline_4_multivariate_specify(locs = locs, m = m, 
#'                                          covmodel = cov_bivariate_flexMatern, 
#'                                          sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                          nu.mat = matrix(0.5, 2, 2), 
#'                                          alpha.mat = matrix(1, 2, 2))$U.prep$revNNarray[, ind]
#' 
#' cond1 == cond2 # similar
#' cond1 == cond3 # different
#' cond3 == cond4 # similar
baseline_1_multivariate_specify <- function(locs, m)
{
  p         <- length(locs)
  n         <- unlist(lapply(locs, nrow))
  n.cumsum  <- cumsum(n)
  
  ord <- list()
  for(i in 1:p) {
    ord[[i]]      <- order_maxmin_euclidean( locs[[i]] )
    if(i > 1) {
      ord.full      <- c(ord.full, ord[[i]] + n.cumsum[i-1])
    } else {
      ord.full      <- ord[[i]]
    }
  } 
  
  locsord     <- do.call(rbind, locs)[ord.full, , drop = FALSE]
  
  cond.sets   <- matrix(NA, nrow = nrow(locsord), ncol = m + 1)
  cond.sets[seq(n[1]), seq(m+1)] <- GpGp::find_ordered_nn(locs = locsord[seq(n[1]), , drop = FALSE], m = m)
  for(i in 2:p) {
    for(j in seq(from = n.cumsum[[i-1]] + 1, to = n.cumsum[[i]], by = 1)) {
      history       <- locsord[seq(j), , drop = FALSE]
      distance      <- as.numeric(fields::rdist(x1 = locsord[j, , drop = FALSE], x2 = history))
      ord.distance  <- order(distance)
      
      siz                 <- min(m + 1, j)
      siz.each            <- rep(NA, times = i)
      siz.each[i]         <- min( ceiling((m + 1) / p), j - n.cumsum[[i-1]] )
      siz.each[seq(i-1)]  <- table(rep_len(seq(i-1), length.out = siz - siz.each[i]))
      
      cond.ind <- list()
      for(k in 1:i) {
        if(k == 1) {
          ind <- seq(n.cumsum[[1]])
        } else if(k == p){
          ind <- seq(from = n.cumsum[[k-1]] + 1, to = j, by = 1)
        } else {
          ind <- seq(from = n.cumsum[[k-1]] + 1, to = n.cumsum[[k]], by = 1)
        }
        
        cond.ind[[k]] <- order(distance[ind])[seq(siz.each[k])]
        if(k != 1) {
          cond.ind[[k]] <- cond.ind[[k]] + as.numeric(n.cumsum[k-1])
        }
      }
      
      cond.sets[j, ] <- unlist(cond.ind)[order(rank(distance[unlist(cond.ind)], ties.method = "last"))]
      # cond.sets[j, ] <- unlist(cond.ind)[order(distance[unlist(cond.ind)])]
    }
  }
  
  Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs         <- rep(TRUE, nrow(locsord))
  U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx  <- list(locsord = locsord, obs = obs, ord = ord.full, ord.z = ord.full, ord.pred = 'general', U.prep = U.prep, cond.yz = 'false', ordering = 'maxmin', ordering.method = 'euclidean', conditioning = 'NN', conditioning.method = 'euclidean')
  
  return(vecchia.approx)
  # return(cond.sets)
}



#' @title The second baseline approximation for multivariate cases
#'
#' @param locs A matrix of locations
#' @param m A size of conditioning sets
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' locs1  <- matrix(runif(20), 10, 2)
#' locs2  <- matrix(runif(30), 15, 2)
#' locs   <- list(locs1 = locs1, locs2 = locs2)
#' m      <- 5
#' 
#' # baseline 1 = separately E-maxmin ordering + half-and-half E-NN conditioning
#' # baseline 2 = separately E-maxmin ordering + unified E-NN conditioning 
#' # baseline 3 = separately E-maxmin ordering + separately E-NN conditioning
#' # baseline 4 = separately E-maxmin ordering + unified C-NN conditioning
#' 
#' # Note 1 = For all baseline approximations, the first process is not relevant to the second one
#' # Note 2 = For the third baseline approximations, two processes are totally irrelevant
#' 
#' ### ordering
#' ord1 <- baseline_1_multivariate_specify(locs = locs, m = m)$ord
#' ord2 <- baseline_2_multivariate_specify(locs = locs, m = m)$ord
#' ord3 <- baseline_3_multivariate_specify(locs = locs, m = m)$ord
#' ord4 <- baseline_4_multivariate_specify(locs = locs, m = m, 
#'                                         covmodel = cov_bivariate_flexMatern, 
#'                                         sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                         nu.mat = matrix(0.5, 2, 2), 
#'                                         alpha.mat = matrix(1, 2, 2))$ord
#' 
#' all.equal(ord1, ord2) ; all.equal(ord1, ord3) ; all.equal(ord1, ord4)
#' 
#' ### conditioning
#' ind <- seq(from = m, to = 1, by = -1)
#' cond1 <- baseline_1_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond2 <- baseline_2_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond3 <- baseline_3_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond4 <- baseline_4_multivariate_specify(locs = locs, m = m, 
#'                                          covmodel = cov_bivariate_flexMatern, 
#'                                          sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                          nu.mat = matrix(0.5, 2, 2), 
#'                                          alpha.mat = matrix(1, 2, 2))$U.prep$revNNarray[, ind]
#' 
#' cond1 == cond2 # similar
#' cond1 == cond3 # different
#' cond3 == cond4 # similar
baseline_2_multivariate_specify <- function(locs, m)
{
  p         <- length(locs)
  n         <- unlist(lapply(locs, nrow))
  n.cumsum  <- cumsum(n)
  
  ord <- list()
  for(i in 1:p) {
    ord[[i]]    <- order_maxmin_euclidean( locs[[i]] )
    
    if(i > 1) {
      ord.full    <- c(ord.full, ord[[i]] + n.cumsum[i-1])
    } else {
      ord.full    <- ord[[i]]
    }
  }
  
  locs.full   <- do.call(rbind, locs)
  locsord     <- locs.full[ord.full, , drop = FALSE]
  cond.sets   <- GpGp::find_ordered_nn(locs = locsord, m = m)
  
  Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs         <- rep(TRUE, nrow(locs.full))
  U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx  <- list(locsord = locsord, obs = obs, ord = ord.full, ord.z = ord.full, ord.pred = 'general', U.prep = U.prep, cond.yz = 'false', ordering = 'maxmin', ordering.method = 'euclidean', conditioning = 'NN', conditioning.method = 'euclidean')
  
  return(vecchia.approx)
  # return(cond.sets)
}



#' @title The third baseline approximation for multivariate cases
#'
#' @param locs A matrix of locations
#' @param m A size of conditioning sets
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' locs1  <- matrix(runif(20), 10, 2)
#' locs2  <- matrix(runif(30), 15, 2)
#' locs   <- list(locs1 = locs1, locs2 = locs2)
#' m      <- 5
#' 
#' # baseline 1 = separately E-maxmin ordering + half-and-half E-NN conditioning
#' # baseline 2 = separately E-maxmin ordering + unified E-NN conditioning 
#' # baseline 3 = separately E-maxmin ordering + separately E-NN conditioning
#' # baseline 4 = separately E-maxmin ordering + unified C-NN conditioning
#' 
#' # Note 1 = For all baseline approximations, the first process is not relevant to the second one
#' # Note 2 = For the third baseline approximations, two processes are totally irrelevant
#' 
#' ### ordering
#' ord1 <- baseline_1_multivariate_specify(locs = locs, m = m)$ord
#' ord2 <- baseline_2_multivariate_specify(locs = locs, m = m)$ord
#' ord3 <- baseline_3_multivariate_specify(locs = locs, m = m)$ord
#' ord4 <- baseline_4_multivariate_specify(locs = locs, m = m, 
#'                                         covmodel = cov_bivariate_flexMatern, 
#'                                         sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                         nu.mat = matrix(0.5, 2, 2), 
#'                                         alpha.mat = matrix(1, 2, 2))$ord
#' 
#' all.equal(ord1, ord2) ; all.equal(ord1, ord3) ; all.equal(ord1, ord4)
#' 
#' ### conditioning
#' ind <- seq(from = m, to = 1, by = -1)
#' cond1 <- baseline_1_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond2 <- baseline_2_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond3 <- baseline_3_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond4 <- baseline_4_multivariate_specify(locs = locs, m = m, 
#'                                          covmodel = cov_bivariate_flexMatern, 
#'                                          sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                          nu.mat = matrix(0.5, 2, 2), 
#'                                          alpha.mat = matrix(1, 2, 2))$U.prep$revNNarray[, ind]
#' 
#' cond1 == cond2 # similar
#' cond1 == cond3 # different
#' cond3 == cond4 # similar
baseline_3_multivariate_specify <- function(locs, m)
{
  p         <- length(locs)
  n         <- unlist(lapply(locs, nrow))
  n.cumsum  <- cumsum(n)
  
  ord <- list()
  for(i in 1:p) {
    ord[[i]]    <- order_maxmin_euclidean( locs[[i]] )
    
    if(i > 1) {
      ord.full    <- c(ord.full, ord[[i]] + n.cumsum[i-1])
    } else {
      ord.full    <- ord[[i]]
    }
  }  
  
  locsord     <- do.call(rbind, locs)
  locsord     <- locsord[ord.full, , drop = FALSE]
  
  cond.each <- list()
  for(i in 1:p) {
    ind             <- seq(from = ifelse(i == 1, 1, n.cumsum[i-1] + 1), to = n.cumsum[i], by = 1)
    cond.each[[i]]  <- GpGp::find_ordered_nn(locs = locsord[ind, , drop = FALSE], m = m)
    
    if(i > 1) {
      cond.sets       <- rbind(cond.sets, cond.each[[i]] + n.cumsum[i-1])
    } else {
      cond.sets       <- cond.each[[i]]
    }
  }
  
  Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs         <- rep(TRUE, nrow(locsord))
  U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx  <- list(locsord = locsord, obs = obs, ord = ord.full, ord.z = ord.full, ord.pred = 'general', U.prep = U.prep, cond.yz = 'false', ordering = 'maxmin', ordering.method = 'euclidean', conditioning = 'NN', conditioning.method = 'euclidean')
  
  return(vecchia.approx)
  # return(cond.sets) 
}



#' @title The fourth baseline approximation for multivariate cases
#'
#' @param locs A matrix of locations
#' @param m A size of conditioning sets
#' @param abs.corr FALSE -> 1-rho or TRUE -> 1-|rho|
#' @param covmodel A covariance function
#' @param ... Covariance parameters
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' locs1  <- matrix(runif(20), 10, 2)
#' locs2  <- matrix(runif(30), 15, 2)
#' locs   <- list(locs1 = locs1, locs2 = locs2)
#' m      <- 5
#' 
#' # baseline 1 = separately E-maxmin ordering + half-and-half E-NN conditioning
#' # baseline 2 = separately E-maxmin ordering + unified E-NN conditioning 
#' # baseline 3 = separately E-maxmin ordering + separately E-NN conditioning
#' # baseline 4 = separately E-maxmin ordering + unified C-NN conditioning
#' 
#' # Note 1 = For all baseline approximations, the first process is not relevant to the second one
#' # Note 2 = For the third baseline approximations, two processes are totally irrelevant
#' 
#' ### ordering
#' ord1 <- baseline_1_multivariate_specify(locs = locs, m = m)$ord
#' ord2 <- baseline_2_multivariate_specify(locs = locs, m = m)$ord
#' ord3 <- baseline_3_multivariate_specify(locs = locs, m = m)$ord
#' ord4 <- baseline_4_multivariate_specify(locs = locs, m = m, 
#'                                         covmodel = cov_bivariate_flexMatern, 
#'                                         sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                         nu.mat = matrix(0.5, 2, 2), 
#'                                         alpha.mat = matrix(1, 2, 2))$ord
#' 
#' all.equal(ord1, ord2) ; all.equal(ord1, ord3) ; all.equal(ord1, ord4)
#' 
#' ### conditioning
#' ind <- seq(from = m, to = 1, by = -1)
#' cond1 <- baseline_1_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond2 <- baseline_2_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond3 <- baseline_3_multivariate_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond4 <- baseline_4_multivariate_specify(locs = locs, m = m, 
#'                                          covmodel = cov_bivariate_flexMatern, 
#'                                          sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                          nu.mat = matrix(0.5, 2, 2), 
#'                                          alpha.mat = matrix(1, 2, 2))$U.prep$revNNarray[, ind]
#' 
#' cond1 == cond2 # similar
#' cond1 == cond3 # different
#' cond3 == cond4 # similar
baseline_4_multivariate_specify <- function(locs, m, abs.corr = FALSE, covmodel, ...)
{
  p         <- length(locs)
  n         <- unlist(lapply(locs, nrow))
  n.cumsum  <- cumsum(n)
  
  ord <- list()
  locsord.each <- list()
  for(i in 1:p) {
    ord[[i]]    <- order_maxmin_euclidean( locs[[i]] )
    locsord.each[[i]] <- locs[[i]][ord[[i]], , drop = FALSE]
    
    if(i > 1) {
      ord.full    <- c(ord.full, ord[[i]] + n.cumsum[i-1])
    } else {
      ord.full    <- ord[[i]]
    }
  }
  
  locs.full   <- do.call(rbind, locs)
  locsord     <- locs.full[ord.full, , drop = FALSE]
  covmat      <- covmodel(locs = locsord.each, ...)
  
  if(abs.corr == FALSE) {
    cond.sets   <- conditioning_nn(m = m, d = 1 - covmat)
  } else if(abs.corr == TRUE) {
    cond.sets   <- conditioning_nn(m = m, d = 1 - abs(covmat))
  } else {
    stop("abs.corr must be logical!")
  }
  
  Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs         <- rep(TRUE, nrow(locs.full))
  U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx  <- list(locsord = locsord, obs = obs, ord = ord.full, ord.z = ord.full, ord.pred = 'general', U.prep = U.prep, cond.yz = 'false', ordering = 'maxmin', ordering.method = 'euclidean', conditioning = 'NN', conditioning.method = 'euclidean')
  
  return(vecchia.approx)
  # return(cond.sets)
}



#' @title The first baseline approximation for spatio-temporal cases
#'
#' @param locs A matrix of spatio-temporal locations
#' @param m A size of conditioning sets
#' @param coordinate integer or vector of integers in \code{1,...,d}
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' locs <- matrix(runif(30), nrow = 10, ncol = 3)
#' m <- 4
#' covparms <- c(1, 0.75, 50, 25)
#' 
#' # baseline 1 = time-based ordering + the m previous obs conditioning
#' # baseline 2 = time-based ordering + E-NN without 
#' #                                      distinction between time and space coordinates
#' # baseline 3 = time-based ordering + C-NN conditioning
#' 
#' ### ordering
#' ord1 <- baseline_1_spacetime_specify(locs = locs, m = m)$ord
#' ord2 <- baseline_2_spacetime_specify(locs = locs, m = m)$ord
#' ord3 <- baseline_3_spacetime_specify(locs = locs, m = m, 
#'                                      covmodel = cov_spacetime_expo, 
#'                                      covparms = covparms)$ord
#' 
#' all.equal(ord1, ord2) ; all.equal(ord1, ord3)
#' 
#' ### conditioning
#' ind <- seq(from = m, to = 1, by = -1)
#' cond1 <- baseline_1_spacetime_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond2 <- baseline_2_spacetime_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond3 <- baseline_3_spacetime_specify(locs = locs, m = m, 
#'                                       covmodel = cov_spacetime_expo, 
#'                                       covparms = covparms)$U.prep$revNNarray[, ind]
#' 
#' cond1
#' cond2 
#' cond3
baseline_1_spacetime_specify <- function(locs, m, coordinate = NULL)
{
  if( any(duplicated(locs)) ) stop("Time location is duplicated.")
  
  ord         <- order_time(locs = locs, coordinate = coordinate)
  locsord     <- locs[ord, , drop = FALSE]
  
  cond.sets   <- matrix(NA, nrow = nrow(locsord), ncol = m + 1)
  for(i in 1:nrow(locsord)) {
    ind                 <- seq(from = 1, to = min(i, m + 1), by = 1)
    cond.sets[i, ind]   <- seq(from = i, by = -1, length.out = length(ind)) 
  }
  
  Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs         <- rep(TRUE, nrow(locsord))
  U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx  <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred = 'general', U.prep = U.prep, cond.yz = 'false', ordering = 'coordinate', ordering.method = 'time', conditioning = 'NN', conditioning.method = 'time')
  
  return(vecchia.approx)
}


#' @title The second baseline approximation for spatio-temporal cases
#'
#' @param locs A matrix of spatio-temporal locations
#' @param m A size of conditioning sets 
#' @param coordinate integer or vector of integers in \code{1,...,d}
#' @param theta space distance + theta * time distance
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' locs <- matrix(runif(30), nrow = 10, ncol = 3)
#' m <- 4
#' covparms <- c(1, 0.75, 50, 25)
#' 
#' # baseline 1 = time-based ordering + the m previous obs conditioning
#' # baseline 2 = time-based ordering + E-NN without 
#' #                                      distinction between time and space coordinates
#' # baseline 3 = time-based ordering + C-NN conditioning
#' 
#' ### ordering
#' ord1 <- baseline_1_spacetime_specify(locs = locs, m = m)$ord
#' ord2 <- baseline_2_spacetime_specify(locs = locs, m = m)$ord
#' ord3 <- baseline_3_spacetime_specify(locs = locs, m = m, 
#'                                      covmodel = cov_spacetime_expo, 
#'                                      covparms = covparms)$ord
#' 
#' all.equal(ord1, ord2) ; all.equal(ord1, ord3)
#' 
#' ### conditioning
#' ind <- seq(from = m, to = 1, by = -1)
#' cond1 <- baseline_1_spacetime_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond2 <- baseline_2_spacetime_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond3 <- baseline_3_spacetime_specify(locs = locs, m = m, 
#'                                       covmodel = cov_spacetime_expo, 
#'                                       covparms = covparms)$U.prep$revNNarray[, ind]
#' 
#' cond1
#' cond2 
#' cond3
baseline_2_spacetime_specify <- function(locs, m, coordinate = NULL, theta = 1)
{
  if( any(duplicated(locs)) ) stop("Time location is duplicated.")
  
  ord         <- order_time(locs = locs, coordinate = coordinate)
  locsord     <- locs[ord, , drop = FALSE]
  
  d.spacetime <- fields::rdist(x1 = locs[, 1:2, drop = FALSE]) + theta * fields::rdist(x1 = locs[, 3, drop = FALSE])
  cond.sets   <- conditioning_nn(m = m, d = d.spacetime)
  
  Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs         <- rep(TRUE, nrow(locsord))
  U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx  <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred = 'general', U.prep = U.prep, cond.yz = 'false', ordering = 'coordinate', ordering.method = 'time', conditioning = 'NN', conditioning.method = 'spacetime')
  
  return(vecchia.approx)
}



#' @title The third baseline approximation for spatio-temporal cases
#'
#' @param locs A matrix of spatio-temporal locations
#' @param m A size of conditioning sets  
#' @param coordinate integer or vector of integers in \code{1,...,d}
#' @param covmodel covariance function
#' @param covparms A numeric vector of covariance parameters
#'
#' @return list
#' 
#' @export
#'
#' @examples
#' locs <- matrix(runif(30), nrow = 10, ncol = 3)
#' m <- 4
#' covparms <- c(1, 0.75, 50, 25)
#' 
#' # baseline 1 = time-based ordering + the m previous obs conditioning
#' # baseline 2 = time-based ordering + E-NN without 
#' #                                      distinction between time and space coordinates
#' # baseline 3 = time-based ordering + C-NN conditioning
#' 
#' ### ordering
#' ord1 <- baseline_1_spacetime_specify(locs = locs, m = m)$ord
#' ord2 <- baseline_2_spacetime_specify(locs = locs, m = m)$ord
#' ord3 <- baseline_3_spacetime_specify(locs = locs, m = m, 
#'                                      covmodel = cov_spacetime_expo, 
#'                                      covparms = covparms)$ord
#' 
#' all.equal(ord1, ord2) ; all.equal(ord1, ord3)
#' 
#' ### conditioning
#' ind <- seq(from = m, to = 1, by = -1)
#' cond1 <- baseline_1_spacetime_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond2 <- baseline_2_spacetime_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
#' cond3 <- baseline_3_spacetime_specify(locs = locs, m = m, 
#'                                       covmodel = cov_spacetime_expo, 
#'                                       covparms = covparms)$U.prep$revNNarray[, ind]
#' 
#' cond1
#' cond2 
#' cond3
baseline_3_spacetime_specify <- function(locs, m, coordinate = NULL, covmodel, covparms)
{
  ord         <- order_time(locs = locs, coordinate = coordinate)
  locsord     <- locs[ord, , drop = FALSE]
  
  covmat      <- covmodel(locs = locs, covparms = covparms) / covparms[1]
  covmatord   <- covmat[ord, ord]
  
  cond.sets   <- conditioning_nn(m = m, d = 1 - covmat)
  
  Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs         <- rep(TRUE, nrow(locsord))
  U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx  <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred = 'general', U.prep = U.prep, cond.yz = 'false', ordering = 'coordinate', ordering.method = 'time', conditioning = 'NN', conditioning.method = 'spacetime')
  
  return(vecchia.approx)
}



#' @title Straightforward time-based ordering of spatio-temporal locations
#'
#' @param locs A matrix of spatio-temporal locations
#' @param coordinate integer or vector of integers in \code{1,...,d}. At \code{NULL} by default
#'
#' @return A vector of indices giving the time-based ordering. If two observations are measured at the same time, then they are ordered with respect to their spatial locations.
#' 
#' @export
#'
#' @examples
#' # all random
#' locs <- matrix(runif(30), 10, 3)
#' locs[order_coordinate(locs, coordinate = 3), ]
#' locs[order_time(locs, coordinate = NULL), ] # different
#' 
#' # repeated measurement
#' locs <- matrix(runif(20), 10, 2)
#' locs <- cbind(locs, sample(1:4, 10, replace = TRUE))
#' locs[order_coordinate(locs, coordinate = 3), ]
#' locs[order_time(locs, coordinate = NULL), ] # different
#' 
#' # In detail,
#' locs <- as.matrix(expand.grid(seq(4), seq(4), seq(4)))[sample(4^3), ]
#' all.equal(do.call(order, data.frame(locs)[, 3:1]), order(locs[, 3], locs[, 2], locs[, 1])) # same
order_time <- function(locs, coordinate = NULL)
{
  d         <- ncol(locs) - 1
  if(d < 1) stop("The number of columns of the locs must be greater than 1.")
  
  if(is.null(coordinate)) {
    coordinate <- seq(from = d+1, to = 1, by = -1)
  } else {
    coordinate <- c(d+1, coordinate)
    coordinate <- unique(coordinate)
  }
  
  ord <- do.call(order, data.frame(locs)[, coordinate, drop = FALSE]) # If d = 2 and coordinate = c(3, 2, 1), then ord = order(locs[, 3], locs[, 2], locs[, 1])
  return(ord)
}


#' @title Straightforward time-based ordering of spatio-temporal locations (old version)
#'
#' @param locs A matrix of spatio-temporal locations
#' @param coordinate integer or vector of integers in \code{1,...,d}. At \code{NULL} by default
#'
#' @return A vector of indices giving the time-based ordering. If two observations are measured at the same time, then they are ordered with respect to their spatial locations.
#' 
#' @export
#'
#' @examples
#' locs <- generate_gp_spacetime(nsim = 1, n = 100, d = 2, t.len = 3, 
#'                               method.locs = 'all.random', 
#'                               covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
#' locs <- locs$sim$sim1$locs
#' 
#' all.equal(order_time(locs), order_time_old(locs)) # same
#' 
#' locs <- generate_gp_spacetime(nsim = 1, n = 100, d = 2, t.len = 3, 
#'                               method.locs = 'space.random.time.grid', 
#'                               covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
#' locs <- locs$sim$sim1$locs
#' 
#' all.equal(order_time(locs), order_time_old(locs)) # same
#' 
#' locs <- generate_gp_spacetime(nsim = 1, n = 100, d = 2, t.len = 3, 
#'                               method.locs = 'all.grid', 
#'                               covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
#' locs <- locs$sim$sim1$locs
#' 
#' all.equal(order_time(locs), order_time_old(locs)) # same
#' 
#' locs <- generate_gp_spacetime(nsim = 1, n = 100, d = 2, t.len = 3, 
#'                               method.locs = 'satellite', 
#'                               covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
#' locs <- locs$sim$sim1$locs
#' 
#' all.equal(order_time(locs), order_time_old(locs)) # same
#' 
#' # However...
#' locs <- as.matrix(expand.grid(seq(4), seq(4)))[sample(4^2), ]
#' all.equal(order_time(locs), order_time_old(locs)) # same
#' 
#' locs <- as.matrix(expand.grid(seq(4), seq(4), seq(4)))[sample(4^3), ]
#' all.equal(order_time(locs), order_time_old(locs)) # different!
order_time_old <- function(locs, coordinate = NULL)
{
  d         <- ncol(locs) - 1
  if(d < 1) stop("The number of columns of the locs must be greater than 1.")
  
  if(is.null(coordinate)) coordinate <- d
  
  t.coord <- locs[, d+1]
  s.coord <- rowSums(locs[, coordinate, drop = FALSE])
  
  return( order(t.coord, s.coord) )
}
