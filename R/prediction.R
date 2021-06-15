####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to obtain prediction.
###
###   Contents: 
###     prediction_spacetime_baseline
###     prediction_spacetime_cvecchia
###
####################################################################################

#' @title Prediction of baseline approximations for spacetime GPs
#'
#' @param approx b1, b2, or b3
#' @param z A realization of the process of interest
#' @param locs A numeric matrix of observed locations
#' @param locs.pred A numeric matrix of unobserved locations
#' @param m The size of conditioning sets
#' @param coordinate coordinate at NULL by default
#' @param covmodel Covariance model
#' @param covparms A numeric vector of covariance parameters
#' @param nuggets nuggets
#' @param predcond.method It must be general for now 
#' @param var.exact should prediction variances be computed exactly, or is a (faster) approximation acceptable
#' @param return.values either 'mean' only, 'meanvar', 'meanmat', or 'all'
#'
#' @return list of vecchia.approx and vecchia.predict
#' 
#' @export
#'
#' @examples
#' 1 + 1
prediction_spacetime_baseline <- function(approx, z, locs, locs.pred, m, coordinate = NULL, covmodel, covparms, nuggets, predcond.method = "general", var.exact, return.values = "all")
{
  ### initialization
  
  ### main
  if( approx == 1 | approx %in% c("b1", "T-MM + T-NN") ) {
    
    ## ordering
    ord.obsv          <- order_time(locs = locs, coordinate = coordinate)
    ord.pred          <- order_time(locs = locs.pred, coordinate = coordinate)
    ord               <- c(ord.obsv, ord.pred + nrow(locs))
    
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    ## conditioning
    cond.sets   <- matrix(NA, nrow = nrow(locsord), ncol = m + 1)
    for(i in 1:nrow(locsord)) {
      ind                 <- seq(from = 1, to = min(i, m + 1), by = 1)
      cond.sets[i, ind]   <- seq(from = i, by = -1, length.out = length(ind)) 
    }
    
    Cond        <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep      <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    ## prediction
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    # cormat            <- .correlation(locs = locsord, covmodel = covmodel, covparms = covparms, abs.corr = FALSE)
    
    vecchia.predict   <- GPvecchia::vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms, nuggets = nuggets, var.exact = var.exact, covmodel = covmodel(locs = locsord, covparms = covparms), return.values = return.values)
    
  } else if( approx == 2 | approx %in% c("b2", "T-MM + E-NN") ) {
    
    ## ordering
    ord.obsv          <- order_time(locs = locs, coordinate = coordinate)
    ord.pred          <- order_time(locs = locs.pred, coordinate = coordinate)
    ord               <- c(ord.obsv, ord.pred + nrow(locs))
    
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    ## conditioning
    d.spacetime       <- fields::rdist(x1 = locsord[, 1:2, drop = FALSE]) + 1 * fields::rdist(x1 = locsord[, 3, drop = FALSE])
    
    cond.sets         <- conditioning_m_Rcpp(m = m, d = d.spacetime) + 1
    
    Cond              <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep            <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    ## prediction
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    # cormat            <- .correlation(locs = locsord, covmodel = covmodel, covparms = covparms, abs.corr = FALSE)
    
    vecchia.predict   <- GPvecchia::vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms, nuggets = nuggets, var.exact = var.exact, covmodel = covmodel(locs = locsord, covparms = covparms), return.values = return.values)
    
  } else if( approx == 3 | approx %in% c("b3", "T-MM + C-NN") ) {
    
    ## ordering
    ord.obsv          <- order_time(locs = locs, coordinate = coordinate)
    ord.pred          <- order_time(locs = locs.pred, coordinate = coordinate)
    ord               <- c(ord.obsv, ord.pred + nrow(locs))
    
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    ## conditioning
    cormat            <- .correlation(locs = locsord, covmodel = covmodel, covparms = covparms, abs.corr = FALSE)
    cond.sets         <- conditioning_m_Rcpp(m = m, d = 1 - cormat) + 1
    Cond              <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep            <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    ## prediction
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    vecchia.predict   <- GPvecchia::vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms, nuggets = nuggets, var.exact = var.exact, covmodel = covmodel(locs = locsord, covparms = covparms), return.values = return.values)
    
  } else {
    
    stop("Check the argument approx! It must be b1, b2, or b3.")
  }
  
  ### return
  return(list(approx = vecchia.approx, predict = vecchia.predict))
}

#' @title Prediction of Vecchia approximations for spacetime GPs
#'
#' @param method euc or cor
#' @param z A realization of the process of interest
#' @param locs A numeric matrix of observed locations
#' @param locs.pred A numeric matrix of unobserved locations
#' @param m The size of conditioning sets 
#' @param initial.pt initial.pt
#' @param covmodel Covariance model
#' @param covparms A numeric vector of covariance parameters
#' @param nuggets nuggets
#' @param predcond.method It must be general for now 
#' @param var.exact should prediction variances be computed exactly, or is a (faster) approximation acceptable
#' @param return.values either 'mean' only, 'meanvar', 'meanmat', or 'all'
#'
#' @return list of vecchia.approx and vecchia.predict
#' 
#' @export
#'
#' @examples
#' 1 + 1
prediction_spacetime_cvecchia <- function(method, z, locs, locs.pred, m, initial.pt = NULL, covmodel, covparms, nuggets, predcond.method = "general", var.exact, return.values = "all")
{
  ### initialization
  if(is.function(covmodel)) {
    
    covname <- deparse(substitute(covmodel))    
    
  } else if(is.character(covmodel)) {
    
    covname <- covmodel
    
  } else {
    
    stop("Please use a function or string to specify the argument covmodel!")
  }

  if( covname %in% c("cov_expo_spacetime", "cov_expo_spacetime_cpp") ) {
    
    covmodel  <- cov_expo_spacetime
    covftn    <- "cov_expo_spacetime"
    
  } else if( covname %in% c("cov_matern_spacetime", "cov_matern_spacetime_cpp") ) {
    
    covmodel  <- cov_matern_spacetime
    covftn    <- "cov_matern_spacetime"
    
  } else if( covname %in% c("cov_matern_scaledim", "cov_matern_scaledim_cpp") ) {
    
    covmodel  <- cov_matern_scaledim
    covftn    <- "cov_matern_scaledim"
    
  } else {
    
    stop("This function only works for cov_expo_spacetime, cov_matern_spacetime, and cov_expo_scaledim for now.")
  }
  
  if( method %in% c("cor", "correlation") ) {
    
    if( is.null(initial.pt) ){
      cormat            <- .correlation(locs = locs, covmodel = covmodel, covparms = covparms, abs.corr = FALSE)
      initial.pt        <-  which.max(rowSums(cormat))
    }
    
  } else if( method %in% c("euc", "euclidean") ) {
    
    if(is.null(initial.pt)) {
      cen               <- t(as.matrix(colMeans(locs)))
      initial.pt        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = nrow(locs), ncol = ncol(locs), byrow = T))^2))
    }
    
  } else {
    
    stop("Check the argument method! It must be euc or cor.")
  }
  
  ### main
  if( method %in% c("cor", "correlation") ) {
    
    ## ordering 
    ord               <- predSortSparse_Rcpp(locs, locs.pred, 2.01, initial.pt-1, "correlation", covftn, covparms)$P + 1
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    ## conditioning
    cormat            <- .correlation(locs = locsord, covmodel = covmodel, covparms = covparms, abs.corr = FALSE)
    cond.sets         <- conditioning_m_Rcpp(m = m, d = 1-cormat) + 1
    Cond              <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep            <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    ## prediction
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    vecchia.predict   <- GPvecchia::vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms, nuggets = nuggets, var.exact = var.exact, covmodel = covmodel(locs = locsord, covparms = covparms), return.values = return.values)
    
  } else if( method %in% c("euc", "euclidean") ) {
    
    ## ordering 
    ord               <- predSortSparse_Rcpp(locs, locs.pred, 2.01, initial.pt-1, "euclidean", covftn, covparms)$P + 1
    locsord           <- rbind(locs, locs.pred)[ord, , drop = FALSE]
    obs               <- c(rep(TRUE, nrow(locs)), rep(FALSE, nrow(locs.pred)))
    
    ord.z             <- ord[seq(nrow(locs))]
    ord.pred          <- "obspred"
    cond.yz           <- "y"
    
    ## conditioning
    cond.sets         <- GpGp::find_ordered_nn(locs = locsord, m = m)
    Cond              <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    U.prep            <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
    
    ## prediction
    vecchia.approx    <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord.z, ord.pred = predcond.method, U.prep = U.prep, cond.yz = cond.yz, ordering = "maxmin", conditioning = "NN", ic0 = FALSE)
    
    vecchia.predict   <- GPvecchia::vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms, nuggets = nuggets, var.exact = var.exact, covmodel = covmodel(locs = locsord, covparms = covparms), return.values = return.values)
    
  } else {
    
    stop("Check the argument method! It must be euc or cor.")
  }
  
  ### return
  return(list(approx = vecchia.approx, predict = vecchia.predict))
}
