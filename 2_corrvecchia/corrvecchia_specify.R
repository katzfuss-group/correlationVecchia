####################################################################
####
####  Author: Myeongjong (MJ) Kang (kmj.stat@gmail.com)
####
####  Description: 
####
####################################################################

####################################################################
#### load required packages
####################################################################

library(GPvecchia)


####################################################################
#### corrvecchia_specify()
####################################################################

#' @param dat: the observed data
#' @param locs: nxd matrix of observed locs
#' @param m: Number of nearby points to condition on
#'
#' @param ordering: options are 'coord' or 'maxmin'
#' @param which.coord: 
#' @param cond.yz: options are 'y', 'z', 'SGV', 'SGVT', 'RVP', 'LK', and 'zy'
#'
#' @param conditioning: conditioning on 'NN' (nearest neighbor) or 'firstm' (fixed set for low rank) or 'mra'
#' @param covmodel: covariance model, 'matern' by default
#' @param covparms: covariance parameters as a vector
#' @param nugget: either a single (constant) nugget or a vector of nugget terms for the observations
#' @param covparms.ini: initial values of covariance parameters
#' @param nugget.ini: initial value of single nugget
#'
#' @param locs.pred: nxd matrix of locations at which to make predictions
#' @param ordering.pred: options are 'obspred' or 'general'
#' @param pred.cond: prediction conditioning, options are 'general' or 'independent'
#'
#' @param mra.options: settings for number of levels and neighbors per level
#' @param verbose: provide more detail when using MRA calculations. Default is false
#'
#' @return An object that specifies the vecchia approximation for later use in likelihood evaluation or prediction
corrvecchia_specify <- function(dat, locs, m=-1, ordering, which.coord=NULL, cond.yz, conditioning, covmodel='matern', covparms=NULL, nugget=NULL, covparms.ini=NULL, nugget.ini=NULL, locs.pred, ordering.pred, pred.cond, mra.options=NULL, verbose=FALSE)
{
  
  ##### Check whether the input "locs" is proper or not #####
  if(!is.matrix(locs)) {
    warning("Locations must be in matrix form")
    return(NA)
  }
  
  ##### Check whether the input "m" is proper or not #####
  if(m==-1 || is.null(m)) {
    if(conditioning=='mra' && !is.null(mra.options) &&  !is.null(mra.options$J) && !is.null(mra.options$r) && !is.null(mra.options$J)) {
      warning("m not defined; using MRA parameters")
    } else if(is.null(mra.options$r)) {
      stop("neither m nor r defined!")
    }
  }
  
  ##### Default values: spatial.dim, n, which.coord #####
  spatial.dim   <- ncol(locs)
  n             <- nrow(locs)
  
  if(is.null(which.coord)) which.coord <- 1:spatial.dim
  
  ##### Check whether the input "locs.preds" intersects the input "locs" or not #####
  if(!missing(locs.pred)){
    locs.all <- rbind(locs, locs.pred)
    if(anyDuplicated(locs.all)>0) stop("Prediction locations contain observed location(s), remove redundancies.")
  }
  
  ##### Check whether the input "m" is smaller than the input "n" or not #####
  if(m>n){
    warning("Conditioning set size m chosen to be larger than n. Changing to m=n-1")
    m=n-1
  }
  
  ##### Case 0: the fully independent case with no conditioning sets (the input "m" equals 0) #####
  if(m==0){
    if(!missing(locs.pred)) cat("Attempting to make predictions with m=0. Prediction ignored")
    
    ord       <- 1:n
    ord.z     <- ord
    locsord   <- locs[ord, , drop=FALSE]
    NNarray   <- matrix(cbind(ord, rep(NA, nrow(locs))), ncol=2)
    Cond      <- matrix(NA, n, 2); Cond[ ,1] <- T
    obs       <- rep(TRUE, n)
    
    ### determine the sparsity structure of U
    U.prep    <- U_sparsity(locsord, NNarray, obs, Cond)
    
    ### object that specifies the vecchia approximation
    vecchia.approx <- list(locsord=locsord, obs=obs, ord=ord, ord.z=ord.z, ord.pred='general', U.prep=U.prep,cond.yz='false',conditioning='NN')
    return(vecchia.approx)
  }
  
  ##### Subsume firstm into mra #####
  if(!missing(conditioning) && conditioning == 'firstm'){
    conditioning  <- 'mra'
    mra.options   <- list(r=c(m,1))
  }
  
  ##### Default options: ordering, pred.cond, conditioning, cond.yz #####
  if(missing(ordering)){
    if(spatial.dim==1) ordering <- 'coord' else ordering <- 'maxmin'
  }
  if(missing(pred.cond)) pred.cond <- 'general'
  if(missing(conditioning)) conditioning <- 'NN'
  if(conditioning %in% c('mra', 'firstm')) {
    if(ordering!='maxmin') warning("ordering for the selected conditioning scheme changed to required 'maxmin'")
    ordering <- 'maxmin'
  }
  if(missing(cond.yz)) {
    if (conditioning %in% c('mra', 'firstm')) { 
      cond.yz <- 'y'
    } else if(missing(locs.pred) | spatial.dim==1 ) { 
      cond.yz <- 'SGV'
    } else {
      cond.yz <- 'zy'
    }
  }
  
  ##### corrvecchia #####
  if(is.null(covparms) & is.null(nugget) & !is.null(covparms.ini) & !is.null(nugget.ini)) {
    # Start with estimating covparms and nugget
    vecchia.approx <- corrvecchia_specify_unknownparams(dat, locs, m, ordering, which.coord, cond.yz, conditioning, covmodel, covparms.ini, nugget.ini, locs.pred, ordering.pred, pred.cond, mra.options, verbose)
  } else if(!is.null(covparms) & !is.null(nugget) & is.null(covparms.ini) & is.null(nugget.ini)) {
    # Start with pre-specified covparms and nugget
    vecchia.approx <- corrvecchia_specify_knownparams(locs, m, ordering, which.coord, cond.yz, conditioning, covmodel, covparms, nugget, locs.pred, ordering.pred, pred.cond, mra.options, verbose)
  } else {
    stop("input error!!!")
  }
  
  return(vecchia.approx)
}

##### not clear!!!! #####
corrvecchia_specify_unknownparams <- function(dat, locs, m, ordering, which.coord=NULL, cond.yz, conditioning, covmodel, covparms.ini, nugget.ini, locs.pred, ordering.pred, pred.cond, mra.options=NULL, verbose=FALSE)
{
  # iteratively...? # Estimate covparams and nugget
  # iteratively...? # Order locs and z
  # iteratively...? # Obtain conditioning sets
  
  # Determine the sparsity structure of U
  # Return an object that specifies the correlation-based Vecchia approximation
}


corrvecchia_specify_knownparams <- function(locs, m, ordering, which.coord=NULL, cond.yz, conditioning, covmodel, covparms, nugget, locs.pred, ordering.pred, pred.cond, mra.options=NULL, verbose=FALSE)
{
  # Order locs and z
  # Obtain conditioning sets
  # Determine the sparsity structure of U
  # Return an object that specifies the correlation-based Vecchia approximation
}