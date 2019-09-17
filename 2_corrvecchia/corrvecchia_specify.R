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
  if(is.null(covparms) & is.null(nugget) & !is.null(covparms.ini) & !is.null(nugget.ini)) {
    # Start with estimating covparms and nugget
    corrvecchia_specify_unknownparams(dat, locs, m, ordering, which.coord, cond.yz, conditioning, covmodel, covparms.ini, nugget.ini, locs.pred, ordering.pred, pred.cond, mra.options, verbose)
  } else if(!is.null(covparms) & !is.null(nugget) & is.null(covparms.ini) & is.null(nugget.ini)) {
    # Start with pre-specified covparms and nugget
    corrvecchia_specify_knownparams(locs, m, ordering, which.coord, cond.yz, conditioning, covmodel, covparms, nugget, locs.pred, ordering.pred, pred.cond, mra.options, verbose)
  } else {
    stop("input error!!!")
  }
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