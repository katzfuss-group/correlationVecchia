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
  
  ##### Order locs and z #####
  if(missing(locs.pred)){  # no prediction
    
    if(ordering=='coord') { 
      ord <- order_coordinate(locs, which.coord)
    } else if(ordering=='maxmin') { 
      ord <- order_maxmin_exact(locs)
    } else if(ordering=='outsidein') {
      ord <- order_outsidein(locs)
    } else if(ordering=='none') {
      ord <- seq(n)
    }
    
    ord.z           <- ord
    locsord         <- locs[ord, , drop=FALSE]
    obs             <- rep(TRUE, n)
    ordering.pred   <- 'general'
    
  } else {    # prediction is desired
    
    n.p                 <- nrow(locs.pred)
    locs.all            <- rbind(locs, locs.pred)
    observed.obspred    <- c(rep(TRUE, n), rep(FALSE, n.p))
    
    if(missing(ordering.pred)) {
      if(spatial.dim==1 & ordering=='coord') ordering.pred <- 'general' else ordering.pred <- 'obspred'
    }
      
    if(ordering.pred=='general'){
      if(ordering=='coord') ord <- order_coordinate(locs.all) else ord <- order_maxmin_exact(locs.all)
      ord.obs <- ord[ord<=n]
    } else {
      if(ordering=='coord') {
        ord.obs   <- order_coordinate(locs, which.coord)
        ord.pred  <- order_coordinate(locs.pred, which.coord)
      } else if(ordering=='none') {
        ord.obs   <- seq(n)
        ord.pred  <- seq(n.p)
      } else {
        temp      <- order_maxmin_exact_obs_pred(locs, locs.pred)
        ord.obs   <- temp$ord
        ord.pred  <- temp$ord_pred
      }
      ord <- c(ord.obs, ord.pred+n)
    }
    
    ord.z     <- ord.obs
    locsord   <- locs.all[ord, , drop=FALSE]
    obs       <- observed.obspred[ord]
  }
  
  ##### Obtain conditioning sets #####
  if( conditioning == 'mra' ) {
    
    NNarray <- findOrderedNN_mra(locsord, mra.options, m, verbose)
    if(!methods::hasArg(m)) m <- ncol(NNarray) - 1
    
  } else if( conditioning %in% c('firstm', 'NN')) {
    
    if(spatial.dim == 1) {
      NNarray <- findOrderedNN_kdtree2(locsord, m)
    } else {
      NNarray <- GpGp::find_ordered_nn(locsord, m)
    }
    
    if(conditioning == 'firstm'){
      first_m   <- NNarray[m+1,2:(m+1)]
      n.all     <- nrow(NNarray)
      
      # if m=n-1, nothing to replace
      if (m < n.all-1) NNarray[(m+2):n.all, 2:(m+1)] <- matrix(rep(first_m, n.all-m-1), byrow = TRUE, ncol = m)
    }
    
  } else {
    stop(paste0("conditioning='", conditioning, "' not defined"))
  }
  
  if(!missing(locs.pred) & pred.cond == 'independent') {
    if(ordering.pred == 'obspred') {
      
      NNarray.pred <- array(dim=c(n.p, m+1))
      for(j in 1:n.p){
        dists             <-fields::rdist(locsord[n+j, , drop=FALSE], locsord[1:n, , drop=FALSE])
        m.nearest.obs     <- sort(order(dists)[1:m], decreasing=TRUE)
        NNarray.pred[j,]  <- c(n+j, m.nearest.obs)
      }
      
      NNarray[n+(1:n.p),] <- NNarray.pred
      
    } else {
      print('indep. conditioning currently only implemented for obspred ordering')
    }
  }
  
  ##### Conditioning on y or z? #####
  if(cond.yz == 'SGV'){
    Cond    <- whichCondOnLatent(NNarray, firstind.pred=n+1)
  } else if(cond.yz == 'SGVT'){
    Cond    <- rbind(whichCondOnLatent(NNarray[1:n,]),matrix(TRUE,nrow=n.p,ncol=m+1))
  } else if(cond.yz == 'y'){
    Cond    <- matrix(NA, nrow(NNarray), ncol(NNarray)); Cond[!is.na(NNarray)] <- TRUE
    #Cond=matrix(NA,nrow(NNarray),m+1); Cond[!is.na(NNarray)]=FALSE; Cond[,1]=TRUE
  } else if(cond.yz == 'z'){
    Cond    <- matrix(NA, nrow(NNarray), m+1); Cond[!is.na(NNarray)] <- FALSE; Cond[,1] <- TRUE
  } else if(cond.yz %in% c('RVP','LK','zy')){
    ### "trick" code into response-latent ('zy') ordering
    
    ## reset variables
    obs       <- c(rep(TRUE, n), rep(FALSE, nrow(locsord)))
    locsord   <- rbind(locsord[1:n, , drop=FALSE], locsord)
    
    ## specify neighbors
    NNs <- FNN::get.knn(locsord[1:n, , drop=FALSE], m-1)$nn.index
    if(cond.yz %in% c('RVP','zy')){
      prev        <- (NNs<matrix(rep(1:n, m-1), nrow=n))
      NNs[prev]   <- NNs[prev] + n # condition on latent y.obs if possible
    }
    
    ## create NN array
    NNarray.z   <- cbind(1:n, matrix(nrow=n, ncol=m))
    NNarray.y   <- cbind((1:n)+n, 1:n, NNs)
    if(missing(locs.pred)){
      NNarray.yp      <- matrix(nrow=0, ncol=m+1)
      ordering.pred   <- 'obspred'
    } else {
      if(ordering.pred!='obspred') print('Warning: ZY only implemented for obspred ordering')
      if(cond.yz=='zy'){
        NNarray.yp                <- NNarray[n+(1:n.p), ] + n
      } else {
        NNarray.yp                <- NNarray[n+(1:n.p), ]
        NNarray.yp[NNarray.yp>n]  <- NNarray.yp[NNarray.yp>n] + n
      }
    }
    NNarray <- rbind(NNarray.z, NNarray.y, NNarray.yp)
    
    ## conditioning
    Cond      <- (NNarray>n); Cond[,1] <- TRUE
    cond.yz   <- 'zy'
    
  } else {
    stop(paste0("cond.yz='", cond.yz, "' not defined"))
  }
  
  ##### Determine the sparsity structure of U #####
  U.prep <- U_sparsity(locsord, NNarray, obs, Cond)
  
  ##### Return an object that specifies the correlation-based Vecchia approximation #####
  vecchia.approx <- list(locsord=locsord, obs=obs, ord=ord, ord.z=ord.z, ord.pred=ordering.pred, U.prep=U.prep, cond.yz=cond.yz, conditioning=conditioning)
  return(vecchia.approx)
}