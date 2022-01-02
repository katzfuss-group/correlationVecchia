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

#' @title prediction using scaled Vecchia, using output from fit_scaled()
#'
#' @param fit object returned from fit_scaled()
#' @param locs_pred n.p x d matrix of test/prediction inputs/locations
#' @param m conditioning-set size (larger is more accurate but slower)
#' @param joint Joint predictions (does not return variances) or separate/independent predictions (does not produce joint samples)
#' @param nsims desired number of samples from predictive distributions. if \code{n.sims=0}, posterior mean is returned.
#' @param predvar return prediction variances? (only if \code{joint=FALSE})
#' @param X_pred n.p x p trend matrix at locs_pred (if missing, will be generated based on fit object)
#' @param scale scaling of inputs for ordering and conditioning. 'parms': by parameter estimates. 'ranges': to the unit interval. 'none': no scaling
#' @param tol numerical tolerence
#'
#' @return Vector of length n.p (\code{n.sims=0}, \code{predvar=FALSE}) or list with entries \code{means} and/or \code{vars} and/or \code{samples}
#'
#' @export
#'
#' @examples
#' library(GpGp)
#'
#' inputs=matrix(runif(200),ncol=2)
#' y=sin(rowSums(inputs*5))
#' inputs.test=matrix(runif(100),ncol=2)
#' fit=fit_scaled(y,inputs)
#' preds=predictions_scaled(fit,inputs.test)
#' plot(rowSums(inputs.test),preds)
predictions_scaled <- function(fit,locs_pred,m=100,joint=TRUE,nsims=0,
                               predvar=FALSE,X_pred,scale='parms',tol=1e-8){

  y_obs = fit$y
  locs_obs = fit$locs
  X_obs = fit$X
  beta = fit$betahat
  covparms = fit$covparms
  covfun_name = fit$covfun_name
  n_obs <- nrow(locs_obs)
  n_pred <- nrow(locs_pred)
  if(is.null(fit$vcf)) vcf=1 else vcf=fit$vcf

  # ## add nugget for numerical stability
  # if(covparms[length(covparms)]==0)
  #   covparms[length(covparms)]=covparms[1]*1e-12

  # specify trend if missing
  if(missing(X_pred)){
    if(fit$trend=='zero'){
      X_pred=as.matrix(rep(0,n_pred))
    } else if(fit$trend=='intercept'){
      X_pred=as.matrix(rep(1,n_pred))
    } else if(fit$trend=='linear'){
      X_pred=cbind(rep(1,n_pred),locs_pred)
    } else stop('X_pred must be specified')
  }

  # specify how to scale input dimensions
  if(scale=='parms'){ scales=1/covparms[1+(1:ncol(locs_obs))]
  } else if(scale=='ranges'){ scales=1/apply(locs_obs,2,function(x) diff(range(x)))
  } else stop(paste0('invalid argument scale=',scale))


  ###
  if(joint){  # joint predictions

    # get orderings
    temp=order_maxmin_pred(t(t(locs_obs)*scales),t(t(locs_pred)*scales))
    ord1=temp$ord
    ord2=temp$ord_pred

    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    y  <- y_obs[ord1] - Xord_obs %*% beta

    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )

    # get nearest neighbor array
    sm = if (n_pred<1e5) 2 else 1.5
    NNarray_all <- find_ordered_nn_pred(t(t(locs_all)*scales),m,
                                        fix.first=n_obs,searchmult=sm)
    NNarray_pred=NNarray_all[-(1:n_obs),-1]

    means=numeric(length=n_pred)
    if(nsims>0) samples=array(dim=c(n_pred,nsims))

    # make predictions sequentially
    for(i in 1:n_pred){

      # NN conditioning sets
      NN=sort(NNarray_pred[i,])
      NN_obs=NN[NN<=n_obs]
      NN_pred=NN[NN>n_obs]-n_obs

      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(locs_all[c(NN,i+n_obs),]))
      cl=t(chol(K + tol * diag(nrow(K))))

      # prediction
      y.NN=y[NN_obs]
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,means[NN_pred]))
      if(nsims>0){ # conditional simulation
        pred.var=cl[m+1,m+1]^2*vcf
        for(s in 1:nsims){
          pm=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,samples[NN_pred,s]))
          samples[i,s]=stats::rnorm(1,pm,sqrt(pred.var))
        }
      }

    }

    # add (prior) mean and return to original ordering
    means[ord2] = means + c(Xord_pred %*% beta)
    if(nsims==0){
      preds=means
    } else {
      samples[ord2,] = samples + c(Xord_pred %*% beta)
      preds=list(means=means,samples=samples)
    }

  } else {  # separate predictions

    if(nsims>0) stop('cannot produce joint samples when joint=FALSE')

    y  = y_obs - X_obs %*% beta

    # find the NNs
    m=min(m,nrow(locs_obs))
    NNarray=FNN::get.knnx(t(t(locs_obs)*scales),
                          t(t(locs_pred)*scales),m)$nn.index


    means=vars=numeric(length=n_pred)
    for(i in 1:n_pred){

      # NN conditioning sets
      NN=NNarray[i,]

      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(rbind(locs_obs[NN,],locs_pred[i,])))
      cl=t(chol(K + tol * diag(nrow(K))))

      # prediction
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],y[NN])
      vars[i]=cl[m+1,m+1]^2*vcf

    }
    means=means+c(X_pred %*% beta)

    if(predvar==FALSE){
      preds=means
    } else {
      preds=list(means=means,vars=vars)
    }

  }

  return(preds)
}

#' @title Maximum-minimum (maxmin) ordering of observed and unobserved locations in terms of the Euclidean distance
#'
#' @param locs observed location matrix
#' @param locs_pred unobserved location matrix
#' @param refine logical
#'
#' @return list
#'
#' @export
#'
#' @examples
#' 1 + 1
order_maxmin_pred <- function(locs, locs_pred, refine = FALSE){

  #######   obs.pred maxmin ordering   ########

  ord<-1:nrow(locs) #GPvecchia::order_maxmin_exact(locs)
  ord_pred <-GPvecchia::order_maxmin_exact(locs_pred)

  if(refine){

    locs_all = rbind(locs, locs_pred)

    n <- nrow(locs)
    m <- min( round(sqrt(n)), 200 )

    n_pred <- nrow(locs_pred)
    # next is to find 'ord_pred', a maxmin reordering of prediction locations
    NN <- FNN::get.knn( locs_all, k = m )$nn.index
    #NN_pred <- FNN::get.knnx( locs, locs_pred, k = 1 )$nn.dist
    # use ord, then order by NN_pred
    index_in_position <- c( ord, n + ord_pred, rep(NA,n_pred) )
    position_of_index <- order(index_in_position[1:(n+n_pred)])

    # move an index to the end if it is a
    # near neighbor of a previous location
    curlen <- n + n_pred
    nmoved <- 0
    for(j in (n+1):(n+2*n_pred) ){
      # nneigh tells us how many neighbors to look at
      # in order to decide whether the current point
      # has a previously ordered neighbor
      nneigh <- round( min(m,1*(n+n_pred)/(j-nmoved+1)) )
      neighbors <- NN[index_in_position[j],1:nneigh]
      if( min( position_of_index[neighbors], na.rm = TRUE ) < j ){
        nmoved <- nmoved+1
        curlen <- curlen + 1
        position_of_index[ index_in_position[j] ] <- curlen
        index_in_position[curlen] <- index_in_position[j]
        index_in_position[j] <- NA
      }
    }

    ord_pred <- index_in_position[ !is.na( index_in_position ) ][(n+1):(n+n_pred)] - n

  }

  return(list(ord=ord, ord_pred=ord_pred))
}

#' @title Constructing conditioning sets for observed and unobserved locations
#'
#' @param locs location matrix
#' @param m conditioning set size
#' @param fix.first fix.first
#' @param searchmult searchmult
#'
#' @return NNarray
#'
#' @export
#'
#' @examples
#' 1 + 1
find_ordered_nn_pred <- function(locs, m, fix.first = 0, searchmult = 2){

  #######   find NN for prediction locations   ########

  # if locs is a vector, convert to matrix
  if( is.null(ncol(locs)) ){
    locs <- as.matrix(locs)
  }

  # number of locations
  n <- nrow(locs)
  m <- min(m,n-1)
  mult <- 2

  # FNN::get.knnx has strange behavior for exact matches
  # so add a small amount of noise to each location
  ee <- min(apply( locs, 2, stats::sd ))
  locs <- locs + matrix( ee*1e-6*stats::rnorm(n*ncol(locs)), n, ncol(locs) )

  # to store the nearest neighbor indices
  NNarray <- matrix(NA,n,m+1)

  # to the first mult*m+1 by brute force
  maxval <- min( mult*m + 1, n )
  if(fix.first<=maxval){
    NNarray[1:maxval,] <- GpGp::find_ordered_nn_brute(locs[1:maxval,,drop=FALSE],m)
  } else {
    maxval=fix.first
    NNarray[1:(m+1),] <- GpGp::find_ordered_nn_brute(locs[1:(m+1),,drop=FALSE],m)
    NNarray[1:maxval,1]=1:maxval
    NNarray[(m+1):maxval,1+(1:m)]=matrix(rep(1:m,maxval-m),byrow=TRUE,ncol=m)
  }
  query_inds <- min( maxval+1, n):n
  data_inds <- 1:n

  msearch <- m

  while( length(query_inds) > 0 ){
    msearch <- min( max(query_inds), round(searchmult*msearch) )
    data_inds <- 1:min( max(query_inds), n )
    NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
    less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
    sum_less_than_k <- apply(less_than_k,1,sum)
    ind_less_than_k <- which(sum_less_than_k >= m+1)

    NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))

    NNarray[ query_inds[ind_less_than_k], ] <- NN_m

    query_inds <- query_inds[-ind_less_than_k]

  }

  return(NNarray)
}

#' @title Line search for variance correction factor
#'
#' @param fit fit
#' @param m.pred m.pred
#' @param n.test n.test
#' @param scale scale
#' @param scorefun scorefun
#'
#' @return output of optimize()
#'
#' @export
#'
#' @examples
#' 1 + 1
fit_vcf = function(fit, m.pred = 140, n.test = min(1e3, round(nrow(fit$locs)/5)), scale = 'parms', scorefun = lscore)
{
  ##########   line search for variance correction factor   ###

  # remove test data from fit object
  fitsearch=fit
  inds.test=sample(1:nrow(fit$locs),n.test)
  fitsearch$y=fit$y[-inds.test]
  fitsearch$locs=fit$locs[-inds.test,,drop=FALSE]
  fitsearch$X=fit$X[-inds.test,,drop=FALSE]

  # make predictions
  preds=predictions_scaled(fitsearch,locs_pred=fit$locs[inds.test,,drop=FALSE],
                           m=m.pred,joint=FALSE,predvar=TRUE,scale=scale,
                           X_pred=fit$X[inds.test,,drop=FALSE])

  # optimize correction factor
  y.test=fit$y[inds.test]
  objfun=function(vcf) scorefun(y.test,preds$means,preds$vars*vcf)
  vcf=optimize(objfun,c(1e-6,1e6))$minimum

  return(vcf)
}

#' @title log score
#'
#' @param dat dat
#' @param mu mu
#' @param sig2 sig2
#'
#' @return log score
#'
#' @export
#'
#' @examples
#' 1 + 1
lscore = function(dat, mu, sig2)
{
  ### log score

  return( -mean(dnorm(dat,mu,sqrt(sig2),log=TRUE)) )
}

