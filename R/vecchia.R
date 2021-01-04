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



#' @title specify Euclidean/Correlation-based Vecchia approximation with known parameters
#'
#' @param locs A matrix of locations
#' @param m Number of nearby points to condition on
#' @param ordering 'maxmin' or 'coord'
#' @param ordering.method either 'euclidean' or 'correlation'
#' @param coordinate integer or vector of integers in \code{1,...,d}. This argument is used only when ordering is coord
#' @param corr.dist 1-rho, 1-abs(rho), 1-rho^2, sqrt(1-rho), sqrt(1-abs(rho)), or sqrt(1-rho^2)
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
#' ### What you should know ...
#' 
#' locs <- matrix(runif(400), 200, 2) ; m <- 10
#' 
#' rho <- correlationVecchia:::.correlation(locs = locs, covmodel = cov_expo_iso, 
#'                                          covparms = c(1, 0.1), abs.corr = FALSE)
#' 
#' ord <- GPvecchia::order_maxmin_exact(locs)
#' all.equal(ord, 
#'           order_maxmin_correlation_straightforward(locs = locs, 
#'                                                    d = 1 - rho, 
#'                                                    initial.pt = ord[1]))
#' all.equal(ord, 
#'           order_maxmin_correlation_inverseDist(locs = locs, 
#'                                                d = rho, 
#'                                                initial.pt = ord[1]))
#' 
#' locsord <- locs[ord, , drop = FALSE] ; rhoord <- rho[ord, ord]
#' all.equal(GpGp::find_ordered_nn_brute(locs = locsord, m = m), 
#'           conditioning_nn(m = m, d = 1 - rhoord))
#' all.equal(GpGp::find_ordered_nn(locs = locsord, m = m), 
#'           conditioning_nn(m = m, d = 1 - rhoord))
#' 
#' ### Example
#' 
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
#' out01 <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m,
#'                                              ordering = "maxmin", ordering.method = "euclidean",
#'                                              coordinate = NULL, corr.dist = "1-rho", 
#'                                              initial.pt = NULL,
#'                                              conditioning = "NN", 
#'                                              conditioning.method = "euclidean",
#'                                              covmodel = cov_expo_aniso, covparms = covparms)
#' out02 <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m, ordering = "maxmin",
#'                                              ordering.method = "euclidean",
#'                                              coordinate = NULL, corr.dist = "1-rho", 
#'                                              initial.pt = NULL,
#'                                              conditioning = "NN", 
#'                                              conditioning.method = "correlation",
#'                                              covmodel = cov_expo_aniso, covparms = covparms)
#' out03 <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m, ordering = "maxmin",
#'                                              ordering.method = "correlation",
#'                                              coordinate = NULL, corr.dist = "1-rho", 
#'                                              initial.pt = NULL,
#'                                              conditioning = "NN", 
#'                                              conditioning.method = "euclidean",
#'                                              covmodel = cov_expo_aniso, covparms = covparms)
#' out04 <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m,
#'                                              ordering = "maxmin", ordering.method = "correlation",
#'                                              coordinate = NULL, corr.dist = "1-rho", 
#'                                              initial.pt = NULL,
#'                                              conditioning = "NN", 
#'                                              conditioning.method = "correlation",
#'                                              covmodel = cov_expo_aniso, covparms = covparms)
#' out05 <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m,
#'                                              ordering = "coord", ordering.method = "euclidean",
#'                                              coordinate = NULL, corr.dist = "1-rho", 
#'                                              initial.pt = NULL,
#'                                              conditioning = "NN", 
#'                                              conditioning.method = "euclidean",
#'                                              covmodel = cov_expo_aniso, covparms = covparms)
#' out06 <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m,
#'                                              ordering = "coord", ordering.method = "euclidean",
#'                                              coordinate = c(1), corr.dist = "1-rho", 
#'                                              initial.pt = NULL,
#'                                              conditioning = "NN", 
#'                                              conditioning.method = "euclidean",
#'                                              covmodel = cov_expo_aniso, covparms = covparms)
#' out07 <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m,
#'                                              ordering = "coord", ordering.method = "euclidean",
#'                                              coordinate = c(2), corr.dist = "1-rho", 
#'                                              initial.pt = NULL,
#'                                              conditioning = "NN", 
#'                                              conditioning.method = "euclidean",
#'                                              covmodel = cov_expo_aniso, covparms = covparms)
#' 
#' kls.euc.euc    <- performance(out = out01, locs = locs,
#'                                covmodel = cov_expo_aniso, covparms = covparms)
#' kls.cor.euc    <- performance(out = out02, locs = locs,
#'                                covmodel = cov_expo_aniso, covparms = covparms)
#' kls.euc.cor    <- performance(out = out03, locs = locs,
#'                                covmodel = cov_expo_aniso, covparms = covparms)
#' kls.cor.cor    <- performance(out = out04, locs = locs,
#'                                covmodel = cov_expo_aniso, covparms = covparms)
#' kls.coord.xy   <- performance(out = out05, locs = locs,
#'                                covmodel = cov_expo_aniso, covparms = covparms)
#' kls.coord.x    <- performance(out = out06, locs = locs,
#'                                covmodel = cov_expo_aniso, covparms = covparms)
#' kls.coord.y    <- performance(out = out07, locs = locs,
#'                                covmodel = cov_expo_aniso, covparms = covparms)
#' 
#' c(kls.euc.euc, kls.euc.cor, kls.cor.euc, kls.cor.cor, kls.coord.xy, kls.coord.x, kls.coord.y)
corrvecchia_specify_knownCovparms_2 <- function(locs, m, ordering = "maxmin", ordering.method = "correlation", coordinate = NULL, corr.dist = "1-rho", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel, covparms) 
{
  locs  <- as.matrix(locs)
  p     <- ncol(locs)
  n     <- nrow(locs)
  
  if(ordering.method == "correlation" | conditioning.method == "correlation") {
    
    rho <- .correlation(locs = locs, covmodel = covmodel, covparms = covparms, abs.corr = FALSE)
    
  }
  
  if(ordering == "coord") {
    
    # ordering
    if(is.null(coordinate)) coordinate <- seq(p)
    ord       <- order_coordinate(locs = locs, coordinate = coordinate)
    
    # to order locations
    locsord   <- locs[ord, , drop = FALSE]
    
    if(is.matrix(covmodel)) covmodel <- covmodel[ord, ord]
    
    # conditioning
    if(conditioning.method == "euclidean") {                                                    # ?-Coord + E-NN
      
      cond.sets   <- GpGp::find_ordered_nn(locs = locsord, m = m)
      
    } else if(conditioning.method == "correlation") {                                           # ?-Coord + C-NN
      
      rho         <- rho[ord, ord]
      cond.sets   <- .find_ordered_cnn(m = m, rho = rho, corr.dist = corr.dist)
      
    } else {
      
      stop("The argument conditioning.method must be euclidean or correlation. Please try again!")
      
    }
    
  } else if(ordering == "maxmin") {
    
    if(ordering.method == "euclidean" & conditioning.method == "euclidean") {                   # E-maxmin + E-NN
      
      # ordering
      ord       <- order_maxmin_euclidean(locs = locs)
      
      # to order locations
      locsord   <- locs[ord, , drop = FALSE]
      
      if(is.matrix(covmodel)) covmodel <- covmodel[ord, ord]
      
      # conditioning
      cond.sets <- GpGp::find_ordered_nn_brute(locs = locsord, m = m)
      
    } else if(ordering.method == "correlation" & conditioning.method == "euclidean") {          # C-maxmin + E-NN
      
      # ordering
      ord       <- .order_cmaxmin(locs = locs, rho = rho, initial.pt = initial.pt, corr.dist = corr.dist)
      
      # to order locations
      locsord   <- locs[ord, , drop = FALSE]
      
      if(is.matrix(covmodel)) covmodel <- covmodel[ord, ord]
      
      # conditioning
      cond.sets <- GpGp::find_ordered_nn_brute(locs = locsord, m = m)
      
    } else if(ordering.method == "euclidean" & conditioning.method == "correlation") {          # E-maxmin + C-NN
      
      # ordering
      ord       <- order_maxmin_euclidean(locs = locs)
      
      # to order locations
      locsord   <- locs[ord, , drop = FALSE]
      
      if(is.matrix(covmodel)) covmodel <- covmodel[ord, ord]
      
      # conditioning
      rho         <- rho[ord, ord]
      cond.sets   <- .find_ordered_cnn(m = m, rho = rho, corr.dist = corr.dist)
      
    } else if(ordering.method == "correlation" & conditioning.method == "correlation") {        # C-maxmin + C-NN
      
      # ordering
      ord       <- .order_cmaxmin(locs = locs, rho = rho, initial.pt = initial.pt, corr.dist = corr.dist)
      
      # to order locations
      locsord   <- locs[ord, , drop = FALSE]
      
      if(is.matrix(covmodel)) covmodel <- covmodel[ord, ord]
      
      # conditioning
      rho         <- rho[ord, ord]
      cond.sets   <- .find_ordered_cnn(m = m, rho = rho, corr.dist = corr.dist)
      
    } else {
      
      stop("Both the arguments ordering.method and conditioning.method must be euclidean or correlation. Please try again!")
      
    }
    
  } else {
    
    stop("The argument ordering must be coord or maxmin. Please try again!")
    
  }
  
  ### return
  Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)
  
  vecchia.approx <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = ordering, ordering.method = ordering.method, conditioning = 'NN', conditioning.method = conditioning.method)
  return(vecchia.approx)
}

.order_cmaxmin <- function(locs, rho, initial.pt, corr.dist)
{
  if(corr.dist == "1-rho") {
    # ord       <- order_maxmin_correlation_straightforward(locs = locs, d = 1 - rho, initial.pt = initial.pt)
    ord       <- order_maxmin_correlation_inverseDist(locs = locs, d.inv = rho, initial.pt = initial.pt)
  } else if(corr.dist == "1-abs(rho)") {
    # ord       <- order_maxmin_correlation_straightforward(locs = locs, d = 1 - abs(rho), initial.pt = initial.pt)
    ord       <- order_maxmin_correlation_inverseDist(locs = locs, d.inv = abs(rho), initial.pt = initial.pt)
  } else if(corr.dist == "1-rho^2") {
    # ord       <- order_maxmin_correlation_straightforward(locs = locs, d = 1 - rho^2, initial.pt = initial.pt)
    ord       <- order_maxmin_correlation_inverseDist(locs = locs, d.inv = abs(rho), initial.pt = initial.pt)
  } else if(corr.dist == "sqrt(1-rho)") {
    # ord       <- order_maxmin_correlation_straightforward(locs = locs, d = sqrt(1 - rho), initial.pt = initial.pt)
    ord       <- order_maxmin_correlation_inverseDist(locs = locs, d.inv = rho, initial.pt = initial.pt)
  } else if(corr.dist == "sqrt(1-abs(rho))") {
    # ord       <- order_maxmin_correlation_straightforward(locs = locs, d = sqrt(1 - abs(rho)), initial.pt = initial.pt)
    ord       <- order_maxmin_correlation_inverseDist(locs = locs, d.inv = abs(rho), initial.pt = initial.pt)
  } else if(corr.dist == "sqrt(1-rho^2)") {
    # ord       <- order_maxmin_correlation_straightforward(locs = locs, d = sqrt(1 - rho^2), initial.pt = initial.pt)
    ord       <- order_maxmin_correlation_inverseDist(locs = locs, d.inv = abs(rho), initial.pt = initial.pt)
  } else {
    stop("The argument corr.dist must be one of the followings: 1-rho, 1-abs(rho), 1-rho^2, sqrt(1-rho), sqrt(1-abs(rho)), and sqrt(1-rho^2).")
  }
  
  return(ord)
}

.find_ordered_cnn <- function(m, rho, corr.dist) 
{
  if(corr.dist == "1-rho") {
    # cond.sets <- conditioning_nn(m = m, d = 1 - rho)
    cond.sets <- conditioning_nn_Rcpp(m = m, d = 1- rho) + 1
  } else if(corr.dist == "1-abs(rho)") {
    # cond.sets <- conditioning_nn(m = m, d = 1 - abs(rho))
    cond.sets <- conditioning_nn_Rcpp(m = m, d = 1- abs(rho)) + 1
  } else if(corr.dist == "1-rho^2") {
    # cond.sets <- conditioning_nn(m = m, d = 1 - rho^2)
    # cond.sets <- conditioning_nn(m = m, d = 1 - abs(rho))
    cond.sets <- conditioning_nn_Rcpp(m = m, d = 1- abs(rho)) + 1
  } else if(corr.dist == "sqrt(1-rho)") {
    # cond.sets <- conditioning_nn(m = m, d = sqrt(1 - rho))
    # cond.sets <- conditioning_nn(m = m, d = 1 - rho)
    cond.sets <- conditioning_nn_Rcpp(m = m, d = 1- rho) + 1
  } else if(corr.dist == "sqrt(1-abs(rho))") {
    # cond.sets <- conditioning_nn(m = m, d = sqrt(1 - abs(rho)))
    # cond.sets <- conditioning_nn(m = m, d = 1 - abs(rho))
    cond.sets <- conditioning_nn_Rcpp(m = m, d = 1- abs(rho)) + 1
  } else if(corr.dist == "sqrt(1-rho^2)") {
    # cond.sets <- conditioning_nn(m = m, d = sqrt(1 - rho^2))
    # cond.sets <- conditioning_nn(m = m, d = 1 - abs(rho))
    cond.sets <- conditioning_nn_Rcpp(m = m, d = 1- abs(rho)) + 1
  } else {
    stop("The argument corr.dist must be one of the followings: 1-rho, 1-abs(rho), 1-rho^2, sqrt(1-rho), sqrt(1-abs(rho)), and sqrt(1-rho^2).")
  }
  
  return(cond.sets)
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
#' @param initial.pt If \code{NULL} then the most centered location is selected as the first location. If an integer then it is used as the index of the first location
#'
#' @return A vector of indices giving the euclidean-based maxmin ordering 
#' @export
#'
#' @examples
#' locs <- matrix(runif(100 * 2), 100, 2)
#' 
#' identical(order_maxmin_euclidean(locs), GPvecchia::order_maxmin_exact(locs))
order_maxmin_euclidean <- function(locs, initial.pt = NULL)
{
  n             <- nrow(locs)
  p             <- ncol(locs)
  ord           <- rep(NA, n)
  cen           <- t(as.matrix(colMeans(locs)))
  
  if(is.null(initial.pt)) {
    ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
    cand.argmax   <- seq(n)[seq(n) != ord[1]]
  } else if(is.numeric(initial.pt) & length(initial.pt) == 1) {
    ord[1]        <- initial.pt
    cand.argmax   <- seq(n)[seq(n) != ord[1]]
  } else {
    stop("The argument initial.pt must be an integer.")
  }
  
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
#' rho <- correlationVecchia:::.correlation(locs, cov_expo_iso, covparms, FALSE)
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
  
  # second step
  ind           <- as.matrix(expand.grid(cand.argmax, ord[1]))
  cdist         <- d[ind]
  ord[2]        <- cand.argmax[which.max(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]
  
  # third and so on ...
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
#' rho <- correlationVecchia:::.correlation(locs, cov_expo_iso, covparms, FALSE)
#' identical(order_maxmin_correlation_inverseDist(locs, rho, 
#'                        initial.pt = GPvecchia::order_maxmin_exact(locs)[1]), 
#'           order_maxmin_correlation_straightforward(locs, 1 - rho, 
#'                        initial.pt = GPvecchia::order_maxmin_exact(locs)[1]))
#' 
#' covparms <- c(1, 0.1, 10)
#' locs.trans <- cbind(locs[ ,1] * covparms[3], locs[,2]) / covparms[2]
#' 
#' rho.new <- correlationVecchia:::.correlation(locs, cov_expo_aniso, covparms, FALSE)
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
  
  # second step
  ind           <- as.matrix(expand.grid(cand.argmax, ord[1]))
  cdist         <- d.inv[ind]
  ord[2]        <- cand.argmax[which.min(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]
  
  # third and so on ...
  for(j in 3:(n-1)){
    ind               <- as.matrix(expand.grid(cand.argmax, ord[seq(j-1)]))
    cdist             <- matrix(d.inv[ind], nrow = length(cand.argmax), ncol = j-1, byrow = F)
    cdist             <- Rfast::rowMaxs(cdist, value = TRUE)
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
#' 
#' n <- 10 ; m <- 3
#' 
#' locs <- matrix(runif(n * 2), n, 2)
#' d <- as.matrix(dist(locs))
#' 
#' conditioning_nn(m, d)
#' conditioning_nn_Rcpp(m, d) + 1
#' 
#' require(microbenchmark)
#' microbenchmark(conditioning_nn(m, d), conditioning_nn_Rcpp(m, d), times = 2)
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
