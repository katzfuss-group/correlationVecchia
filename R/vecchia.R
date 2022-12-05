####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes R functions to implement vecchia approximations.
###
###   Contents:
###     cvecchia_specify / cvecchia_m_specify
###     order_coordinate / order_maxmin_euclidean / order_maxmin_correlation
###     cvecchia_rho_specify / rhoBased_sortSparse / rhoBased_sortSparse_rev
###
###   Note: These functions are heavily based on Rcpp functions in vecchia_cpp.cpp.
###
####################################################################################

#' @title specify Euclidean/Correlation-based Vecchia approximation with known parameters
#'
#' @param locs A matrix of locations
#' @param m Number of nearby points to condition on
#' @param rho Hyperparameter rho (instead of m)
#' @param initial.pt Check vecchia_m_specify() and vecchia_rho_specify() functions' description
#' @param covmodel Covariance model (or matrix)
#' @param covparms Covariance parameters as a vector. The first entry must be its overall variance
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|rho|. If \code{FALSE} then distane = 1-rho
#' @param coordinate integer or vector of integers in \code{1,...,d}. This argument is used only when ordering is coord
#' @param ordering 'MM' or 'coord'
#' @param ordering.method 'euc' or 'cor'
#' @param conditioning 'NN' (nearest neighbor)
#' @param conditioning.method 'euc' or 'cor'
#'
#' @return An object that specifies the vecchia approximation for later use in likelihood evaluation or prediction.
#'
#' @export
#'
#' @examples
#' n             <- 20^2
#' m             <- 20
#' locs          <- matrix(runif(n * 2, 0, 1), n, 2)
#' covparms      <- c(1, 0.1, 10)
#'
#' # true cov matrix
#' covmat <- cov_expo_aniso(locs, covparms)
#'
#' # Visualize the process
#' y <- as.numeric(t(chol(covmat)) %*% rnorm(n))
#' fields::quilt.plot(locs[,1], locs[,2], y)
#'
#' out01 <- cvecchia_specify(locs = locs, m = 18, initial.pt = NULL,
#'                           covmodel = cov_expo_aniso, covparms = covparms,
#'                           abs.corr = FALSE, coordinate = NULL,
#'                           ordering = "MM", ordering.method = "euc",
#'                           conditioning = "NN", conditioning.method = "euc")
#'
#' out02 <- cvecchia_specify(locs = locs, rho = 2.1, initial.pt = NULL,
#'                           covmodel = cov_expo_aniso, covparms = covparms,
#'                           abs.corr = FALSE, coordinate = NULL,
#'                           ordering = "MM", ordering.method = "cor",
#'                           conditioning = "NN", conditioning.method = "cor")
#'
#' kls.mbased    <- performance(vecchia.approx = out01, locs = locs,
#'                              covmodel = cov_expo_aniso, covparms = covparms)
#' kls.rhobased  <- performance(vecchia.approx = out02, locs = locs,
#'                              covmodel = cov_expo_aniso, covparms = covparms)
#'
#' barplot(log10(c(kls.mbased, kls.rhobased)),
#'         names.arg = c("m-based CVecchia", "rho-based CVecchia"),
#'         main = "Vecchia Approximations", ylab = "log10-scale KL divergence")
#'
#' sum(!is.na(out01$U.prep$revNNarray)) / n^2
#' sum(!is.na(out02$U.prep$revNNarray)) / n^2
cvecchia_specify <- function(locs, m = NULL, rho = NULL, initial.pt = NULL, covmodel, covparms, abs.corr = FALSE, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor")
{
  if(!is.null(m) & is.null(rho)) {

    vecchia.approx  <- cvecchia_m_specify(locs = locs, m = m, initial.pt = initial.pt, covmodel = covmodel, covparms = covparms, abs.corr = abs.corr, coordinate = coordinate, ordering = ordering, ordering.method = ordering.method, conditioning = conditioning, conditioning.method = conditioning.method)

  } else if(is.null(m) & !is.null(rho)) {

    if(is.matrix(covmodel)) stop("If you want to use rho-based Vecchia, you must specify your covariance model for now!")

    if(is.function(covmodel)) covmodel <- deparse(substitute(covmodel))

    # fname <- function(covmodel)
    # {
    #   covmodel <- deparse(substitute(covmodel))
    #   return(covmodel)
    # }
    #
    # fname(cov_expo_aniso)

    vecchia.approx  <- cvecchia_rho_specify(locs = locs, rho = rho, initial.pt = initial.pt, covmodel = covmodel, covparms = covparms, abs.corr = abs.corr, coordinate = coordinate, ordering = ordering, ordering.method = ordering.method, conditioning = conditioning, conditioning.method = conditioning.method)

  } else {

    stop("You must choose m or rho!")
  }

  return(vecchia.approx)
}

####################################################################################

#' @title specify Euclidean/Correlation-based Vecchia approximation with known parameters
#'
#' @param locs A matrix of locations
#' @param m Number of nearby points to condition on
#' @param initial.pt NULL = which.min(rowMeans(d)), center = euclidean-based center, integer = specify the first obs, 'random' = at random, and else = which.min(rowMeans(d))
#' @param covmodel Covariance model (or matrix)
#' @param covparms Covariance parameters as a vector. The first entry must be its overall variance
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|rho|. If \code{FALSE} then distane = 1-rho
#' @param coordinate integer or vector of integers in \code{1,...,d}. This argument is used only when ordering is coord
#' @param ordering 'MM' or 'coord'
#' @param ordering.method 'euc' or 'cor'
#' @param conditioning 'NN' (nearest neighbor)
#' @param conditioning.method 'euc' or 'cor'
#'
#' @return An object that specifies the vecchia approximation for later use in likelihood evaluation or prediction.
#'
#' @export
#'
#' @examples
#' locs    <- matrix(runif(400), 200, 2)
#' m       <- 10
#' cormat  <- cov_expo_iso(locs = locs, covparms = c(1, 0.1))
#'
#' ord     <- GPvecchia::order_maxmin_exact(locs)
#' all.equal(ord, order_maxmin_euclidean(locs = locs, initial.pt = ord[1]))
#' all.equal(ord, order_maxmin_correlation(locs = locs, d.inv = cormat, initial.pt = ord[1]))
#'
#' locsord <- locs[ord, , drop = FALSE]
#' corord  <- cormat[ord, ord]
#'
#' all.equal(GpGp::find_ordered_nn_brute(locs = locsord, m = m),
#'           conditioning_m_Rcpp(m = m, d = 1 - corord) + 1)
#' all.equal(GpGp::find_ordered_nn(locs = locsord, m = m),
#'           conditioning_m_Rcpp(m = m, d = 1 - corord) + 1)
#'
#' ### Example
#'
#' n             <- 20^2
#' m             <- 20
#' locs          <- matrix(runif(n * 2, 0, 1), n, 2)
#' covparms      <- c(1, 0.1, 10)
#'
#' # true cov matrix
#' covmat <- cov_expo_aniso(locs, covparms)
#'
#' # Visualize the process
#' y <- as.numeric(t(chol(covmat)) %*% rnorm(n))
#' fields::quilt.plot(locs[,1], locs[,2], y)
#'
#' out01 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
#'                             covmodel = covmat, covparms = covparms,
#'                             abs.corr = FALSE, coordinate = c(1),
#'                             ordering = "coord", ordering.method = "euc",
#'                             conditioning = "NN", conditioning.method = "euc")
#'
#' out02 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
#'                             covmodel = covmat, covparms = covparms,
#'                             abs.corr = FALSE, coordinate = c(2),
#'                             ordering = "coord", ordering.method = "euc",
#'                             conditioning = "NN", conditioning.method = "euc")
#'
#' out03 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
#'                             covmodel = covmat, covparms = covparms,
#'                             abs.corr = FALSE, coordinate = NULL,
#'                             ordering = "MM", ordering.method = "euc",
#'                             conditioning = "NN", conditioning.method = "euc")
#'
#' out04 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
#'                             covmodel = covmat, covparms = covparms,
#'                             abs.corr = FALSE, coordinate = NULL,
#'                             ordering = "MM", ordering.method = "euc",
#'                             conditioning = "NN", conditioning.method = "cor")
#'
#' out05 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
#'                             covmodel = covmat, covparms = covparms,
#'                             abs.corr = FALSE, coordinate = NULL,
#'                             ordering = "MM", ordering.method = "cor",
#'                             conditioning = "NN", conditioning.method = "euc")
#'
#' out06 <- cvecchia_m_specify(locs = locs, m = m, initial.pt = NULL,
#'                             covmodel = covmat, covparms = covparms,
#'                             abs.corr = FALSE, coordinate = NULL,
#'                             ordering = "MM", ordering.method = "cor",
#'                             conditioning = "NN", conditioning.method = "cor")
#'
#' kls.coord.x    <- performance(vecchia.approx = out01, locs = locs,
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.coord.y    <- performance(vecchia.approx = out02, locs = locs,
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.euc.euc    <- performance(vecchia.approx = out03, locs = locs,
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.euc.cor    <- performance(vecchia.approx = out04, locs = locs,
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.cor.euc    <- performance(vecchia.approx = out05, locs = locs,
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.cor.cor    <- performance(vecchia.approx = out06, locs = locs,
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#'
#' barplot(log10(c(kls.coord.x, kls.coord.y, kls.euc.euc,
#'                 kls.euc.cor, kls.cor.euc, kls.cor.cor)),
#'         names.arg = c("X-coord + E-NN", "Y-coord + E-NN",
#'                       "E-MM + E-NN", "E-MM + C-NN",
#'                       "C-MM + E-NN", "C-MM + C-NN"),
#'         main = "Vecchia Approximations", ylab = "log10-scale KL divergence")
cvecchia_m_specify <- function(locs, m, initial.pt = NULL, covmodel, covparms, abs.corr = FALSE, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor")
{
  ### checkargs
  if(conditioning != "NN") stop("You must use NN conditioning for now!")

  ### initialization
  locs  <- as.matrix(locs)
  p     <- ncol(locs)
  n     <- nrow(locs)

  if(ordering.method %in% c("cor", "correlation") | conditioning.method %in% c("cor", "correlation")) {

    covmodel <- .correlation(locs = locs, covmodel = covmodel, covparms = covparms, abs.corr = abs.corr)
  }

  ### body
  if( ordering %in% c("coord", "coordinate") ) {

    # ordering
    if(is.null(coordinate)) coordinate <- seq(p)
    ord       <- order_coordinate(locs = locs, coordinate = coordinate)

    # to order locations
    locsord   <- locs[ord, , drop = FALSE]

    # conditioning
    if( conditioning.method %in% c("euclidean", "euc") ) {                                                                # ?-Coord + E-NN

      cond.sets   <- .find_ordered_NN(m = m, locsord = locsord, corord = NULL, method = "euc", corr.dist = NULL)

    } else if( conditioning.method %in% c("correlation", "cor") ) {                                                       # ?-Coord + C-NN

      covmodel    <- covmodel[ord, ord]
      cond.sets   <- .find_ordered_NN(m = m, locsord = NULL, corord = covmodel, method = "cor", corr.dist = "sqrt(1-cor)")

    } else {

      stop("The argument conditioning.method must be euc (euclidean) or cor (correlation). Please try again!")
    }

  } else if( ordering %in% c("MM", "maxmin") ){

    if( ordering.method %in% c("euc", "euclidean") & conditioning.method %in% c("euc", "euclidean") ) {                   # E-MM + E-NN

      # ordering
      ord       <- .find_MM(locs = locs, initial.pt = initial.pt, method = "euc", cormat = NULL, corr.dist = NULL)

      # to order locations
      locsord   <- locs[ord, , drop = FALSE]

      # conditioning
      cond.sets <- .find_ordered_NN(m = m, locsord = locsord, corord = NULL, method = "euc", corr.dist = NULL)

    } else if( ordering.method %in% c("cor", "correlation") & conditioning.method %in% c("euc", "euclidean") ) {          # C-MM + E-NN

      # ordering
      ord       <- .find_MM(locs = locs, initial.pt = initial.pt, method = "cor", cormat = covmodel, corr.dist = "sqrt(1-cor)")

      # to order locations
      locsord   <- locs[ord, , drop = FALSE]

      # conditioning
      cond.sets <- .find_ordered_NN(m = m, locsord = locsord, corord = NULL, method = "euc", corr.dist = NULL)

    } else if( ordering.method %in% c("euc", "euclidean") & conditioning.method %in% c("cor", "correlation") ) {          # E-MM + C-NN

      # ordering
      ord       <- .find_MM(locs = locs, initial.pt = initial.pt, method = "euc", cormat = NULL, corr.dist = NULL)

      # to order locations
      locsord   <- locs[ord, , drop = FALSE]
      covmodel  <- covmodel[ord, ord]

      # conditioning
      cond.sets <- .find_ordered_NN(m = m, locsord = NULL, corord = covmodel, method = "cor", corr.dist = "sqrt(1-cor)")

    } else if( ordering.method %in% c("cor", "correlation") & conditioning.method %in% c("cor", "correlation") ) {        # C-MM + C-NN

      # ordering
      ord       <- .find_MM(locs = locs, initial.pt = initial.pt, method = "cor", cormat = covmodel, corr.dist = "sqrt(1-cor)")

      # to order locations
      locsord   <- locs[ord, , drop = FALSE]
      covmodel  <- covmodel[ord, ord]

      # conditioning
      cond.sets <- .find_ordered_NN(m = m, locsord = NULL, corord = covmodel, method = "cor", corr.dist = "sqrt(1-cor)")

    } else {

      stop("Both the arguments ordering.method and conditioning.method must be euc (euclidean) or cor (correlation). Please try again!")
    }

  } else {

    stop("The argument ordering must be coord (coordinate) or MM (maxmin). Please try again!")
  }

  ### return
  Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)

  vecchia.approx <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = ordering, ordering.method = ordering.method, conditioning = 'NN', conditioning.method = conditioning.method)
  return(vecchia.approx)
}

.find_MM <- function(locs, initial.pt = NULL, method = "cor", cormat = NULL, corr.dist = NULL)
{
  # # example
  #
  # locs        <- matrix(runif(30), 15, 2)
  # covparms    <- c(1, 0.1, 1)
  # cormat      <- cov_expo_aniso(locs, covparms)
  # initial.pt  <- 1
  #
  # ord1        <- correlationVecchia:::.find_MM(locs, initial.pt, "euc", cormat, NULL) ; ord1
  # ord2        <- correlationVecchia:::.find_MM(locs, initial.pt, "cor", cormat, "1-cor") ; ord2
  #
  # locsord     <- locs[ord2, ]
  # corord      <- cormat[ord2, ord2]
  #
  # cond1       <- correlationVecchia:::.find_ordered_NN(5, locsord, corord, "euc", NULL) ; tail(cond1)
  # cond2       <- correlationVecchia:::.find_ordered_NN(5, locsord, corord, "cor", "1-cor") ; tail(cond2)
  #
  # cond3       <- GpGp::find_ordered_nn_brute(locsord, 5)
  # cond4       <- GpGp::find_ordered_nn(locsord, 5)
  #
  # cond1 == cond2
  # cond1 == cond3
  # cond1 == cond4

  if(method %in% c("euclidean", "euc", "euclidean-based")) {

    ord       <- order_maxmin_euclidean(locs = locs, initial.pt = initial.pt)

  } else if(method %in% c("correlation", "cor", "correlation-based")) {

    if(is.null(corr.dist)) {

      stop("If method is cor, then the argument corr.dist must be provided!")

    } else if( corr.dist %in% c("1-cor", "sqrt(1-cor)") ) {

      ord       <- order_maxmin_correlation(locs = locs, d.inv = cormat, initial.pt = initial.pt, method = "inverseDist")

    } else if( corr.dist %in% c("1-abs(cor)", "1-cor^2", "sqrt(1-abs(cor))", "sqrt(1-cor^2)") ) {

      ord       <- order_maxmin_correlation(locs = locs, d.inv = abs(cormat), initial.pt = initial.pt, method = "inverseDist")

    } else {

      stop("The argument corr.dist must be one of the followings: 1-cor, 1-abs(cor), 1-cor^2, sqrt(1-cor), sqrt(1-abs(cor)), and sqrt(1-cor^2).")
    }

  } else {

    stop("Check the argument method! It must be euc (euclidean) or cor (correlation).")
  }

  return(ord)
}

.find_ordered_NN <- function(m, locsord = NULL, corord = NULL, method = "cor", corr.dist = NULL)
{
  # # example
  #
  # locs        <- matrix(runif(30), 15, 2)
  # covparms    <- c(1, 0.1, 1)
  # cormat      <- cov_expo_aniso(locs, covparms)
  # initial.pt  <- 1
  #
  # ord1        <- correlationVecchia:::.find_MM(locs, initial.pt, "euc", cormat, NULL) ; ord1
  # ord2        <- correlationVecchia:::.find_MM(locs, initial.pt, "cor", cormat, "1-cor") ; ord2
  #
  # locsord     <- locs[ord2, ]
  # corord      <- cormat[ord2, ord2]
  #
  # cond1       <- correlationVecchia:::.find_ordered_NN(5, locsord, corord, "euc", NULL) ; tail(cond1)
  # cond2       <- correlationVecchia:::.find_ordered_NN(5, locsord, corord, "cor", "1-cor") ; tail(cond2)
  #
  # cond3       <- GpGp::find_ordered_nn_brute(locsord, 5)
  # cond4       <- GpGp::find_ordered_nn(locsord, 5)
  #
  # cond1 == cond2
  # cond1 == cond3
  # cond1 == cond4

  if(method %in% c("euclidean", "euc", "euclidean-based")) {

    if( is.null(locsord) ) stop("If method = euc, then locsord must be provided!")

    cond.sets <- GpGp::find_ordered_nn_brute(locs = locsord, m = m)

  } else if(method %in% c("correlation", "cor", "correlation-based")) {

    if( is.null(corord) ) stop("If emthod = cor, then corord must be provided!")

    if(is.null(corr.dist)) {

      stop("If method is cor, then the argument corr.dist must be provided!")

    } else if( corr.dist %in% c("1-cor", "sqrt(1-cor)") ) {

      cond.sets <- conditioning_m_Rcpp(m = m, d = 1 - corord) + 1

    } else if( corr.dist %in% c("1-abs(cor)", "1-cor^2", "sqrt(1-abs(cor))", "sqrt(1-cor^2)") ) {

      cond.sets <- conditioning_m_Rcpp(m = m, d = 1 - abs(corord)) + 1

    } else {

      stop("The argument corr.dist must be one of the followings: 1-cor, 1-abs(cor), 1-cor^2, sqrt(1-cor), sqrt(1-abs(cor)), and sqrt(1-cor^2).")
    }

  } else {

    stop("Check the argument method! It must be euc (euclidean) or cor (correlation).")
  }

  return(cond.sets)
}

.cov2cor <- function(covmodel, locs = NULL, covparms = NULL)
{
  if(is.function(covmodel)) {

    if(is.null(locs) | is.null(covparms)) stop("If covmodel is a function, then both locs and covparms must be provided!")

    covmodel  <- covmodel(locs = locs, covparms = covparms)

  } else if(is.character(covmodel)) {

    covmodel  <- get(covmodel)

    if(is.null(locs) | is.null(covparms)) stop("If covmodel is a character, then both locs and covparms must be provided!")

    covmodel  <- covmodel(locs = locs, covparms = covparms)

  } else if(is.matrix(covmodel)) {

    # desirable!!!

  } else {

    stop("Check the argument covmodel!")
  }

  n         <- nrow(covmodel)
  diagvec   <- diag(covmodel)

  covmodel  <- covmodel / sqrt( matrix(diagvec, n, n, byrow = TRUE) * matrix(diagvec, n, n, byrow = FALSE) )

  return( covmodel )
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

#' @title Maximum-minimum-distance (maxmin) ordering of locations with respect to a user-defined distance matrix
#'
#' @param locs A matrix of locations
#' @param d A matrix of distances between locations (for instance, 1-cor)
#' @param d.inv A matrix of inverse of distances between locations (for instance, cor)
#' @param initial.pt NULL = which.min(rowMeans(d)), center = euclidean-based center, integer = specify the first obs, 'random' = at random, and else = which.min(rowMeans(d))
#' @param method "straightforward" or "inverseDist"
#'
#' @return  A vector of indices giving the maxmin ordering with respect to the user-defined distance matrix
#'
#' @export
#'
#' @examples
#' locs        <- matrix(runif(15 * 2), 15, 2)
#' covparms    <- c(1, 0.1)
#' cormat      <- cov_expo_iso(locs, covparms)
#'
#' identical(order_maxmin_euclidean(locs),
#'           order_maxmin_correlation(locs = locs, d = 1 - cormat,
#'                                    initial.pt = order_maxmin_euclidean(locs)[1],
#'                                    method = "straightforward"),
#'           order_maxmin_correlation(locs = locs, d.inv = cormat,
#'                                    initial.pt = order_maxmin_euclidean(locs)[1],
#'                                    method = "inverseDist"))
#'
#' locs        <- matrix(runif(15 * 2), 15, 2)
#' covparms    <- c(1, 0.1, 10)
#' locs.trans  <- cbind(locs[ ,1] * covparms[3], locs[,2]) / covparms[2]
#' cormat      <- cov_expo_aniso(locs = locs, covparms = covparms)
#'
#' identical(order_maxmin_euclidean(locs.trans),
#'           order_maxmin_correlation(locs = locs, d.inv = cormat,
#'                                    initial.pt = order_maxmin_euclidean(locs.trans)[1],
#'                                    method = "inverseDist"))
#'
#' # CAUTION: This function can cause numerical issue
#' identical(order_maxmin_euclidean(locs.trans),
#'           order_maxmin_correlation(locs = locs, d = 1 - cormat,
#'                                    initial.pt = order_maxmin_euclidean(locs.trans)[1],
#'                                    method = "straightforward"))
#' # Compare to sortSparse_Rcpp
#' locs         <- matrix(runif(30), 15, 2)
#' cormat       <- cov_expo_aniso(locs, c(1, 0.1, 1.0))
#' initial.pt   <- 1
#'
#' order_maxmin_correlation(locs = locs, d.inv = cormat,
#'                          initial.pt = initial.pt, method = "inverseDist")
#' sortSparse_Rcpp(locs, 2.01, initial.pt-1, "correlation",
#'                 "cov_expo_aniso", c(1, 0.1, 1.0))$P + 1
order_maxmin_correlation <- function(locs, d = NULL, d.inv = NULL, initial.pt = NULL, method = "inverseDist")
{
  if(method %in% c("inverseDist", "d.inv", "distance", "dist", "inv", "inverse")) {

    if(is.null(d.inv)) stop("If method = inverseDist, then the argument d must be provided!")

    ord <- .order_maxmin_correlation_inverseDist(locs = locs, d.inv = d.inv, initial.pt = initial.pt)

  } else if(method %in% c("straightforward", "d", "covariance", "cov", "correlation", "cor")) {

    if(is.null(d)) stop("If emthod = straightforward, then the argument d must be provided!")

    ord <- .order_maxmin_correlation_straightforward(locs = locs, d = d, initial.pt = initial.pt)

  } else {

    stop("Check the argument method!")
  }

  return( ord )
}

.order_maxmin_correlation_straightforward <- function(locs, d, initial.pt)
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

.order_maxmin_correlation_inverseDist <- function(locs, d.inv, initial.pt)
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

.distance_correlation <- function(locs, covmodel, covparms = NULL, abs.corr)
{
  ### CAUTION: This function can cause numerical issue. Please use the '.correlation()' function, instead.

  if(is.function(covmodel)) {

    if(is.null(covparms)) stop("If covmodel is a function, then covparms must be provided!")

    if(abs.corr == FALSE) {

      dist.matrix <- 1 - .cov2cor( covmodel(locs, covparms) ) # 1 - cor

    } else if(abs.corr == TRUE) {

      dist.matrix <- 1 - abs( .cov2cor( covmodel(locs, covparms) ) ) # 1-|cor|

    } else {

      stop("abs.corr must be logical.")
    }

  } else if(is.matrix(covmodel)) {

    if(abs.corr == FALSE) {

      dist.matrix <- 1 - .cov2cor( covmodel ) # 1 - cor

    } else if(abs.corr == TRUE) {

      dist.matrix <- 1 - abs( .cov2cor(covmodel) ) # 1-|cor|

    } else {

      stop("abs.corr must be logical.")
    }

  } else {

    stop("covmodel must be either function or matrix.")

  }

  return(dist.matrix)
}

.correlation <- function(locs, covmodel, covparms = NULL, abs.corr)
{
  if(is.function(covmodel)) {

    if(is.null(covparms)) stop("If covmodel is a function, then covparms must be provided!")

    if(abs.corr == FALSE) {

      corr.matrix <- .cov2cor( covmodel(locs, covparms) ) # cor

    }  else if(abs.corr == TRUE) {

      corr.matrix <- abs( .cov2cor( covmodel(locs, covparms) ) ) # |cor|

    } else {

      stop("abs.corr must be logical.")
    }

  } else if(is.matrix(covmodel)) {

    if(abs.corr == FALSE) {

      corr.matrix <- .cov2cor( covmodel ) # cor

    }  else if(abs.corr == TRUE) {

      corr.matrix <- abs( .cov2cor(covmodel) )# |cor|

    } else {

      stop("abs.corr must be logical.")
    }

  } else {

    stop("covmodel must be either function or matrix.")
  }

  return(corr.matrix)
}

####################################################################################

#' @title specify a general vecchia approximation using the Schaefer's algorithm with known parameters
#'
#' @param locs A matrix of locations
#' @param rho Hyperparameter rho (instead of m)
#' @param initial.pt NULL = euclidean-based center, integer = specify the first obs
#' @param covmodel Covariance model
#' @param covparms Covariance parameters as a vector. The first entry must be its overall variance
#' @param abs.corr Logical. If \code{TRUE} then distance = 1-|rho|. If \code{FALSE} then distane = 1-rho
#' @param coordinate integer or vector of integers in \code{1,...,d}. This argument is used only when ordering is coord
#' @param ordering 'MM' or 'coord'
#' @param ordering.method 'euc' or 'cor'
#' @param conditioning 'NN' (nearest neighbor)
#' @param conditioning.method 'euc' or 'cor'
#'
#' @return An object that specifies the vecchia approximation for later use in likelihood evaluation or prediction
#'
#' @export
#'
#' @examples
#' n             <- 20^2
#' m             <- 20
#' locs          <- matrix(runif(n * 2, 0, 1), n, 2)
#' covparms      <- c(1, 0.1, 10)
#'
#' # true cov matrix
#' covmat <- cov_expo_aniso(locs, covparms)
#'
#' # Visualize the process
#' y <- as.numeric(t(chol(covmat)) %*% rnorm(n))
#' fields::quilt.plot(locs[,1], locs[,2], y)
#'
#' out01 <- cvecchia_rho_specify(locs = locs, rho = 2.5, initial.pt = NULL,
#'                               covmodel = cov_expo_aniso, covparms = covparms,
#'                               abs.corr = FALSE, coordinate = NULL,
#'                               ordering = "MM", ordering.method = "euc",
#'                               conditioning = "NN", conditioning.method = "euc")
#'
#' out02 <- cvecchia_rho_specify(locs = locs, rho = 2.5, initial.pt = NULL,
#'                               covmodel = cov_expo_aniso, covparms = covparms,
#'                               abs.corr = FALSE, coordinate = NULL,
#'                               ordering = "MM", ordering.method = "cor",
#'                               conditioning = "NN", conditioning.method = "cor")
#'
#' kls.euc.euc    <- performance(vecchia.approx = out01, locs = locs,
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#' kls.cor.cor    <- performance(vecchia.approx = out02, locs = locs,
#'                               covmodel = cov_expo_aniso, covparms = covparms)
#'
#' barplot(log10(c(kls.euc.euc,kls.cor.cor)),
#'         names.arg = c("E-MM + E-NN", "C-MM + C-NN"),
#'         main = "Vecchia Approximations", ylab = "log10-scale KL divergence")
cvecchia_rho_specify <- function(locs, rho, initial.pt = NULL, covmodel, covparms, abs.corr = FALSE, coordinate = NULL, ordering = "MM", ordering.method = "cor", conditioning = "NN", conditioning.method = "cor")
{
  ### checkargs
  if(abs.corr) stop("If you want to use rho-based Vecchia, you must use abs.corr = FALSE for now!")

  if(ordering != "MM") stop("If you want to use rho-based Vecchia, you must use MM ordering for now!")

  if(conditioning != "NN") stop("If you want to use rho-based Vecchia, you must use NN conditioning for now!")

  if(ordering.method != conditioning.method) stop("If you want to use rho-based Vecchia, you must be the same method for ordering and conditioning for now!")

  if(is.matrix(covmodel)) stop("If you want to use rho-based Vecchia, you must specify your covariance model for now!")

  if(is.function(covmodel)) covmodel <- deparse(substitute(covmodel))

  if(covmodel %in% c("cov_expo_iso", "cov_expo_aniso", "cov_expo_spacetime", "cov_expo_squared", "cov_matern_iso", "cov_matern_aniso", "cov_matern_scaledim", "cov_matern_spacetime")) {

    # desirable!

  } else if(covmodel %in% c("cov_expo_iso_cpp", "cov_expo_aniso_cpp", "cov_expo_spacetime_cpp", "cov_expo_squared_cpp", "cov_matern_iso_cpp", "cov_matern_aniso_cpp", "cov_matern_scaledim_cpp", "cov_matern_spacetime_cpp")) {

    # desirable!

  } else {

    stop("This function does not work for your covmodel. Please use one of these: cov_expo_iso, cov_expo_aniso, cov_expo_spacetime, cov_expo_squared, cov_matern_iso, cov_matern_aniso, cov_matern_scaledim, and cov_matern_spacetime.")
  }

  ### initialization
  locs  <- as.matrix(locs)
  p     <- ncol(locs)
  n     <- nrow(locs)

  if(is.null(initial.pt)) {

    cen           <- t(as.matrix(colMeans(locs)))
    initial.pt    <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))

  } else if(is.numeric(initial.pt) & length(initial.pt) == 1) {

    # desirable!

  } else {

    stop("The argument initial.pt must be an integer.")
  }

  ### body
  if( ordering.method %in% c("euc", "euclidean") ) {

    output        <- rhoBased_sortSparse(locs = locs, rho = rho, initial.pt = initial.pt, distype = "euclidean", covmodel = covmodel, covparms = covparms)

    output$sparsity   <- methods::as(output$sparsity, "TsparseMatrix")
    maxsize           <- max(table(output$sparsity@j))
    output$sparsity   <- cbind(rev(output$sparsity@j), rev(output$sparsity@i))

    cond.sets     <- conditioning_rho_Rcpp(output$sparsity[, 1], output$sparsity[, 2], output$P-1, maxsize, locs, "euclidean", covmodel, covparms) + 1

  } else if( ordering.method %in% c("cor", "correlation") ) {

    output        <- rhoBased_sortSparse(locs = locs, rho = rho, initial.pt = initial.pt, distype = "correlation", covmodel = covmodel, covparms = covparms)

    output$sparsity   <- methods::as(output$sparsity, "TsparseMatrix")
    maxsize           <- max(table(output$sparsity@j))
    output$sparsity   <- cbind(rev(output$sparsity@j), rev(output$sparsity@i))

    cond.sets     <- conditioning_rho_Rcpp(output$sparsity[, 1], output$sparsity[, 2], output$P-1, maxsize, locs, "correlation", covmodel, covparms) + 1

  } else {

    stop("Check the arguments ordering.method and conditioning.method! They must be euc (euclidean) or cor (correlation).")
  }

  locsord     <- locs[output$P, , drop = FALSE]

  # return
  Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- GPvecchia:::U_sparsity(locsord, cond.sets, obs, Cond)

  ### return
  vecchia.approx <- list(locsord = locsord, obs = obs, ord = output$P, ord.z = output$P, ord.pred = 'general', U.prep = U.prep, cond.yz = 'false', ordering = 'MM', ordering.method = ordering.method, conditioning = 'NN', conditioning.method = conditioning.method)
  return(vecchia.approx)
}

#' @title Schaefer's algorithm with maxmin ordering implemented in C++
#'
#' @param locs A numeric matrix of locations
#' @param rho rho (instead of m)
#' @param initial.pt initial index
#' @param distype "euclidean" or "correlation"
#' @param covmodel covariance function
#' @param covparms A numeric vector of covariance parameters
#'
#' @return A list. Check that out!
#' @export
#'
#' @examples
#' locs <- matrix(c(0.630680164, 0.785343212, 0.725967403, 0.384746804, 0.066332777,
#'                  0.133667622, 0.684830531, 0.642028648, 0.708973463, 0.156014202,
#'                  0.859438927, 0.369073307, 0.030589286, 0.249472289, 0.766846143,
#'                  0.762733708, 0.020763421, 0.676813798, 0.235655257, 0.278117884,
#'                  0.596838115, 0.61879513, 0.840685027, 0.488689207, 0.030872465,
#'                  0.130888217, 0.124606567, 0.050117557, 0.138916279, 0.615490804,
#'                  0.91239613, 0.313120618, 0.050610625, 0.990313355, 0.283005201,
#'                  0.029293502, 0.970113917, 0.837657897, 0.170631242, 0.038255007,
#'                  0.511757776, 0.174439784, 0.565152799, 0.355351824, 0.204757612,
#'                  0.513890368, 0.517269295, 0.791975927, 0.387010792, 0.776910901,
#'                  0.398383379, 0.516435312, 0.556239369, 0.418524799, 0.46378566,
#'                  0.681082833, 0.29392978, 0.067728284, 0.254978952, 0.337992825),
#'                30, 2)
#' rho <- 2.01
#' initial.pt <- 1
#'
#' out1 <- rhoBased_sortSparse(locs = locs, rho = rho, distype = "correlation",
#'                             covmodel = "cov_expo_iso", covparms = c(1, 0.1))
#' out1$P
#'
#' out2 <- cvecchia_m_specify(locs, m = 5, initial.pt = 1,
#'                            covmodel = cov_expo_iso, covparms = c(1, 0.1))
#' out2$ord
rhoBased_sortSparse <- function(locs, rho, initial.pt = 1, distype = "euclidean", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
{
  locs        <- as.matrix(locs)
  initial.pt  <- initial.pt - 1 # 0-based indexing vs. 1-based indexing

  ##### 0-based indexing #####

  output <- sortSparse_Rcpp(locs, rho, initial.pt, distype, covmodel, covparms)

  output$P <- as.integer(output$P)
  output$revP <- as.integer(output$revP)
  output$colptr <- as.integer(output$colptr)
  output$rowval <- as.integer(output$rowval)

  sparseMat <- methods::as(Matrix::sparseMatrix(i = output$rowval, p = output$colptr, dims = rep(nrow(locs), 2), index1 = FALSE), "TsparseMatrix")
  ind_row <- sparseMat@i
  ind_col <- sparseMat@j

  NNchk <- as.numeric(NNcheck_Rcpp(ind_row, ind_col, output$P, output$distances, locs, rho, distype, covmodel, covparms))

  ##### 1-based indexing #####

  output$P <- as.integer(output$P + 1)
  output$revP <- as.integer(output$revP + 1)

  ind_row <- seq(nrow(locs))[ind_row + 1]
  ind_col <- seq(nrow(locs))[ind_col + 1]

  sparseMat <- Matrix::drop0(Matrix::sparseMatrix(i = ind_col, j = ind_row, x = NNchk, dims = rep(nrow(locs), 2), index1 = TRUE))

  ##### return #####

  return(list(P = output$P, revP = output$revP, sparsity = sparseMat, distances = output$distances))
}

#' @title Schaefer's algorithm with reverse-maxmin ordering implemented in C++
#'
#' @param locs A numeric matrix of locations
#' @param rho rho (instead of m)
#' @param initial.pt initial index
#' @param distype "euclidean" or "correlation"
#' @param covmodel covariance function
#' @param covparms A numeric vector of covariance parameters
#'
#' @return A list. Check that out!
#'
#' @export
#'
#' @examples
#' locs <- matrix(c(0.630680164, 0.785343212, 0.725967403, 0.384746804, 0.066332777,
#'                  0.133667622, 0.684830531, 0.642028648, 0.708973463, 0.156014202,
#'                  0.859438927, 0.369073307, 0.030589286, 0.249472289, 0.766846143,
#'                  0.762733708, 0.020763421, 0.676813798, 0.235655257, 0.278117884,
#'                  0.596838115, 0.61879513, 0.840685027, 0.488689207, 0.030872465,
#'                  0.130888217, 0.124606567, 0.050117557, 0.138916279, 0.615490804,
#'                  0.91239613, 0.313120618, 0.050610625, 0.990313355, 0.283005201,
#'                  0.029293502, 0.970113917, 0.837657897, 0.170631242, 0.038255007,
#'                  0.511757776, 0.174439784, 0.565152799, 0.355351824, 0.204757612,
#'                  0.513890368, 0.517269295, 0.791975927, 0.387010792, 0.776910901,
#'                  0.398383379, 0.516435312, 0.556239369, 0.418524799, 0.46378566,
#'                  0.681082833, 0.29392978, 0.067728284, 0.254978952, 0.337992825),
#'                30, 2)
#' rho <- 2.01
#' initial.pt <- 1
#'
#' out1 <- rhoBased_sortSparse_rev(locs = locs, rho = rho, distype = "correlation",
#'                                 covmodel = "cov_expo_iso", covparms = c(1, 0.1))
#' rev(out1$P)
#'
#' out2 <- cvecchia_m_specify(locs, m = 5, initial.pt = 1,
#'                            covmodel = cov_expo_iso, covparms = c(1, 0.1))
#' out2$ord
rhoBased_sortSparse_rev <- function(locs, rho, initial.pt = 1, distype = "euclidean", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
{
  locs        <- as.matrix(locs)
  initial.pt  <- initial.pt - 1 # 0-based indexing vs. 1-based indexing

  ##### 0-based indexing #####

  output <- sortSparse_Rcpp(locs, rho, initial.pt, distype, covmodel, covparms)

  output$P <- as.integer(output$P)
  output$revP <- as.integer(output$revP)
  output$colptr <- as.integer(output$colptr)
  output$rowval <- as.integer(output$rowval)

  sparseMat <- methods::as(Matrix::sparseMatrix(i = output$rowval, p = output$colptr, dims = rep(nrow(locs), 2), index1 = FALSE), "TsparseMatrix")
  ind_row <- sparseMat@i
  ind_col <- sparseMat@j

  NNchk <- as.numeric(NNcheck_Rcpp(ind_row, ind_col, output$P, output$distances, locs, rho, distype, covmodel, covparms))

  ##### 1-based indexing #####

  output$P <- as.integer(output$P + 1)
  output$revP <- as.integer(output$revP + 1)

  ind_row <- rev(seq(nrow(locs)))[ind_row + 1]
  ind_col <- rev(seq(nrow(locs)))[ind_col + 1]

  sparseMat <- Matrix::drop0(Matrix::sparseMatrix(i = ind_col, j = ind_row, x = NNchk, dims = rep(nrow(locs), 2), index1 = TRUE))

  output$P <- rev(output$P)
  output$revP[output$P] <- seq(nrow(locs))
  output$distances <- rev(output$distances)

  ##### return #####

  return(list(P = output$P, revP = output$revP, sparsity = sparseMat, distances = output$distances))
}
