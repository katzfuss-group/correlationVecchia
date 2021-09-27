####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
###   Contents:
###     cvecchia_specify_noneuc / order_maxmin_correlation_noneuc / find_ordered_cnn_noneuc
###
####################################################################################

#' @title Specify Correlation-based Vecchia approximation without locations
#'
#' @param cormat Correlation matrix
#' @param initial.pt NULL = which.min(rowMeans(d)), integer = specify the first obs, 'random' = at random, and else = which.min(rowMeans(d))
#' @param m Number of nearby points to condition on
#'
#' @return An object that specifies the vecchia approximation for later use in likelihood evaluation or prediction.
#'
#' @export
#'
#' @examples
#' 1 + 1
cvecchia_specify_noneuc <- function(cormat, initial.pt = NULL, m)
{
  ### setting
  n             <- nrow(cormat)

  ### ordering and conditioning
  ord           <- order_maxmin_correlation_noneuc(cormat, initial.pt)
  corord        <- cormat[ord, ord]
  cond.sets     <- find_ordered_cnn_noneuc(corord, m)

  ### return
  Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- GPvecchia:::U_sparsity(matrix(NA, n, 1), cond.sets, obs, Cond)

  vecchia.approx <- list(locsord = matrix(NA, n, 1), obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = "MM", conditioning = 'NN')
  return(vecchia.approx)
}

#' @title Maximum-minimum (maxmin) ordering of locations in terms of the correlation-based distance
#'
#' @param cormat Correlation matrix
#' @param initial.pt NULL = which.min(rowMeans(d)), integer = specify the first obs, 'random' = at random, and else = which.min(rowMeans(d))
#'
#' @return A vector of indices giving the correlation-based maxmin ordering
#'
#' @export
#'
#' @examples
#' 1 + 1
order_maxmin_correlation_noneuc <- function(cormat, initial.pt)
{
  n             <- nrow(cormat)
  ord           <- rep(NA, n)

  # first step
  if( is.null(initial.pt) ){
    ord[1]        <-  which.max(rowSums(cormat))
  } else if( is.numeric(initial.pt) ) {
    ord[1]        <- initial.pt
  } else if( initial.pt == 'random') { # at random
    ord[1]        <- sample(1:n, 1)
  } else {
    ord[1]        <-  which.max(rowSums(cormat))
  }
  cand.argmax   <- seq(n)[seq(n) != ord[1]]

  # second step
  ind           <- as.matrix(expand.grid(cand.argmax, ord[1]))
  cdist         <- cormat[ind]
  ord[2]        <- cand.argmax[which.min(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]

  # third and so on ...
  for(j in 3:(n-1)){
    ind               <- as.matrix(expand.grid(cand.argmax, ord[seq(j-1)]))
    cdist             <- matrix(cormat[ind], nrow = length(cand.argmax), ncol = j-1, byrow = F)
    cdist             <- Rfast::rowMaxs(cdist, value = TRUE)
    ord[j]            <- cand.argmax[which.min(cdist)]
    cand.argmax       <- cand.argmax[cand.argmax != ord[j]]
  }

  ord[n]        <- cand.argmax
  return(ord)
}

#' @title Constructing conditioning sets using only a given ordered correlation matrix
#'
#' @param corord Correlation matrix of ordered inputs
#' @param m Number of nearby points to condition on
#'
#' @return Check \code{conditioning_m_Rcpp()}.
#'
#' @export
#'
#' @examples
#' 1 + 1
find_ordered_cnn_noneuc <- function(corord, m)
{
  return( conditioning_m_Rcpp(m = m, d = 1 - corord) + 1 )
}

#' @title Specify Lexicographic Vecchia approximation without locations
#'
#' @param cormat Correlation matrix
#' @param m Number of nearby points to condition on
#'
#' @return An object that specifies the vecchia approximation for later use in likelihood evaluation or prediction.
#'
#' @export
#'
#' @examples
#' 1 + 1
lvecchia_specify_noneuc <- function(cormat, m)
{
  ### setting
  n             <- nrow(cormat)

  ### ordering and conditioning
  ord           <- 1:n
  cond.sets     <- find_ordered_cnn_noneuc(cormat, m)

  ### return
  Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- GPvecchia:::U_sparsity(matrix(NA, n, 1), cond.sets, obs, Cond)

  vecchia.approx <- list(locsord = matrix(NA, n, 1), obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = 'MM', conditioning = 'NN')
  return(vecchia.approx)
}

#' @title Constructing conditioning sets randomly
#'
#' @param n Number of inputs
#' @param m Number of nearby points to condition on
#'
#' @return Check \code{conditioning_m_Rcpp()}
#'
#' @export
#'
#' @examples
#' 1 + 1
find_ordered_rn <- function(n, m)
{
  cond.sets <- matrix(NA, nrow = n, ncol = m + 1)
  cond.sets[, 1] <- 1:n

  for(i in 2:n) {

    idx <- sample( x = seq(i-1), size = min(m, i-1) )

    if(length(idx) < m) idx <- c(idx, rep(NA, m - length(idx)))

    cond.sets[i, 2:ncol(cond.sets)] <- idx
  }

  return(cond.sets)
}

#' @title Specify Random Vecchia approximation without locations
#'
#' @param cormat Correlation matrix
#' @param m Number of nearby points to condition on
#'
#' @return An object that specifies the vecchia approximation for later use in likelihood evaluation or prediction.
#'
#' @export
#'
#' @examples
#' 1 + 1
rvecchia_specify_noneuc <- function(cormat, m)
{
  ### setting
  n             <- nrow(cormat)

  ### ordering and conditioning
  ord           <- sample(1:n)
  cond.sets     <- find_ordered_rn(n, m)

  ### return
  Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- GPvecchia:::U_sparsity(matrix(NA, n, 1), cond.sets, obs, Cond)

  vecchia.approx <- list(locsord = matrix(NA, n, 1), obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = 'MM', conditioning = 'NN')
  return(vecchia.approx)
}

#' @title Specify Random Vecchia approximations with different values of m but without locations
#'
#' @param cormat Correlation matrix
#' @param ms A numeric vector of numbers of nearby points to condition on
#'
#' @return An list of objects that specify the vecchia approximation for later use in likelihood evaluation or prediction.
#'
#' @export
#'
#' @examples
#' 1 + 1
rvecchias_specify_noneuc <- function(cormat, ms)
{
  ### setting
  n             <- nrow(cormat)

  ### ordering and conditioning
  ord           <- sample(1:n)
  cond.sets.lar <- find_ordered_rn(n, max(ms))

  vecchia.approxs <- list()
  for(i in 1:length(ms)) {

    cond.sets     <- cond.sets.lar[, seq(ms[i]+1)]

    Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
    obs           <- rep(TRUE, n)
    U.prep        <- GPvecchia:::U_sparsity(matrix(NA, n, 1), cond.sets, obs, Cond)

    vecchia.approxs[[i]] <- list(locsord = matrix(NA, n, 1), obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = 'MM', conditioning = 'NN')
  }

  return(vecchia.approxs)
}
