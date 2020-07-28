####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions to implement vecchia approximations using C++.
###
###   Contents:
###
####################################################################################



#' @title Schaefer's algorithm with maxmin ordering implemented in C++ 
#'
#' @param locs location 
#' @param rho rho
#' @param initial.pt initial index
#' @param covmodel covmodel
#' @param covparms covariance parameters
#'
#' @return
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
#' rho <- 1.5
#' initial.pt <- 1
#' 
#' out <- rho_based_sortSparse(locs, rho, covmodel = "cov_expo_iso", 
#'                             covparms = c(1, 0.1))
#' out$P
#' corrvecchia_specify_knownCovparms(locs, m = 5, initial.pt = 1, 
#'                                   covmodel = cov_expo_iso, 
#'                                   covparms = c(1, 0.1))$ord
rho_based_sortSparse <- function(locs, rho, initial.pt = 1, covmodel = "cov_expo_iso", covparms = c(1, 0.1))
{
  locs        <- as.matrix(locs)
  initial.pt  <- initial.pt - 1 # 0-based indexing vs. 1-based indexing
  
  ##### 0-based indexing #####
  
  output <- sortSparse_Rcpp(locs, rho, initial.pt, covmodel, covparms)
  
  output$P <- as.integer(output$P)
  output$revP <- as.integer(output$revP)
  output$colptr <- as.integer(output$colptr)
  output$rowval <- as.integer(output$rowval)
  
  sparseMat <- as(Matrix::sparseMatrix(i = output$rowval, p = output$colptr, dims = rep(nrow(locs), 2), index1 = FALSE), "TsparseMatrix")
  ind_row <- sparseMat@i
  ind_col <- sparseMat@j
  
  NNchk <- as.numeric(NNcheck_Rcpp(ind_row, ind_col, output$P, output$distances, locs, rho, covmodel, covparms))
  
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
#' @param locs location 
#' @param rho rho
#' @param initial.pt initial index
#' @param covmodel covmodel 
#' @param covparms covariance parameters
#'
#' @return
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
#' rho <- 1.5
#' initial.pt <- 1
#' 
#' out <- rho_based_sortSparse_reverse(locs, rho, covmodel = "cov_expo_iso", 
#'                                     covparms = c(1, 0.1))
#' out$P
#' corrvecchia_specify_knownCovparms(locs, m = 5, initial.pt = 1, 
#'                                   covmodel = cov_expo_iso, 
#'                                   covparms = c(1, 0.1))$ord
rho_based_sortSparse_reverse <- function(locs, rho, initial.pt = 1, covmodel = "cov_cpp", covparms = c(1, 0.1, 10)) 
{
  locs        <- as.matrix(locs)
  initial.pt  <- initial.pt - 1 # 0-based indexing vs. 1-based indexing
  
  ##### 0-based indexing #####
  
  output <- sortSparse_Rcpp(locs, rho, initial.pt, covmodel, covparms)
  
  output$P <- as.integer(output$P)
  output$revP <- as.integer(output$revP)
  output$colptr <- as.integer(output$colptr)
  output$rowval <- as.integer(output$rowval)
  
  sparseMat <- as(Matrix::sparseMatrix(i = output$rowval, p = output$colptr, dims = rep(nrow(locs), 2), index1 = FALSE), "TsparseMatrix")
  ind_row <- sparseMat@i
  ind_col <- sparseMat@j
  
  NNchk <- as.numeric(NNcheck_Rcpp(ind_row, ind_col, output$P, output$distances, locs, rho, covmodel, covparms))
  
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



