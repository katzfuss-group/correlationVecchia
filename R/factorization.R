####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes the factorize() function to compute a cholesky factor of a user-defined covariance matrix stably.
###             This script also includes the modify() function to compute modified covariance matrix (Simply speaking, modify() = factorize() * factorize()^t). 
###
###   Contents:
###       FTN - factorize() consisting of .algorithm_modchol(), .algorithm_modldl(), .algorithm_GMW81(), and .algorithm_SE99()
###       FTN - modify()
###
####################################################################################



#' @title Cholesky factorization of  a covariane matrix with respect to a selected correction method for positive definiteness
#'
#' @param covmat covariance matrix
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix.  At \code{FALSE} by default
#' @param method An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'               If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'               If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'               If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'               Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'               Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param return.err Logical at \code{TRUE} by default. If return.err is true, then this function returns a decomposition error of the covariance matrix in the Frobenius norm.
#' @param verbose Logical at \code{FALSE} by default. If verbose is true, then this function prints out all messages.
#' @param ... only for the function nearPD()
#'
#' @return \code{factorize} returns
#'     \itemize{
#'         \item{\code{covfactor}: } Cholesky factor of the modified covariance matrix
#'         \item{\code{method}: } An argument specifying a correction method for the cholesky factorization of a covariance matrix.
#'         \item{\code{tol}: } Used Numerical tolerance 
#'         \item{\code{decomp.err}: } (optional) a decomposition error of the covariance matrix
#'     }
#'     
#' @references Fang, Haw-ren, and Dianne P. O’leary. "Modified Cholesky algorithms: a catalog with new approaches." Mathematical Programming 115.2 (2008): 319-349.
#' 
#' @export
#'
#' @examples
#' # Example 1: trivial non positive definite matrix
#' covmat <- matrix(1, 3, 3)
#' 
#' out1 <- factorize(covmat = covmat, pivot = FALSE, method = "qr", tol = .Machine$double.eps, 
#'                   return.err = TRUE, verbose = TRUE)
#' out2 <- factorize(covmat = covmat, pivot = FALSE, method = "diag", tol = .Machine$double.eps, 
#'                   return.err = TRUE, verbose = TRUE)
#' out3 <- factorize(covmat = covmat, pivot = FALSE, method = "eigen-I", tol = .Machine$double.eps, 
#'                   return.err = TRUE, verbose = TRUE)
#' out4 <- factorize(covmat = covmat, pivot = FALSE, method = "eigen-II", tol = .Machine$double.eps, 
#'                   return.err = TRUE, verbose = TRUE)
#' out5 <- factorize(covmat = covmat, pivot = FALSE, method = "GMW81", tol = .Machine$double.eps, 
#'                   return.err = TRUE, verbose = TRUE)
#' out6 <- factorize(covmat = covmat, pivot = TRUE, method = "GMW81", tol = .Machine$double.eps, 
#'                   return.err = TRUE, verbose = TRUE)
#' out7 <- factorize(covmat = covmat, pivot = TRUE, method = "SE99", tol = .Machine$double.eps, 
#'                   return.err = TRUE, verbose = TRUE)
#' out8 <- factorize(covmat = covmat, pivot = FALSE, method = "nearPD", tol = .Machine$double.eps, 
#'                   return.err = TRUE, verbose = TRUE)
#' out9 <- factorize(covmat = covmat, pivot = TRUE, method = "nearPD", tol = .Machine$double.eps, 
#'                   return.err = TRUE, verbose = TRUE)
#' 
#' err <- c(out1$decomp.err, out2$decomp.err, out3$decomp.err, out4$decomp.err, out5$decomp.err, 
#'          out6$decomp.err, out7$decomp.err, out8$decomp.err, out9$decomp.err)
#' err.table <- data.frame(order(err * 1e6), err) ; colnames(err.table) <- c("rank", "decomp.err")
#' err.table
#' 
#' # Example 2: trivial positive definite matrix 
#' covmat <- matrix(c(1, 0.9, 0.7, 0.9, 1, 0.4, 0.7, 0.4, 1), 3, 3)
#' 
#' out01 <- factorize(covmat = covmat, pivot = FALSE, method = NULL, tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out02 <- factorize(covmat = covmat, pivot = TRUE, method = NULL, tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out03 <- factorize(covmat = covmat, pivot = FALSE, method = "qr", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out04 <- factorize(covmat = covmat, pivot = TRUE, method = "qr", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out05 <- factorize(covmat = covmat, pivot = FALSE, method = "diag", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out06 <- factorize(covmat = covmat, pivot = TRUE, method = "diag", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out07 <- factorize(covmat = covmat, pivot = FALSE, method = "type-I", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out08 <- factorize(covmat = covmat, pivot = TRUE, method = "type-I", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out09 <- factorize(covmat = covmat, pivot = FALSE, method = "type-II", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out10 <- factorize(covmat = covmat, pivot = TRUE, method = "type-II", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out11 <- factorize(covmat = covmat, pivot = FALSE, method = "eigen-I", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out12 <- factorize(covmat = covmat, pivot = TRUE, method = "eigen-I", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out13 <- factorize(covmat = covmat, pivot = FALSE, method = "eigen-II", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out14 <- factorize(covmat = covmat, pivot = TRUE, method = "eigen-II", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out15 <- factorize(covmat = covmat, pivot = FALSE, method = "GMW81", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out16 <- factorize(covmat = covmat, pivot = TRUE, method = "GMW81", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out17 <- factorize(covmat = covmat, pivot = FALSE, method = "SE99", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out18 <- factorize(covmat = covmat, pivot = TRUE, method = "SE99", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out19 <- factorize(covmat = covmat, pivot = FALSE, method = "nearPD", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' out20 <- factorize(covmat = covmat, pivot = TRUE, method = "nearPD", tol = .Machine$double.eps, 
#'                    return.err = TRUE, verbose = TRUE)
#' 
#' err <- c(out01$decomp.err, out02$decomp.err, out03$decomp.err, out04$decomp.err, out05$decomp.err, 
#'          out06$decomp.err, out07$decomp.err, out08$decomp.err, out09$decomp.err, out10$decomp.err, 
#'          out11$decomp.err, out12$decomp.err, out13$decomp.err, out14$decomp.err, out15$decomp.err, 
#'          out16$decomp.err, out17$decomp.err, out18$decomp.err, out19$decomp.err, out20$decomp.err)
#' err.table <- data.frame(order(err * 1e6), err) ; colnames(err.table) <- c("rank", "decomp.err")
#' err.table
#' 
#' # Example 3: non-trivial non positive definite matrix 
#' # (from the description of the function nearPD())
#' 
#' covmat <- matrix(c(1,     0.477, 0.644, 0.478, 0.651, 0.826,
#'                    0.477, 1,     0.516, 0.233, 0.682, 0.75,
#'                    0.644, 0.516, 1,     0.599, 0.581, 0.742,
#'                    0.478, 0.233, 0.599, 1,     0.741, 0.8,
#'                    0.651, 0.682, 0.581, 0.741, 1,     0.798,
#'                    0.826, 0.75,  0.742, 0.8,   0.798, 1),
#'                  nrow = 6, ncol = 6)
#' 
#' out001 <- factorize(covmat = covmat, pivot = FALSE, method = "eigen-I", tol = .Machine$double.eps, 
#'                     return.err = TRUE, verbose = TRUE)
#' out002 <- factorize(covmat = covmat, pivot = FALSE, method = "eigen-II", tol = .Machine$double.eps, 
#'                     return.err = TRUE, verbose = TRUE)
#' out003 <- factorize(covmat = covmat, pivot = FALSE, method = "nearPD", tol = .Machine$double.eps, 
#'                     return.err = TRUE, verbose = TRUE)
#' out004 <- factorize(covmat = covmat, pivot = TRUE, method = "nearPD", tol = .Machine$double.eps, 
#'                     return.err = TRUE, verbose = TRUE)
#' 
#' err <- c(out001$decomp.err, out002$decomp.err, out003$decomp.err, out004$decomp.err)
#' err.table <- data.frame(order(err * 1e6), err) ; colnames(err.table) <- c("rank", "decomp.err")
#' err.table
factorize <- function(covmat, pivot = FALSE, method = NULL, tol = .Machine$double.eps, return.err = TRUE, verbose = FALSE, ...) {
  
  n <- nrow(covmat)
  
  # correction
  if( is.null(method) ) { # qr()
    
    covfactor             <- tryCatch(chol(x = covmat, pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the covariance matrix. Please use one of correction methods.")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else if( method == "qr" ) { # eigen()
    
    if(verbose == TRUE & pivot == TRUE) message("This method does not use pivoting.")
    
    eigendecomp           <- eigen(covmat)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors), tol = tol )), error = function(e) "Oops")
    
    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use another correction method.")
    
  } else if( method == "diag" ) { # factorize covmat + diag(tol, n) with chol()
    
    covfactor             <- tryCatch(chol(covmat + diag(tol, n), pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else if( method == "type-I") { # replace diag(covmat) with pmax(abs(diag(covmat)), tol) -> chol()
    
    diag(covmat)         <- pmax(abs(diag(covmat)), tol)
    covfactor             <- tryCatch(chol(covmat, pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else if( method == "type-II") { # replace diag(covmat) with pmax(diag(covmat), tol) -> chol()
    
    diag(covmat)         <- pmax(diag(covmat), tol)
    covfactor             <- tryCatch(chol(covmat, pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else if( method == "eigen-I") { # replace eigen values with pmax(abs(eigen values), tol) -> qr()
    
    if(verbose == TRUE & pivot == TRUE) message("This method does not use pivoting.")
    
    eigendecomp           <- eigen(covmat)
    eigendecomp$values    <- pmax(abs(eigendecomp$values), tol) # ifelse(abs(eigendecomp$values) > tol, abs(eigendecomp$values), tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")
    
    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")
    
  } else if( method == "eigen-II") { # replace eigen values with pmax(eigen values, tol) -> qr()
    
    if(verbose == TRUE & pivot == TRUE) message("This method does not use pivoting.")
    
    eigendecomp           <- eigen(covmat)
    eigendecomp$values    <- pmax(eigendecomp$values, tol) # ifelse(eigendecomp$values > tol, eigendecomp$values, tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")
    
    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")
    
  } else if( method == "GMW81" ) { # please refer .algorithm_GMW81()
    
    out <- tryCatch(.algorithm_GMW81(A = covmat, beta = NULL, pivot = pivot, tol = tol), error = function(e) "Oops")
    
    if(identical(out, "Oops")) stop("GMW81 algorithm cannnot be applied. Please use larger tolerance value (tol).")
    
    if(pivot == FALSE) {
      covfactor <- sqrt(out$Dvec) * t(out$Lmat) # U = t(out$Lmat %*% diag(sqrt(out$Dvec)))
    } else if(pivot == TRUE) {
      covfactor <- sqrt(out$Dvec) * t(out$Lmat)
      covfactor <- covfactor[, out$Pvec]
    } else {
      stop("The argument pivot must be logical: TRUE or FALSE.")
    }
    
  } else if( method == "SE99" ) { # please refer .algorithm_SE99()
    
    if(verbose == TRUE & pivot == FALSE) message("SE99 algorithm must use pivoting.")
    
    out <- tryCatch(.algorithm_SE99(A = covmat, tol = tol), error = function(e) "Oops")
    
    if(identical(out, "Oops")) stop("SE99 algorithm cannnot be applied.")
    
    covfactor <- t(out$Lmat)
    covfactor <- covfactor[, out$Pvec]
    
  } else if( method == "nearPD") {
    
    if(verbose == TRUE) warning("This method only use covmat, pivot, and ellipsis arguments. Plese refer the description for the function nearPD() in the R package Matrix.")
    
    out <- tryCatch(Matrix::nearPD(x = covmat, ...), error = function(e) "Oops")
    
    if(identical(out, "Oops")) stop("nearPD() function cannnot be applied.")
    
    covfactor             <- tryCatch(chol(x = as.matrix(out$mat), pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the covariance matrix. Please use one of correction methods.")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else {
    
    stop("Not yet so far.")
    
  }
  
  # return 
  if(return.err == TRUE) {
    
    covmat.modified   <- t(covfactor) %*% covfactor
    decomp.err        <- sqrt(sum( (covmat - covmat.modified)^2 ))
    
    if(verbose == TRUE) message("The decomposition error of the covariance matrix is ", decomp.err, " in the Frobenius norm.")
    
    return(list(covfactor = covfactor, method = method, tol = tol, decomp.err = decomp.err))
    
  } else {
    
    return(list(covfactor = covfactor, method = method, tol = tol))
    
  }
}

.algorithm_modchol <- function(A, delta, pivot) {
  
  n <- nrow(A)
  
  if(pivot == FALSE) {
    
    P <- NULL
    for(k in 1:(n-1)) {
      A[k, k] <- sqrt(A[k, k] + delta[k])
      A[(k+1):n, k] <- A[(k+1):n, k] / A[k, k]
      A[(k+1):n, (k+1):n] <- A[(k+1):n, (k+1):n] - tcrossprod(A[(k+1):n, k])
    }
    
    A[n, n] <- sqrt(A[n, n] + delta[n])
    A <- A * !upper.tri(A)
    
  } else {
    
    P <- seq(n)
    for(k in 1:(n-1)) {
      pvt             <- which.max(diag(A)[k:n])
      P[k:n]          <- c(P[k:n][pvt], P[k:n][-pvt])
      A               <- A[P, P]
      
      A[k, k] <- sqrt(A[k, k] + delta[k])
      A[(k+1):n, k] <- A[(k+1):n, k] / A[k, k]
      A[(k+1):n, (k+1):n] <- A[(k+1):n, (k+1):n] - tcrossprod(A[(k+1):n, k])
    }
    
    A[n, n] <- sqrt(A[n, n] + delta[n])
    A <- A * !upper.tri(A)
  }
  
  return(list(Lmat = A, Pvec = P))
}

.algorithm_modldl <- function(A, delta, pivot) {
  
  n <- nrow(A)
  
  L <- diag(n)
  D <- rep(NA, n)
  
  if(pivot == FALSE) {
    
    P     <- NULL
    A.k   <- A
    for(k in 1:(n-1)) {
      D[k]            <- A.k[1, 1] + delta[k]
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }
    
    D[n] <- A.k[nrow(A.k), nrow(A.k)] + delta[n]
    
  } else {
    
    P     <- seq(n)
    A.k   <- A
    for(k in 1:(n-1)) {
      pvt             <- which.max(diag(A.k))
      P[k:n]          <- c(P[k:n][pvt], P[k:n][-pvt])
      A.k             <- A.k[c(pvt, (1:nrow(A.k))[-pvt]), c(pvt, (1:nrow(A.k))[-pvt])]
      L               <- L[P, P]
      
      D[k]            <- A.k[1, 1] + delta[k]
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }
    
    D[n] <- A.k[nrow(A.k), nrow(A.k)] + delta[n]
    
  }
  
  return(list(Lmat = L, Dvec = D, Pvec = P))
}

.algorithm_GMW81 <- function(A, beta = NULL, pivot, tol = .Machine$double.eps) {
  
  n <- nrow(A)
  
  if(is.null(beta)) {
    eta <- max(abs( diag(A) ))
    xi <- max(abs( A * (1-diag(n)) ))
    beta <- sqrt(max( eta, xi / sqrt(n^2 - 1), .Machine$double.eps))
  }
  
  L   <- diag(n)
  D   <- rep(NA, n)
  P   <- seq(n)
  
  if(pivot == FALSE) {
    
    P     <- NULL
    A.k   <- A
    for(k in 1:(n-1)) {
      D[k]            <- max(tol, abs(A.k[1, 1]), max(A.k[2:nrow(A.k), 1]^2 / beta^2))
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }
    
    D[n]  <- max(tol, abs(A.k[nrow(A.k), nrow(A.k)]))
    
  } else {
    
    P     <- seq(n)
    A.k   <- A
    for(k in 1:(n-1)) {
      pvt             <- which.max(diag(A.k))
      P[k:n]          <- c(P[k:n][pvt], P[k:n][-pvt])
      A.k             <- A.k[c(pvt, (1:nrow(A.k))[-pvt]), c(pvt, (1:nrow(A.k))[-pvt])]
      L               <- L[P, P]
      
      D[k]            <- max(tol, abs(A.k[1, 1]), max(A.k[2:nrow(A.k), 1]^2 / beta^2))
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
      
    }
    
    D[n]  <- max(tol, abs(A.k[nrow(A.k), nrow(A.k)]))
    
  }
  
  return(list(beta = beta, Lmat = L, Dvec = D, Pvec = P))
}

.algorithm_SE99 <- function(A, tol) { # always use pivoting strategy
  
  outmodchol <- mvLSW:::modchol(A, tol = tol)
  
  return(list(Lmat = outmodchol$L, Pvec = outmodchol$P))
  
}



#' @title Computing a modified covariane matrix with respect to a selected correction method for positive definiteness
#'
#' @param covmat covariance matrix
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix.  At \code{FALSE} by default
#' @param method An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'               If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'               If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'               If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'               Correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#'               Correction method \code{"nearPD"} use a built-in function nearPD() in the R package Matrix.
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param return.err Logical at \code{TRUE} by default. If return.err is true, then this function returns a decomposition error of the covariance matrix in the Frobenius norm.
#' @param verbose Logical at \code{FALSE} by default. If verbose is true, then this function prints out all messages.
#' @param ... only for the function nearPD() 
#'
#' @return \code{factorize} returns
#'     \itemize{
#'         \item{\code{covmat.modified}: } Modified covariance matrix
#'         \item{\code{method}: } An argument specifying a correction method for the cholesky factorization of a covariance matrix.
#'         \item{\code{tol}: } Used Numerical tolerance 
#'         \item{\code{decomp.err}: } (optional) a decomposition error of the covariance matrix
#'     }
#'     
#' @export
#'
#' @examples
#' covmat <- matrix(c(1,     0.477, 0.644, 0.478, 0.651, 0.826,
#'                    0.477, 1,     0.516, 0.233, 0.682, 0.75,
#'                    0.644, 0.516, 1,     0.599, 0.581, 0.742,
#'                    0.478, 0.233, 0.599, 1,     0.741, 0.8,
#'                    0.651, 0.682, 0.581, 0.741, 1,     0.798,
#'                    0.826, 0.75,  0.742, 0.8,   0.798, 1),
#'                  nrow = 6, ncol = 6)
#' 
#' modify(covmat = covmat, pivot = FALSE, method = "eigen-I", tol = .Machine$double.eps, 
#'        return.err = TRUE, verbose = TRUE)
#' modify(covmat = covmat, pivot = FALSE, method = "eigen-II", tol = .Machine$double.eps, 
#'        return.err = TRUE, verbose = TRUE)
#' modify(covmat = covmat, pivot = FALSE, method = "nearPD", tol = .Machine$double.eps, 
#'        return.err = TRUE, verbose = TRUE)
#' modify(covmat = covmat, pivot = TRUE, method = "nearPD", tol = .Machine$double.eps, 
#'        return.err = TRUE, verbose = TRUE)
#' 
modify <- function(covmat, pivot = FALSE, method = NULL, tol = .Machine$double.eps, return.err = TRUE, verbose = FALSE, ...) {
  
  n <- nrow(covmat)
  
  # correction
  if( is.null(method) ) { # qr()
    
    covfactor             <- tryCatch(chol(x = covmat, pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the covariance matrix. Please use one of correction methods.")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else if( method == "qr" ) { # eigen()
    
    if(verbose == TRUE & pivot == TRUE) message("This method does not use pivoting.")
    
    eigendecomp           <- eigen(covmat)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors), tol = tol )), error = function(e) "Oops")
    
    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use another correction method.")
    
  } else if( method == "diag" ) { # factorize covmat + diag(tol, n) with chol()
    
    covfactor             <- tryCatch(chol(covmat + diag(tol, n), pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else if( method == "type-I") { # replace diag(covmat) with pmax(abs(diag(covmat)), tol) -> chol()
    
    diag(covmat)         <- pmax(abs(diag(covmat)), tol)
    covfactor             <- tryCatch(chol(covmat, pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else if( method == "type-II") { # replace diag(covmat) with pmax(diag(covmat), tol) -> chol()
    
    diag(covmat)         <- pmax(diag(covmat), tol)
    covfactor             <- tryCatch(chol(covmat, pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else if( method == "eigen-I") { # replace eigen values with pmax(abs(eigen values), tol) -> qr()
    
    if(verbose == TRUE & pivot == TRUE) message("This method does not use pivoting.")
    
    eigendecomp           <- eigen(covmat)
    eigendecomp$values    <- pmax(abs(eigendecomp$values), tol) # ifelse(abs(eigendecomp$values) > tol, abs(eigendecomp$values), tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")
    
    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")
    
  } else if( method == "eigen-II") { # replace eigen values with pmax(eigen values, tol) -> qr()
    
    if(verbose == TRUE & pivot == TRUE) message("This method does not use pivoting.")
    
    eigendecomp           <- eigen(covmat)
    eigendecomp$values    <- pmax(eigendecomp$values, tol) # ifelse(eigendecomp$values > tol, eigendecomp$values, tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")
    
    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")
    
  } else if( method == "GMW81" ) { # please refer .algorithm_GMW81()
    
    out <- tryCatch(.algorithm_GMW81(A = covmat, beta = NULL, pivot = pivot, tol = tol), error = function(e) "Oops")
    
    if(identical(out, "Oops")) stop("GMW81 algorithm cannnot be applied. Please use larger tolerance value (tol).")
    
    if(pivot == FALSE) {
      covfactor <- sqrt(out$Dvec) * t(out$Lmat) # U = t(out$Lmat %*% diag(sqrt(out$Dvec)))
    } else if(pivot == TRUE) {
      covfactor <- sqrt(out$Dvec) * t(out$Lmat)
      covfactor <- covfactor[, out$Pvec]
    } else {
      stop("The argument pivot must be logical: TRUE or FALSE.")
    }
    
  } else if( method == "SE99" ) { # please refer .algorithm_SE99()
    
    if(verbose == TRUE & pivot == FALSE) message("SE99 algorithm must use pivoting.")
    
    out <- tryCatch(.algorithm_SE99(A = covmat, tol = tol), error = function(e) "Oops")
    
    if(identical(out, "Oops")) stop("SE99 algorithm cannnot be applied.")
    
    covfactor <- t(out$Lmat)
    covfactor <- covfactor[, out$Pvec]
    
  } else if( method == "nearPD") {
    
    if(verbose == TRUE) warning("This method only use covmat, pivot, and ellipsis arguments. Plese refer the description for the function nearPD() in the R package Matrix.")
    
    out <- tryCatch(Matrix::nearPD(x = covmat, ...), error = function(e) "Oops")
    
    if(identical(out, "Oops")) stop("nearPD() function cannnot be applied.")
    
    covfactor             <- tryCatch(chol(x = as.matrix(out$mat), pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the covariance matrix. Please use one of correction methods.")
    
    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }
    
  } else {
    
    stop("Not yet so far.")
    
  }
  
  covmat.modified <- t(covfactor) %*% covfactor

  # return 
  if(return.err == TRUE) {
    
    decomp.err        <- sqrt(sum( (covmat - covmat.modified)^2 ))
    
    if(verbose == TRUE) message("The decomposition error of the covariance matrix is ", decomp.err, " in the Frobenius norm.")
    
    return(list(covmat.modified = covmat.modified, method = method, tol = tol, decomp.err = decomp.err))
    
  } else {
    
    return(list(covmat.modified = covmat.modified, method = method, tol = tol))
    
  }
}
