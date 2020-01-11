####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for implementing covariance functions.
###
####################################################################################



cov_expo_iso <- function(locs, covparms) {
  
  # checkargs: locs
  if(!is.matrix(locs) & !is.data.frame(locs)) stop("The argument locs is neither matrix nor data.frame.")
  if(!is.numeric(locs)) stop("The argument locs is not numeric.")
  
  # checkargs: covparms: covparms
  if(!is.numeric(covparms)) stop("The argument covparms is not numeric.")
  if(length(covparms) != 2) stop("The argument covparms is not compatible. Please check the length.")
  
  covparms[1] * exp(-fields::rdist(locs) / covparms[2])
}



cov_expo_aniso <- function(locs, covparms, a = NULL) {
  
  # checkargs: locs
  if(!is.matrix(locs) & !is.data.frame(locs)) stop("The argument locs is neither matrix nor data.frame.")
  if(!is.numeric(locs)) stop("The argument locs is not numeric.")
  
  # checkargs: covparms
  if(!is.numeric(covparms)) stop("The argument covparms is not numeric.")
  
  if(length(covparms) == 3 & is.null(a)) {
    
    return( covparms[1] * exp(-fields::rdist(cbind(locs[, 1] * covparms[3], locs[, -1])) / covparms[2]) )
    
  } else if(length(covparms) == 2 & is.numeric(a)) {
    
    return( covparms[1] * exp(-fields::rdist( t(a * t(locs)) ) / covparms[2]) )
    
    # # Comparison in time performance
    # vec <- seq(5)
    # mat <- matrix(seq(100000), 20000, 5)
    # identical(t(vec * t(mat)), sweep(mat, 2, vec, FUN="*"), mapply("*", as.data.frame(mat), vec), as.matrix(mat) %*% diag(vec))
    # microbenchmark::microbenchmark(t(vec * t(mat)), sweep(mat, 2, vec, FUN="*"), mapply(`*`, as.data.frame(mat), vec), as.matrix(mat) %*% diag(vec), times = 10)
    
  } else {
    
    stop("The arguments are invalid. Please refer to the description.")
    
  }
}



cov_matern_ns_bruteforce <- function(locs1, locs2 = NULL, sigma, smoothness, kernel) {
  
  # checkargs: locs
  if(!is.matrix(locs1) & !is.data.frame(locs1)) stop("The argument locs1 is neither matrix nor data.frame.")
  if(!is.numeric(locs1)) stop("The argument locs1 is not numeric.")
  
  if(!is.null(locs2)) {
    if(!is.matrix(locs2) & !is.data.frame(locs2)) stop("The argument locs2 is neither matrix nor data.frame.")
    if(!is.numeric(locs2)) stop("The argument locs2 is not numeric.")
  }
  
  # checkargs: sigma, smootheness, kernel
  if(!is.function(sigma)) stop("The argument sigma must be a function.")
  if(!is.function(smoothness)) stop("The argument smoothness must be a function.")
  if(!is.function(kernel)) stop("The argument kernel must be a function.")
  
  # If locs2 = NULL, locs2 = locs1
  if(is.null(locs2)) locs2 = locs1
  
  n1 <- nrow(locs1) ; n2 <- nrow(locs2) ; d <- ncol(locs1)
  
  mat.cov <- matrix(NA, nrow = n1, ncol = n2)
  for(i in 1:n1){
    for(j in 1:n2){
      sigma.ij      <- sigma(locs1[i, ]) * sigma(locs2[j, ])
      kernel.ij     <- ( kernel(locs1[i, ]) + kernel(locs2[j, ]) ) / 2 
      smooth.ij     <- ( smoothness(locs1[i, ]) + smoothness(locs2[j, ]) ) / 2
      
      q.ij          <- as.numeric(crossprod( locs1[i, ] - locs2[j, ], base::solve(kernel.ij, locs1[i, ] - locs2[j, ]) ))
      
      mat.cov[i,j]  <- sigma.ij * fields::Matern( sqrt(q.ij), nu = smooth.ij) / sqrt( base::determinant(kernel.ij, logarithm = FALSE)[[1]][1] )
    }
  }
  
  return(mat.cov)
}



cov_wave <- function(locs, covparms, method = "Dampedsine", tol = 1e-8) {
  
  # checkargs: locs
  if(!is.matrix(locs) & !is.data.frame(locs)) stop("The argument locs is neither matrix nor data.frame.")
  if(!is.numeric(locs)) stop("The argument locs is not numeric.")
  
  # checkargs: covparms
  if(!is.numeric(covparms)) stop("The argument covparms is not numeric.")
  
  # checkargs: tol
  if(!is.numeric(tol)) stop("The argument tol is not numeric.")
  
  # considerargs: method
  if(method == "Dampedsine") {
    
    return(.cov_wave_dampedsine(locs, covparms, tol))
    
  } else if(method == "Dampedcosine") {
    
    return(.cov_wave_dampedcosine(locs, covparms))
    
  } else if(method == "BesselJ") {
    
    return(.cov_wave_besselJ(locs, covparms))
    
  } else {
    
    stop("Select one of the following covariance models: Dampedsine, Dampedcosine, or besselJ.")
    
  }
}

.cov_wave_dampedsine <- function(locs, covparms, tol) { # covparms = c(sigma2, period)
  
  h       <- fields::rdist(locs)
  ind     <- which(h < tol, arr.ind = T)
  
  covmat        <- covparms[1] * sin(h / covparms[2]) * (covparms[2] / h)
  covmat[ind]   <- covparms[1]
  
  return(covmat)
}

.cov_wave_dampedcosine <- function(locs, covparms) { # covparms = c(sigma2, range, period)
  
  d       <- ncol(locs)
  
  if( d == 2 & covparms[2] > covparms[3] ) stop("Covariance model is not positive definite because of the covariance parameters.")
  if( d == 3 & sqrt(3) * covparms[2] > covparms[3] ) stop("Covariance model is not positive definite because of the covariance parameters.")
  
  h       <- fields::rdist(locs)
  
  return( covparms[1] * exp(- h / covparms[2]) * cos(h / covparms[3]) )
}

.cov_wave_besselJ <- function(locs, covparms) { # covparms = c(sigma2, nu, period)
  
  h       <- fields::rdist(locs)
  
  covmat        <- gamma(covparms[2] + 1) * ( 2 * covparms[3] / c(h) )^covparms[2] * Bessel::BesselJ(c(h) / covparms[3], covparms[2])
  covmat[which(is.nan(covmat))] <- 1
  
  return( covparms[1] * matrix(covmat, nrow = nrow(h), ncol = ncol(h)) )
}



cov_multivariate_flexMatern <- function(locs, sigma.mat = NULL, p = NULL, rho = NULL, nu.mat = NULL, nu.vec = NULL, alpha.mat = NULL, alpha = NULL) 
{
  ### check arguments
  
  # checkargs: locs
  if(!( is.matrix(locs) | is.data.frame(locs) | is.list(locs) )) stop("The argument locs must be either a matrix or list of matrices.") 
  if(is.list(locs)) if(!all( unlist(Map(class, locs)) == "matrix" | unlist(Map(class, locs)) == "data.frame" )) stop("The argument is a list but not a list of matrices.") 
  
  # checkargs: sigma.mat, p, rho
  if( !is.null(sigma.mat) & !is.null(rho) & !is.null(p) ) {
    
    stop("Please use the arguments sigma.mat, rho, p as follows: One way is to specify sigma.mat and leave both rho and p set to NULL by default. The other way is to specify both rho and p and leave sigma.mat et to NULL by default.")
    
  } else if( !is.null(sigma.mat) & !is.null(rho) & is.null(p) ) {
    
    stop("Please use the arguments sigma.mat, rho, p as follows: One way is to specify sigma.mat and leave both rho and p set to NULL by default. The other way is to specify both rho and p and leave sigma.mat et to NULL by default.")
    
  } else if( !is.null(sigma.mat) & is.null(rho) & !is.null(p) ) {
    
    stop("Please use the arguments sigma.mat, rho, p as follows: One way is to specify sigma.mat and leave both rho and p set to NULL by default. The other way is to specify both rho and p and leave sigma.mat et to NULL by default.")
    
  } else if( is.null(sigma.mat) & !is.null(rho) & !is.null(p) ) { # use rho and p
    
    sigma.mat <- rho^fields::rdist(seq(p))
    
  } else if( is.null(sigma.mat) & !is.null(rho) & is.null(p) ) {
    
    stop("Please use the arguments sigma.mat, rho, p as follows: One way is to specify sigma.mat and leave both rho and p set to NULL by default. The other way is to specify both rho and p and leave sigma.mat et to NULL by default.")
    
  } else if( is.null(sigma.mat) & is.null(rho) & !is.null(p) ) {
    
    stop("Please use the arguments sigma.mat, rho, p as follows: One way is to specify sigma.mat and leave both rho and p set to NULL by default. The other way is to specify both rho and p and leave sigma.mat et to NULL by default.")
    
  } else if( is.null(sigma.mat) & is.null(rho) & is.null(p) ) {
    
    stop("Please use the arguments sigma.mat, rho, p as follows: One way is to specify sigma.mat and leave both rho and p set to NULL by default. The other way is to specify both rho and p and leave sigma.mat et to NULL by default.")
    
  } # else = use sigma.mat # !is.null(sigma.mat) & is.null(rho) & is.null(p) 
  
  if(nrow(sigma.mat) != ncol(sigma.mat)) stop("sigma.mat must be a square matrix.")
  if(!is.numeric(sigma.mat)) stop("sigma.mat must be a numeric matrix.")
  p <- nrow(sigma.mat)
  
  # checkargs: compatibility of locs and sigma.mat
  if(is.list(locs)) if(length(locs) != p) stop("locs is not compatible with sigma.mat.")
  
  # checkargs: nu.mat
  if( !is.null(nu.mat) & !is.null(nu.vec) ) {
    
    stop("Please use the arguments nu.mat and nu.vec as follows: One way is to specify nu.mat and leave nu.vec set to NULL by default. The other way is to specify nu.vec and leave nu.mat set to NULL by default.")
    
  } else if( is.null(nu.mat) & !is.null(nu.vec) ) {
    
    if(length(nu.vec != p)) stop("nu.vec is not compatible with sigma.mat.")
    nu.mat <- (matrix(nu.vec, p, p, byrow = T) + matrix(nu.vec, p, p, byrow = F)) / 2
    
  } else if( is.null(nu.mat) & is.null(nu.vec) ) {
    
    stop("Please use the arguments nu.mat and nu.vec as follows: One way is to specify nu.mat and leave nu.vec set to NULL by default. The other way is to specify nu.vec and leave nu.mat set to NULL by default.")
    
  } # else = use nu.mat # !is.null(nu.mat) & is.null(nu.vec)
  
  if(!is.numeric(nu.mat)) stop("nu.mat must be a numeric matrix.")
  
  # checkargs: compatibility of sigma.mat and nu.mat
  if(nrow(nu.mat) != p | ncol(nu.mat) != p) stop("nu.mat is not compatible with sigma.mat.")
  
  # checkargs: alpha.mat, alpha
  if( !is.null(alpha.mat) & !is.null(alpha) ) {
    
    stop("Please use the arguments alpha.mat and alpha as follows: One way is to specify alpha.mat and leave alpha set to NULL by default. The other way is to specify alpha and leave alpha.mat set to NULL by default.")
    
  } else if( is.null(alpha.mat) & !is.null(alpha) ) {
    
    alpha.mat <- matrix(alpha, p, p)
    
  } else if( is.null(alpha.mat) & is.null(alpha) ) {
    
    stop("Please use the arguments alpha.mat and alpha as follows: One way is to specify alpha.mat and leave alpha set to NULL by default. The other way is to specify alpha and leave alpha.mat set to NULL by default.")
    
  } # else = use alpha.mat # !is.null(alpha.mat) & is.null(alpha)
  
  if(!is.numeric(alpha.mat)) stop("alpha.mat must be a numeric matrix.")
  
  
  ### compute a covariance matrix
  if(is.matrix(locs)) {
    
    n           <- nrow(locs)
    
    covmat      <- matrix(NA, n * p, n * p)
    for(i in 1:p) {
      for(j in 1:p) {
        ind.row       <- seq(from = 1 + n * (i - 1), to = n * i, by = 1)
        ind.col       <- seq(from = 1 + n * (j - 1), to = n * j, by = 1)
        
        covmat[ind.row, ind.col]  <- sigma.mat[i, j] * fields::Matern(fields::rdist(x1 = locs, x2 = NULL), alpha = alpha.mat[i, j], nu = nu.mat[i, j])
      }
    }
    
  } else if(is.list(locs)) {
    
    n           <- sapply(locs, nrow)
    n.cumsum    <- cumsum(n) 
    
    covmat      <- matrix(NA, n.cumsum[p], n.cumsum[p])
    for(i in 1:p) {
      for(j in 1:p) {
        ind.row       <- seq(from = 1 + ifelse(i == 1, 0, n.cumsum[i-1]), to = n.cumsum[i], by = 1)
        ind.col       <- seq(from = 1 + ifelse(j == 1, 0, n.cumsum[j-1]), to = n.cumsum[j], by = 1)
        
        covmat[ind.row, ind.col]  <- sigma.mat[i, j] * fields::Matern(fields::rdist(x1 = locs[[i]], x2 = locs[[j]]), alpha = alpha.mat[i, j], nu = nu.mat[i, j])
      }
    }
    
  } else {
    
    stop("locs must be either a matrix of locations or a list of matrices of locations.")
    
  }
  
  return(covmat)
}