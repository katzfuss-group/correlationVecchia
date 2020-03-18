####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for implementing covariance functions.
###
###   Contents:
###       FTN - cov_expo_iso()
###       FTN - cov_expo_aniso()
###
###       FTN - cov_matern_ns_bruteforce()
###
###       FTN - cov_wave() consisting of .cov_wave_dampedsine(), .cov_wave_dampedcosine(), and .cov_wave_besselJ()
###
###       FTN - cov_bivariate_flexMatern() 
###       FTN - cov_bivariate_expo_latDim()
###
###       FTN - cov_multivariate_flexMatern()
###       FTN - cov_multivariate_flexMatern_GK()
###       FTN - cov_multivariate_flexMatern_GK_bruteforce() consisting of .c_ij_bruteforce(), .c_ij_fields_1(), and .c_ij_fields_2()
###
###       FTN - cov_spacetime_expo()
###
####################################################################################

#' @title Isotropic exponential covariance function
#'
#' @description From a location matrix \code{locs} and a vector with covariance parameters \code{covparms}, this function returns an isotropic exponential covariance matrix which is one of the simplest covariance matrices.
#'
#' @param locs A numerical matrix with \code{n} rows and \code{p} columns. Each row of locs gives a point in R^d
#' @param covparms A numerical vector with covariance parameters in the form (variance, range)
#'
#' @section Parametrization: The covariance parameter vector is (variance, range) = \eqn{(\sigma^2 , r)}. The form of the covariance is \deqn{ C(x, y) = \sigma^2 exp( || x - y || / r)} where \eqn{x} and \eqn{y} are locations in R^d.
#'
#' @return A matrix with \code{n} rows and \code{n} columns, with the \code{(i, j)} entry containing the isotropic exponenital covariance between observations \code{locs[i, ]} and \code{locs[j, ]}
#' @export
#'
#' @examples
#' # grid locations
#' cov_expo_iso(locs = as.matrix(expand.grid(c(0.25, 0.75), c(0.25, 0.75))), covparms = c(1, 0.1))
#'
#' # randomly selected locations
#' cov_expo_iso(locs = matrix(runif(8), 4, 2), covparms = c(1, 0.1))
cov_expo_iso <- function(locs, covparms) {
  
  # checkargs: locs
  if(!is.matrix(locs) & !is.data.frame(locs)) stop("The argument locs is neither matrix nor data.frame.")
  if(!is.numeric(locs)) stop("The argument locs is not numeric.")
  
  # checkargs: covparms: covparms
  if(!is.numeric(covparms)) stop("The argument covparms is not numeric.")
  if(length(covparms) != 2) stop("The argument covparms is not compatible. Please check the length.")
  
  return( covparms[1] * exp(-fields::rdist(locs) / covparms[2]) )
}



#' @title Anisotropic exponential covariance function
#' 
#' @description From a location matrix \code{locs} and a vector with covariance parameters \code{covparms}, this function returns an anisotropic exponential covariance matrix.
#' 
#' @param locs A numerical matrix with \code{n} rows and \code{p} columns. Each row of locs gives a point in R^d
#' @param covparms A numerical vector with covariance parameters. It can be of two different forms: One is the 3-dim'l vector (variance, range, degree of anisotropy) and the other is simply 2-dim'l vector (variance, range). In the latter case, a vector of anisotropy must be specified
#' @param a A vector of anisotropy. At \code{NULL} by default
#'
#' @section Parametrization: The covariance parameters are a variance \eqn{\sigma^2}, range \eqn{r}, and degree (vector) of anisotropy \eqn{\alpha}. The general form of the covariance is \deqn{ C(x, y) = \sigma^2 exp( || A ( x - y ) || / r)} where \eqn{x} and \eqn{y} are locations in R^d. If \eqn{\alpha} is a number, then A is a diagonal matrix \eqn{diag( \sqrt \alpha , 1 , ... , 1 )}. On the other hand, if \eqn{\alpha} is a vector, then A is \eqn{diag( \sqrt \alpha_1, \sqrt \alpha_2, ... , \sqrt \alpha_d )}. 
#'
#' @return A matrix with \code{n} rows and \code{n} columns, with the \code{(i, j)} entry containing the anisotropic exponenital covariance between observations \code{locs[i, ]} and \code{locs[j, ]}
#' @export
#'
#' @examples
#' # grid locations
#' cov_expo_aniso(locs = as.matrix(expand.grid(c(0.25, 0.75), c(0.25, 0.75))), 
#'                covparms = c(1, 0.1, 10))
#' cov_expo_aniso(locs = as.matrix(expand.grid(c(0.25, 0.75), c(0.25, 0.75))), 
#'                covparms = c(1, 0.1), a = c(10, 1)) # equal to above
#'
#' # randomly selected locations
#' cov_expo_aniso(locs = matrix(runif(8), 4, 2), covparms = c(1, 0.1, 10))
#' cov_expo_aniso(locs = matrix(runif(8), 4, 2), covparms = c(1, 0.1), a = c(10, 1)) # equal to above
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



#' @title Brute-force nonstationary Matern covariance function
#'
#' @param locs1 A numerical matrix with \code{n} rows and \code{p} columns. Each row of locs1 gives a point of the first set in R^d
#' @param locs2 A numerical matrix of second set of locations assigned \code{NULL} by default. If this is given as \code{NULL} locs1 is used
#' @param sigma A numeric function for spatially varying standard deviation \eqn{\sigma( loc )}
#' @param smoothness A numeric function for spatially varying smoothness \eqn{\nu( loc )}
#' @param kernel A matrix-valued function for spatially varying (local) geometric anisotropy \eqn{\Sigma( loc )} 
#'
#' @return A matrix with \code{n} rows and \code{n} columns, with the \code{(i, j)} entry containing the nonstationary Matern covariance between observations \code{locs[i, ]} and \code{locs[j, ]}
#' 
#' @references 
#' \itemize{
#'   \item (Main) Katzfuss, Matthias. "Bayesian nonstationary spatial modeling for very large datasets." Environmetrics 24.3 (2013): 189-200.
#'   \item Paciorek, Christopher J., and Mark J. Schervish. "Spatial modelling using a new class of nonstationary covariance functions." Environmetrics: The official journal of the International Environmetrics Society 17.5 (2006): 483-506.
#'   \item Stein, Michael L. "Nonstationary spatial covariance functions." Unpublished technical report (2005).
#' }
#'
#' @export
#'
#' @examples
#' ### compare to the isotropic exponential covariance function
#' aniso <- function(loc) {
#'
#'   d         <- length(loc)
#' 
#'   a         <- function(loc) 1 # please use your own a function for the first coordinate
#'   b         <- function(loc) 1 # please use your own a function for the second coordinate
#'   angle     <- function(loc) 0 # please use your own spatially varying rotation angle
#' 
#'   eta       <- angle(loc)
#'   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
#' 
#'   range     <- c(a(loc)^(-2), b(loc)^(-2))
#' 
#'   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#' }
#' 
#' sigma <- function(loc) determinant(aniso(loc), logarithm = FALSE)[[1]][1]^0.25
#' 
#' smoothness <- function(loc) 0.5
#' 
#' locs <- matrix(runif(6), 3, 2)
#' cov_matern_ns_bruteforce(locs1 = locs, locs2 = NULL, 
#'                          sigma = sigma, smoothness = smoothness, kernel = aniso)
#' cov_expo_iso(locs, covparms = c(1, 1)) # equal to above
#' 
#' 
#' ### compare to the anisotropic exponential covariance function
#' aniso <- function(loc) {
#'
#'   d         <- length(loc)
#' 
#'   a         <- function(loc) 10 
#'   b         <- function(loc) 1 
#'   angle     <- function(loc) 0 
#' 
#'   eta       <- angle(loc)
#'   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
#' 
#'   range     <- c(a(loc)^(-2), b(loc)^(-2))
#' 
#'   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#' }
#' 
#' sigma <- function(loc) determinant(aniso(loc), logarithm = FALSE)[[1]][1]^0.25
#' 
#' smoothness <- function(loc) 0.5
#' 
#' locs <- matrix(runif(6), 3, 2)
#' cov_matern_ns_bruteforce(locs1 = locs, locs2 = NULL, 
#'                          sigma = sigma, smoothness = smoothness, kernel = aniso)
#' cov_expo_aniso(locs, covparms = c(1, 1, 10)) # equal to above
#' 
#' 
#' ### compare to the anisotropic exponential covariance function 2
#' aniso <- function(loc) {
#'
#'   d         <- length(loc)
#' 
#'   a         <- function(loc) 1
#'   b         <- function(loc) 10
#'   angle     <- function(loc) pi/2
#' 
#'   eta       <- angle(loc)
#'   rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = d, ncol = d, byrow = TRUE)
#' 
#'   range     <- c(a(loc)^(-2), b(loc)^(-2))
#' 
#'   return( t(rot.mat) %*% diag(range, nrow = d) %*% rot.mat )
#' }
#' 
#' sigma <- function(loc) determinant(aniso(loc), logarithm = FALSE)[[1]][1]^0.25
#' 
#' smoothness <- function(loc) 0.5
#' 
#' locs <- matrix(runif(6), 3, 2)
#' cov_matern_ns_bruteforce(locs1 = locs, locs2 = NULL, 
#'                          sigma = sigma, smoothness = smoothness, kernel = aniso)
#' cov_expo_aniso(locs, covparms = c(1, 1, 10)) # equal to above
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



#' @title Wave covariance function
#'
#' @param locs A numerical matrix with \code{n} rows and \code{p} columns. Each row of locs gives a point of the first set in R^d
#' @param covparms A numerical vector with covariance parameters. Its form relies on the next argument method. If you use the damped sine function then covparms = (\eqn{\sigma^2}, period). If you use the damped cosine function then covparms = (\eqn{\sigma^2}, range, period). If you use the Bessel J function then covparms = (\eqn{\sigma^2}, nu, period)
#' @param method A method (function) to construct a wave model. It can be "Dampedsine," "Dampedcosine," or "BesselJ." At "Dampedsine" by default
#' @param tol A numerical tolerance for the damped sine function. At \code{1e-8} by default
#'
#' @return A matrix with \code{n} rows and \code{n} columns, with the \code{(i, j)} entry containing the wave covariance between observations \code{locs[i, ]} and \code{locs[j, ]}
#' 
#' @references Abrahamsen, Petter. "A review of Gaussian random fields and correlation functions." (1997).
#' 
#' @export
#'
#' @examples
#' n     <- 15^2
#' locs  <- matrix(runif(n * 2, 0, 1), n, 2)
#' 
#' h     <- fields::rdist(locs)
#' ind   <- order(h)
#' hord  <- h[ind]
#' 
#' covparms  <- c(1, 1/10) # covparms = c(sigma2, period) for the damped sine function
#' covh      <- cov_wave(locs = locs, covparms = covparms, method = "Dampedsine") 
#' covhord   <- covh[ind]
#' plot(hord, covhord, type = 'l', ylim = c(-0.5, 1), 
#'      col = 1, lwd = 2, xlab = "distance", ylab = "covariance")
#' 
#' covparms  <- c(1, 1/10, 1/10) # covparms = c(sigma2, range, period) for the damped cosine function
#' covh      <- cov_wave(locs = locs, covparms = covparms, method = "Dampedcosine") 
#' covhord   <- covh[ind]
#' lines(hord, covhord, type = 'l', ylim = c(-0.5, 1), col = 2, lwd = 2)
#' 
#' covparms  <- c(1, 1/10, 1/10) # covparms = c(sigma2, nu, period) for the Bessel J function
#' covh      <- cov_wave(locs = locs, covparms = covparms, method = "BesselJ")
#' covhord   <- covh[ind]
#' lines(hord, covhord, type = 'l', ylim = c(-0.5, 1), col = 3, lwd = 2)
#' 
#' abline(v = 0, h = 0, col = 'gray', lwd = 2)
#' legend("topright", legend = c("Damped Sine", "Damped Cosine", "Bessel J"), 
#'        col = 1:3, lty = 1, lwd = 2)
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



#' @title Flexible Matern covariance function for bivariate Gaussian Processes
#'
#' @param locs At \code{NULL} by default. When used, it must be a list of two location matrices which is list(locs1 = locs1, locs2 = locs2)
#' @param locs1 At \code{NULL} by default. When used, it must be a numerical matrix with \code{n1} rows and \code{d} columns. Each row of locs1 gives a point of the first set in R^d
#' @param locs2 At \code{NULL} by default. When used, it must be a numerical matrix with \code{n2} rows and \code{d} columns. Each row of locs2 gives a point of the second set in R^d
#' @param sigma.mat A numerical 2 by 2 matrix of collocated covariance coefficients
#' @param nu.mat A numerical 2 by 2 matrix of smoothness parameters
#' @param alpha.mat A numerical 2 by 2 matrix of scale parameters
#'
#' @return A flexible Matern covariance matrix with \code{n1 + n2} rows and \code{n1 + n2} columns
#' 
#' @export
#'
#' @examples
#' locs1 <- matrix(runif(8, 0, 1), 4, 2)
#' locs2 <- matrix(runif(6, 0, 1), 3, 2)
#' 
#' covmat1 <- cov_bivariate_flexMatern(locs1 = locs1, locs2 = locs2, 
#'                                     sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                     nu.mat = matrix(0.5, 2, 2), 
#'                                     alpha.mat = matrix(1, 2, 2))
#' 
#' covmat2 <- cov_bivariate_flexMatern(locs = list(locs1 = locs1, locs2 = locs2), 
#'                                     sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2), 
#'                                     nu.mat = matrix(0.5, 2, 2), 
#'                                     alpha.mat = matrix(1, 2, 2))
#' 
#' sqrt(sum((covmat1 - covmat2)^2))
cov_bivariate_flexMatern <- function(locs = NULL, locs1 = NULL, locs2 = NULL, sigma.mat, nu.mat, alpha.mat)
{
  ### checkargs: locs
  if(is.null(locs) & is.null(locs1) & is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")
  if(is.null(locs) & !is.null(locs1) & is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")
  if(is.null(locs) & is.null(locs1) & !is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")
  if(!is.null(locs) & !is.null(locs1) & is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")
  if(!is.null(locs) & is.null(locs1) & !is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")

  if(is.null(locs) & !is.null(locs1) & !is.null(locs2)) locs <- list(locs1 = locs1, locs2 = locs2)
  if(!is.null(locs) & is.null(locs1) & is.null(locs2)) {
    locs1 <- locs[[1]]
    locs2 <- locs[[2]]
  }
  
  if(!is.matrix(locs1) & !is.data.frame(locs1)) stop("The first location object is neither matrix nor data.frame.")
  if(!is.numeric(locs1)) stop("TThe first location object is not numeric.")
  
  if(!is.matrix(locs2) & !is.data.frame(locs2)) stop("The second location object is neither matrix nor data.frame.")
  if(!is.numeric(locs2)) stop("The second location object is not numeric.")
    
  ### checkargs: sigma.mat
  if(!is.matrix(sigma.mat)) stop("The argument sigma.mat must be a matrix.")
  if(!all(dim(sigma.mat) == c(2, 2))) stop("The argument sigma.mat must be a 2 by 2 matrix.")
  if(!is.numeric(sigma.mat)) stop("The argument sigma.mat must be numeric.")
  
  ### checkargs: nu.mat
  if(!is.matrix(nu.mat)) stop("The argument nu.mat must be a matrix.")
  if(!all(dim(nu.mat) == c(2, 2))) stop("The argument nu.mat must be a 2 by 2 matrix.")
  if(!is.numeric(nu.mat)) stop("The argument nu.mat must be numeric.")
  
  ### checkargs: alpha.mat
  if(!is.matrix(alpha.mat)) stop("The argument alpha.mat must be a matrix.")
  if(!all(dim(alpha.mat) == c(2, 2))) stop("The argument alpha.mat must be a 2 by 2 matrix.")
  if(!is.numeric(alpha.mat)) stop("The argument alpha.mat must be numeric.")
  
  n         <- c(nrow(locs1), nrow(locs2))
  p         <- nrow(sigma.mat)
  
  if(p != 2) stop("This functions is only for computing covariance matrices of bivariate processes.")
  
  covmat <- matrix(NA, sum(n), sum(n))
  for(i in 1:p) {
    for(j in 1:p) {
      ind.row   <- seq(from = 1 + n[1] * (i - 1), to = n[1] + n[2] * (i - 1), by = 1)
      ind.col   <- seq(from = 1 + n[1] * (j - 1), to = n[1] + n[2] * (j - 1), by = 1)
      
      covmat[ind.row, ind.col] <- sigma.mat[i, j] * fields::Matern(fields::rdist(x1 = locs[[i]], x2 = locs[[j]]), alpha = alpha.mat[i, j], nu = nu.mat[i, j])
    }
  }
  
  return(covmat)
}



#' @title Isotropic exponontial covariance function for bivariate Gaussian Processes using latent coordinate
#'
#' @param locs At \code{NULL} by default. When used, it must be a list of two location matrices which is list(locs1 = locs1, locs2 = locs2)
#' @param locs1 At \code{NULL} by default. When used, it must be a numerical matrix with \code{n1} rows and \code{d} columns. Each row of locs1 gives a point of the first set in R^d
#' @param locs2 At \code{NULL} by default. When used, it must be a numerical matrix with \code{n2} rows and \code{d} columns. Each row of locs2 gives a point of the second set in R^d
#' @param covparms A numerical vector with covariance parameters. Its form must be (variance, range, distance between two processes in latent coordinate).
#'
#' @return An isotropic exponential covariance matrix with \code{n1 + n2} rows and \code{n1 + n2} columns
#' 
#' @export
#'
#' @examples
#' locs1 <- matrix(runif(6), 3, 2)
#' locs2 <- matrix(runif(8), 4, 2)
#' covparms <- c(1, 0.1)
#' d.latent <- 1
#' 
#' covmat1 <- cov_bivariate_expo_latDim(locs1 = locs1, locs2 = locs2, 
#'                                      covparms = c(covparms, d.latent))
#' covmat2 <- cov_bivariate_expo_latDim(locs = list(locs1 = locs1, locs2 = locs2), 
#'                                      covparms = c(covparms, d.latent))
#' 
#' sqrt(sum((covmat1 - covmat2)^2))
cov_bivariate_expo_latDim <- function(locs = NULL, locs1 = NULL, locs2 = NULL, covparms)
{
  ### checkargs: locs
  if(is.null(locs) & is.null(locs1) & is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")
  if(is.null(locs) & !is.null(locs1) & is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")
  if(is.null(locs) & is.null(locs1) & !is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")
  if(!is.null(locs) & !is.null(locs1) & is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")
  if(!is.null(locs) & is.null(locs1) & !is.null(locs2)) stop("Please specify the location matrix using locs xor locs1&2.")
  
  if(is.null(locs) & !is.null(locs1) & !is.null(locs2)) locs <- list(locs1 = locs1, locs2 = locs2)
  if(!is.null(locs) & is.null(locs1) & is.null(locs2)) {
    locs1 <- locs[[1]]
    locs2 <- locs[[2]]
  }
  
  if(!is.matrix(locs1) & !is.data.frame(locs1)) stop("The first location object is neither matrix nor data.frame.")
  if(!is.numeric(locs1)) stop("TThe first location object is not numeric.")
  
  if(!is.matrix(locs2) & !is.data.frame(locs2)) stop("The second location object is neither matrix nor data.frame.")
  if(!is.numeric(locs2)) stop("The second location object is not numeric.")
  
  ### Calculation
  for(i in 1:length(locs)) {
    locs[[i]] <- cbind( locs[[i]], (i-1) * covparms[3] )
  }
  
  locs.3d <- do.call(rbind, locs)
  
  return( covparms[1] * exp(-fields::rdist(locs.3d) / covparms[2]) )
}



#' @title Flexible Matern covariance function for Multivariate Processes
#'
#' @param locs A matrix of locations or a list of matrices of locations
#' @param sigma.mat A numerical \code{p} by \code{p} matrix of collocated covariance coefficients at \code{NULL} by default
#' @param p A number of processes at \code{NULL} by default. Please use this argument with the next argument rho
#' @param rho A coefficient for defining a matrix of collocated covariance coefficients as an autoregressive covariance matrix. Please use this argument with the previous argument p
#' @param nu.mat A numerical \code{p} by \code{p} matrix of smoothness parameters at \code{NULL} by default
#' @param nu.vec A numeric vector at \code{NULL} by default. If nu.vec is used instead of nu.mat, the \code{(i, j)}-th entry a matrix of smoothness parametes is defined by \code{(nu_i + nu_j)/2}
#' @param alpha.mat A numerical \code{p} by \code{p} matrix of scale parameters
#' @param alpha Numeric at NULL by default. If alpha is used instead of alpha.mat, a matrix of scale parameters is defined by a constant matrix of alpha
#'
#' @return A flexible Matern covariance matrix
#' 
#' @references Bachoc, François, and Reinhard Furrer. "On the smallest eigenvalues of covariance matrices of multivariate spatial processes." Stat 5.1 (2016): 102-107.
#' 
#' @export
#'
#' @examples
#' # Toy example for the univariate isotropic covariance function
#' n <- 5 ; d <- 2 ; p <- 1
#' locs <- matrix(runif(n * d), n, d)
#' 
#' covmat <- cov_multivariate_flexMatern(locs = locs, 
#'                                       sigma.mat = diag(p), 
#'                                       nu.mat = matrix(0.5, p, p), 
#'                                       alpha.mat = matrix(1, p, p))
#' 
#' covmat
#' exp(-fields::rdist(locs))
#' 
#' isSymmetric(covmat)
#' eigen(covmat)$values
#' 
#' 
#' # Example for the bivariate Matern covariance function
#' n1 <- 4 ; n2 <- 3 ; d <- 2 ; p <- 2
#' locs1 <- matrix(runif(n1 * d), n1, d) ; locs2 <- matrix(runif(n2 * d), n2, d)
#' sigma.mat <- diag(p) ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)
#' 
#' covmat <- cov_multivariate_flexMatern(locs = list(locs1 = locs1, locs2 = locs2), 
#'                                       sigma.mat = sigma.mat, 
#'                                       nu.mat = nu.mat, 
#'                                       alpha.mat = alpha.mat)
#' 
#' isSymmetric(covmat)
#' eigen(covmat)$values
#' 
#' covmat2 <- cov_bivariate_flexMatern(locs1 = locs1, locs2 = locs2, 
#'                                     sigma.mat = sigma.mat, 
#'                                     nu.mat = nu.mat, 
#'                                     alpha.mat = alpha.mat)
#'
#' identical(covmat, covmat2)
#' 
#' # Example for comparison between two differnt versions 
#' #      of covariance matrix of multivariate gaussian process
#' 
#' n <- 3 ; d <- 2 ; p <- 2
#' locs <- matrix(runif(n * d), n, d)
#' sigma.mat <- diag(p) ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)
#' 
#' covmat <- cov_multivariate_flexMatern(locs = locs, 
#'                                       sigma.mat = sigma.mat, 
#'                                       nu.mat = nu.mat, 
#'                                       alpha.mat = alpha.mat)
#' covmat_each <- cov_multivariate_flexMatern_GK(locs = locs, 
#'                                               sigma.mat = sigma.mat, 
#'                                               nu.mat = nu.mat, 
#'                                               alpha.mat = alpha.mat)
#' covmat_each_bruteforce <- cov_multivariate_flexMatern_GK_bruteforce(locs = locs, 
#'                                                                     sigma.mat = sigma.mat, 
#'                                                                     nu.mat = nu.mat, 
#'                                                                     alpha.mat = alpha.mat)
#' 
#' covmat_matched <- matrix(NA, n * p, n * p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'     covmat_matched[seq(from = 1 + n * (i - 1), to = n * i, by = 1), 
#'                    seq(from = 1 + n * (j - 1), to = n * j, by = 1)] <- 
#'                    covmat_each[seq(from = i, to = i + p * (n - 1), by = p), 
#'                                seq(from = j, to = j + p * (n - 1), by = p)]
#'   }
#' }
#' 
#' covmat
#' covmat_each
#' covmat_each_bruteforce
#' 
#' identical(covmat_each, covmat_each_bruteforce)
#' identical(covmat, covmat_matched)
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



#' @title Flexible Matern covariance function for Multivariate Processes (Genton and Kleiber, 2015)
#'
#' @param locs A matrix of locations
#' @param sigma.mat A numerical \code{p} by \code{p} matrix of collocated covariance coefficients at \code{NULL} by default
#' @param p A number of processes at \code{NULL} by default. Please use this argument with the next argument rho
#' @param rho A coefficient for defining a matrix of collocated covariance coefficients as an autoregressive covariance matrix. Please use this argument with the previous argument p
#' @param nu.mat A numerical \code{p} by \code{p} matrix of smoothness parameters at \code{NULL} by default
#' @param nu.vec A numeric vector at \code{NULL} by default. If nu.vec is used instead of nu.mat, the \code{(i, j)}-th entry a matrix of smoothness parametes is defined by \code{(nu_i + nu_j)/2}
#' @param alpha.mat A numerical \code{p} by \code{p} matrix of scale parameters
#' @param alpha Numeric at NULL by default. If alpha is used instead of alpha.mat, a matrix of scale parameters is defined by a constant matrix of alpha
#'
#' @return A flexible Matern covariance matrix
#' 
#' @references Genton, Marc G., and William Kleiber. "Cross-covariance functions for multivariate geostatistics." Statistical Science (2015): 147-163.
#' 
#' @export
#'
#' @examples
#' n <- 3 ; d <- 2 ; p <- 2
#' locs <- matrix(runif(n * d), n, d)
#' sigma.mat <- diag(p) ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)
#' 
#' covmat <- cov_multivariate_flexMatern(locs = locs, 
#'                                       sigma.mat = sigma.mat, 
#'                                       nu.mat = nu.mat, 
#'                                       alpha.mat = alpha.mat)
#' covmat_each <- cov_multivariate_flexMatern_GK(locs = locs, 
#'                                               sigma.mat = sigma.mat, 
#'                                               nu.mat = nu.mat, 
#'                                               alpha.mat = alpha.mat)
#' covmat_each_bruteforce <- cov_multivariate_flexMatern_GK_bruteforce(locs = locs, 
#'                                                                     sigma.mat = sigma.mat, 
#'                                                                     nu.mat = nu.mat, 
#'                                                                     alpha.mat = alpha.mat)
#' 
#' covmat_matched <- matrix(NA, n * p, n * p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'     covmat_matched[seq(from = 1 + n * (i - 1), to = n * i, by = 1), 
#'                    seq(from = 1 + n * (j - 1), to = n * j, by = 1)] <- 
#'                    covmat_each[seq(from = i, to = i + p * (n - 1), by = p), 
#'                                seq(from = j, to = j + p * (n - 1), by = p)]
#'   }
#' }
#' 
#' covmat
#' covmat_each
#' covmat_each_bruteforce
#' 
#' identical(covmat_each, covmat_each_bruteforce)
#' identical(covmat, covmat_matched)
cov_multivariate_flexMatern_GK <- function(locs, sigma.mat = NULL, p = NULL, rho = NULL, nu.mat = NULL, nu.vec = NULL, alpha.mat = NULL, alpha = NULL) 
{
  n <- nrow(locs)
  
  # specifying a matrix of collocated covariance coefficients
  if( !is.null(sigma.mat) & !is.null(rho) ) {
    
    stop("Please use only one argument specifying the collocated covariance coefficient sigma.")
    
  } else if( is.null(sigma.mat) & is.null(rho) ) {
    
    stop("Please specify the collocated covariance coefficient sigma.")
    
  } else if( is.null(sigma.mat) & !is.null(rho) ) {
    
    if(is.null(p)) stop("Please specify the number p of processes if you need to use the argument rho.")
    sigma.mat <- rho^fields::rdist(seq(p))
    
  }
  
  p <- nrow(sigma.mat)
  
  # specifying a matrix of smoothness parameters
  if( !is.null(nu.mat) & !is.null(nu.vec) ) {
    
    stop("Please use only one argument specifying the smoothness parameter nu.")
    
  } else if( is.null(nu.mat) & is.null(nu.vec) ) {
    
    stop("Please specify the smoothness parameter nu.")
    
  } else if( is.null(nu.mat) & !is.null(nu.vec) ) {
    
    nu.mat <- (matrix(nu.vec, p, p, byrow = T) + matrix(nu.vec, p, p, byrow = F)) / 2
    
  }
  
  # specifying a matrix of scale parameters
  if( !is.null(alpha.mat) & !is.null(alpha) ) {
    
    stop("Please use only one argument specifying the scale parameter alpha.")
    
  } else if( is.null(alpha.mat) & is.null(alpha) ) {
    
    stop("Please specify the scale parameter alpha.")
    
  } else if( is.null(alpha.mat) & !is.null(alpha) ) {
    
    alpha.mat <- matrix(alpha, p, p)
    
  }
  
  # compute a covariance matrix
  covmat <- matrix(NA, n * p, n * p)
  for(i in 1:p) {
    for(j in 1:p) {
      ind.row   <- seq(from = i, by = p, length.out = n)
      ind.col   <- seq(from = j, by = p, length.out = n)
      
      covmat[ind.row, ind.col] <- sigma.mat[i, j] * fields::Matern(fields::rdist(x1 = locs, x2 = NULL), alpha = alpha.mat[i, j], nu = nu.mat[i, j])
    }
  }
  
  # return
  return(covmat)
}



#' @title Brute-force flexible Matern covariance function for Multivariate Processes (Genton and Kleiber, 2015)
#'
#' @param locs A matrix of locations
#' @param sigma.mat A numerical \code{p} by \code{p} matrix of collocated covariance coefficients at \code{NULL} by default
#' @param p A number of processes at \code{NULL} by default. Please use this argument with the next argument rho
#' @param rho A coefficient for defining a matrix of collocated covariance coefficients as an autoregressive covariance matrix. Please use this argument with the previous argument p
#' @param nu.mat A numerical \code{p} by \code{p} matrix of smoothness parameters at \code{NULL} by default
#' @param nu.vec A numeric vector at \code{NULL} by default. If nu.vec is used instead of nu.mat, the \code{(i, j)}-th entry a matrix of smoothness parametes is defined by \code{(nu_i + nu_j)/2}
#' @param alpha.mat A numerical \code{p} by \code{p} matrix of scale parameters
#' @param alpha Numeric at NULL by default. If alpha is used instead of alpha.mat, a matrix of scale parameters is defined by a constant matrix of alpha
#'
#' @return A flexible Matern covariance matrix
#' 
#' @references Genton, Marc G., and William Kleiber. "Cross-covariance functions for multivariate geostatistics." Statistical Science (2015): 147-163.
#' 
#' @export
#'
#' @examples
#' n <- 3 ; d <- 2 ; p <- 2
#' locs <- matrix(runif(n * d), n, d)
#' sigma.mat <- diag(p) ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)
#' 
#' covmat <- cov_multivariate_flexMatern(locs = locs, 
#'                                       sigma.mat = sigma.mat, 
#'                                       nu.mat = nu.mat, 
#'                                       alpha.mat = alpha.mat)
#' covmat_each <- cov_multivariate_flexMatern_GK(locs = locs, 
#'                                               sigma.mat = sigma.mat, 
#'                                               nu.mat = nu.mat, 
#'                                               alpha.mat = alpha.mat)
#' covmat_each_bruteforce <- cov_multivariate_flexMatern_GK_bruteforce(locs = locs, 
#'                                                                     sigma.mat = sigma.mat, 
#'                                                                     nu.mat = nu.mat, 
#'                                                                     alpha.mat = alpha.mat)
#' 
#' covmat_matched <- matrix(NA, n * p, n * p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'     covmat_matched[seq(from = 1 + n * (i - 1), to = n * i, by = 1), 
#'                    seq(from = 1 + n * (j - 1), to = n * j, by = 1)] <- 
#'                    covmat_each[seq(from = i, to = i + p * (n - 1), by = p), 
#'                                seq(from = j, to = j + p * (n - 1), by = p)]
#'   }
#' }
#' 
#' covmat
#' covmat_each
#' covmat_each_bruteforce
#' 
#' identical(covmat_each, covmat_each_bruteforce)
#' identical(covmat, covmat_matched)
cov_multivariate_flexMatern_GK_bruteforce <- function(locs, sigma.mat = NULL, p = NULL, rho = NULL, nu.mat = NULL, nu.vec = NULL, alpha.mat = NULL, alpha = NULL) 
{
  n <- nrow(locs)
  
  # specifying a matrix of collocated covariance coefficients
  if( !is.null(sigma.mat) & !is.null(rho) ) {
    
    stop("Please use only one argument specifying the collocated covariance coefficient sigma.")
    
  } else if( is.null(sigma.mat) & is.null(rho) ) {
    
    stop("Please specify the collocated covariance coefficient sigma.")
    
  } else if( is.null(sigma.mat) & !is.null(rho) ) {
    
    if(is.null(p)) stop("Please specify the number p of processes if you need to use the argument rho.")
    sigma.mat <- rho^fields::rdist(seq(p))
    
  }
  
  p <- nrow(sigma.mat)
  
  # specifying a matrix of smoothness parameters
  if( !is.null(nu.mat) & !is.null(nu.vec) ) {
    
    stop("Please use only one argument specifying the smoothness parameter nu.")
    
  } else if( is.null(nu.mat) & is.null(nu.vec) ) {
    
    stop("Please specify the smoothness parameter nu.")
    
  } else if( is.null(nu.mat) & !is.null(nu.vec) ) {
    
    nu.mat <- (matrix(nu.vec, p, p, byrow = T) + matrix(nu.vec, p, p, byrow = F)) / 2
    
  }
  
  # specifying a matrix of scale parameters
  if( !is.null(alpha.mat) & !is.null(alpha) ) {
    
    stop("Please use only one argument specifying the scale parameter alpha.")
    
  } else if( is.null(alpha.mat) & is.null(alpha) ) {
    
    stop("Please specify the scale parameter alpha.")
    
  } else if( is.null(alpha.mat) & !is.null(alpha) ) {
    
    alpha.mat <- matrix(alpha, p, p)
    
  }
  
  # compute a covariance matrix
  covmat <- matrix(NA, n * p, n * p)
  for(r in 1:n) {
    for(s in 1:n) {
      
      block.rs <- matrix(NA, p, p)
      for(i in 1:p) {
        for(j in 1:p) {
          block.rs[i, j] <- .c_ij_fields_2(loc1 = locs[r, ], loc2 = locs[s, ], sigma = sigma.mat[i, j], nu = nu.mat[i, j], alpha = alpha.mat[i, j])
        }
      }
      covmat[seq(from = 1 + p * (r - 1), to = p * r, by = 1), seq(from = 1 + p * (s - 1), to = p * s, by = 1)] <- block.rs
      
    }
  }
  
  # return
  return(covmat)
}

.c_ij_bruteforce <- function(loc1, loc2, sigma, nu, alpha) { # returns NaN when h = 0
  
  h <- sqrt(sum((loc1 - loc2)^2))
  
  return(  sigma * (1 / 2^(nu-1) / gamma(nu)) * (alpha * h)^(nu) * besselK(x = alpha * h, nu = nu)  )
}

.c_ij_fields_1 <- function(loc1, loc2, sigma, nu, alpha) { # returns 1 when h = 0
  
  h <- sqrt(sum((loc1 - loc2)^2))
  
  return(  sigma * fields::Matern(h, range = 1/alpha, smoothness = nu)  )
}

.c_ij_fields_2 <- function(loc1, loc2, sigma, nu, alpha) { # is equivalent to the c_ij_fields_1() function
  
  h <- sqrt(sum((loc1 - loc2)^2))
  
  return(  sigma * fields::Matern(h, alpha = alpha, nu = nu)  )
}



#' @title Non-separable spatio-temporal covariance function
#'
#' @param locs A matrix of locations = (1st coordinate of spatial location, 2nd coordinate of spatial location, temporal location)
#' @param covparms A numeric vector of covariance parameters = (sigma^2, kappa, a, c)
#'
#' @return A matrix with \code{n} rows and \code{n} columns, with the \code{(i, j)} entry containing the non-separable spatio-temporal covariance between observations \code{locs[i, ]} and \code{locs[j, ]} described in the paper below
#' 
#' @references Datta, Abhirup, et al. "Nonseparable dynamic nearest neighbor Gaussian process models for large spatio-temporal data with an application to particulate matter analysis." The annals of applied statistics 10.3 (2016): 1286.
#' 
#' @export
#'
#' @examples
#' locs <- matrix(runif(30), 10, 3)
#' locs
#' 
#' covmat <- cov_spacetime_expo(locs = locs, covparms = c(1, 0.75, 50, 25))
#' fields::image.plot(covmat)
#' 
#' covmat <- cov_spacetime_expo(locs = locs, covparms = c(1, 0.75, 500, 2.5))
#' fields::image.plot(covmat)
cov_spacetime_expo <- function(locs, covparms)
{
  if(ncol(locs) != 3) stop("This covariance function only consider 2 dimensional space and 1 dimensional time.")
  
  h <- fields::rdist(x1 = locs[, 1:2, drop = FALSE])
  u <- fields::rdist(x1 = locs[, 3, drop = FALSE])
  
  return( covparms[1] / (covparms[3] * u^2 + 1)^covparms[2] * exp(- covparms[4] * h / (covparms[3] * u^2 + 1)^(covparms[2]/2) ) )
}
