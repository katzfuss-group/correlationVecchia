cov_flexMatern_bruteforce <- function(locs, sigma.mat, nu.mat, alpha.mat) {
  
  n <- nrow(locs)
  p <- nrow(sigma.mat)
  
  cov.mat <- matrix(NA, n * p, n * p)
  for(r in 1:n) {
    for(s in 1:n) {
      
      block.rs <- matrix(NA, p, p)
      for(i in 1:p) {
        for(j in 1:p) {
          block.rs[i, j] <- c_ij_fields_2(loc1 = locs[r, ], loc2 = locs[s, ], sigma = sigma.mat[i, j], nu = nu.mat[i, j], alpha = alpha.mat[i, j])
        }
      }
      cov.mat[seq(from = 1 + p * (r - 1), to = p * r, by = 1), seq(from = 1 + p * (s - 1), to = p * s, by = 1)] <- block.rs
      
    }
  }
  
  return(cov.mat)
}

c_ij_bruteforce <- function(loc1, loc2, sigma, nu, alpha) { # returns NaN when h = 0
  
  h <- sqrt(sum((loc1 - loc2)^2))
  
  return(  sigma * (1 / 2^(nu-1) / gamma(nu)) * (alpha * h)^(nu) * besselK(x = alpha * h, nu = nu)  )
}

c_ij_fields_1 <- function(loc1, loc2, sigma, nu, alpha) { # returns 1 when h = 0
  
  h <- sqrt(sum((loc1 - loc2)^2))
  
  return(  sigma * fields::Matern(h, range = 1/alpha, smoothness = nu)  )
}

c_ij_fields_2 <- function(loc1, loc2, sigma, nu, alpha) { # is equivalent to the c_ij_fields_1() function
  
  h <- sqrt(sum((loc1 - loc2)^2))
  
  return(  sigma * fields::Matern(h, alpha = alpha, nu = nu)  )
}

# locs <- matrix(runif(4), 2, 2)
# c_ij_bruteforce(loc1 = locs[1, ], loc2 = locs[2, ], sigma = 1, nu = 1, alpha = 1)
# c_ij_fields_1(loc1 = locs[1, ], loc2 = locs[2, ], sigma = 1, nu = 1, alpha = 1)
# c_ij_fields_2(loc1 = locs[1, ], loc2 = locs[2, ], sigma = 1, nu = 1, alpha = 1)
# 
# n <- 5 ; d <- 2 ; p <- 1 
# locs <- matrix(runif(n * d), n, d)
# sigma.mat <- diag(p) 
# nu.mat <- matrix(0.5, p, p)
# alpha.mat <- matrix(1, p, p)
# covmat <- cov_flexMatern_bruteforce(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat)
# covmat
# exp(-fields::rdist(locs))
# isSymmetric(covmat)
# matrixcalc::is.positive.definite(covmat)



cov_flexMatern_old <- function(locs, sigma.mat, nu.mat, alpha.mat) {
  
  n <- nrow(locs)
  p <- nrow(sigma.mat)
  
  cov.mat <- matrix(NA, n*p, n*p)
  for(i in 1:p) {
    for(j in 1:p) {
      ind.row     <- seq(from = i, by = p, length.out = n)
      ind.col     <- seq(from = j, by = p, length.out = n)
      
      cov.mat[ind.row, ind.col]   <- c_ij(locs, sigma = sigma.mat[i, j], nu = nu.mat[i, j], alpha = alpha.mat[i, j])
    }
  }
  
  return(cov.mat)
}

c_ij <- function(locs, sigma, nu, alpha) {
  return(  sigma * fields::Matern(fields::rdist(locs), alpha = alpha, nu = nu)  )
}


cov_flexMatern <- function(locs, locs2 = NULL, sigma.mat, nu.vec = NULL, nu.mat = NULL, alpha = NULL, alpha.mat = NULL) {
  
  n         <- nrow(locs)
  m         <- ifelse(is.null(locs2), n, nrow(locs2))
  p         <- nrow(sigma.mat)
  
  if( !is.null(nu.mat) & !is.null(nu.vec) ) {
    stop("Please use only one argument specifying the smoothness parameter nu.")
  } else if( is.null(nu.mat) & is.null(nu.vec) ) {
    stop("Please specify the smoothness parameter nu.")
  } else if( is.null(nu.mat) & !is.null(nu.vec) ) {
    nu.mat <- (matrix(nu.vec, p, p, byrow = T) + matrix(nu.vec, p, p, byrow = F)) / 2
  }
  
  if( !is.null(alpha.mat) & !is.null(alpha) ) {
    stop("Please use only one argument specifying the scale parameter alpha.")
  } else if( is.null(alpha.mat) & is.null(alpha) ) {
    stop("Please specify the scale parameter alpha.")
  } else if( is.null(alpha.mat) & !is.null(alpha) ) {
    alpha.mat <- matrix(alpha, p, p)
  }
  
  cov.mat   <- matrix(NA, n * p, m * p)
  for(i in 1:p) {
    for(j in 1:p) {
      ind.row   <- seq(from = i, by = p, length.out = n)
      ind.col   <- seq(from = j, by = p, length.out = m)
      
      cov.mat[ind.row, ind.col] <- sigma.mat[i, j] * fields::Matern(fields::rdist(x1 = locs, x2 = locs2), alpha = alpha.mat[i, j], nu = nu.mat[i, j])
    }
  }
  
  return(cov.mat)
}


# n <- 30 ; d <- 2 ; p <- 5
# locs <- matrix(runif(n * d), n, d)
# sigma.mat <- diag(p) ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)
# 
# covmat_bruteforce <- cov_flexMatern_bruteforce(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat)
# covmat_old <- cov_flexMatern_old(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat)
# covmat_new <- cov_flexMatern(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat)
# covmat_cross <- cov_flexMatern(locs = locs, locs2 = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat)
# sqrt(sum((covmat_old - covmat_bruteforce)^2))
# sqrt(sum((covmat_new - covmat_bruteforce)^2))
# sqrt(sum((covmat_cross - covmat_bruteforce)^2))
# 
# microbenchmark::microbenchmark(
#   cov_flexMatern_bruteforce(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat),
#   cov_flexMatern_old(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat),
#   cov_flexMatern(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat),
#   cov_flexMatern(locs = locs, locs2 = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat),
#   times = 10
# )
# 
# 
# n1 <- 3 ; n2 <- 4 ; d <- 2 ; p <- 2
# locs1 <- matrix(runif(n1 * d), n1, d) ; locs2 <- matrix(runif(n2 * d), n2, d) ; locs <- rbind(locs1, locs2)
# sigma.mat <- diag(p) ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)
# 
# round(cov_flexMatern(locs = locs, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat), 3)[seq(1, n1*p), seq(n1*p+1, (n1+n2)*p)]
# round(cov_flexMatern(locs = locs1, locs2 = locs2, sigma.mat = sigma.mat, nu.mat = nu.mat, alpha.mat = alpha.mat), 3)
