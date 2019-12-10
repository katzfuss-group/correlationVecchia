cov_flexMatern_bruteforce <- function(locs, sigma.mat, nu.mat, alpha.mat) {
  
  n <- nrow(locs)
  p <- nrow(sigma.mat)
  
  cov.mat <- matrix(NA, n*p, n*p)
  for(r in 1:n) {
    for(s in 1:n) {
      
      block.rs <- matrix(NA, p, p)
      for(i in 1:p) {
        for(j in 1:p) {
          block.rs[i, j] <- c_ij_fields_1(loc1 = locs[r, ], loc2 = locs[s, ], sigma = sigma.mat[i, j], nu = nu.mat[i, j], alpha = alpha.mat[i, j])
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




