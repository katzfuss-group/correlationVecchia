####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of approximations in various settings.
###
####################################################################################

rm(list = ls())

# set.seed(11082021)

### setting ########################################################################

generate_example <- function(locs = NULL, covmat = NULL, idx1 = NULL, idx2 = NULL, noise = NULL)
{
  ### chkargs
  if( is.null(locs) ) locs <- matrix(c(0.75, 0.25, 0.25, 0.75, 0.25, 0.25), nrow = 3, ncol = 2, byrow = TRUE)

  if( is.null(covmat) ) covmat <- matrix(c(4, 0.125, 0.5, 0.125, 0.25, 0.125, 0.5, 0.125, 1), nrow = 3, ncol = 3, byrow = TRUE)

  if( is.null(idx1) ) idx1 <- 1

  if( is.null(idx2) ) idx2 <- 2

  if( is.null(noise) ) noise <- rep(0, nrow(locs))

  ### cormat
  n         <- nrow(covmat)
  diagvec   <- diag(covmat)

  cormat    <- covmat / sqrt( matrix(diagvec, n, n, byrow = TRUE) * matrix(diagvec, n, n, byrow = FALSE) )

  ### condl var
  var_last_given_idx1 <- covmat[n, n, drop = FALSE] - covmat[n, idx1, drop = FALSE] %*% solve(covmat[idx1, idx1, drop = FALSE] + noise[idx1] * diag(rep(1, length(idx1)))) %*% covmat[idx1, n, drop = FALSE] # covmat[3, 3] * (1 - cormat[3, 1]^2)
  var_last_given_idx2 <- covmat[n, n, drop = FALSE] - covmat[n, idx2, drop = FALSE] %*% solve(covmat[idx2, idx2, drop = FALSE] + noise[idx2] * diag(rep(1, length(idx2)))) %*% covmat[idx2, n, drop = FALSE] # covmat[3, 3] * (1 - cormat[3, 2]^2)

  return( list(locs = locs, covmat = covmat, vars = diagvec, cormat = cormat, noise = noise, idx1 = idx1, idx2 = idx2, var_last_given_idx1 = var_last_given_idx1, var_last_given_idx2 = var_last_given_idx2) )
}

### case 0 ########################################################################

rm(list = ls()[!(ls() %in% c("generate_example"))])

temp      <- generate_example()
temp$var_last_given_idx1
temp$var_last_given_idx2

##

temp      <- generate_example(noise = c(0.1, 0.2))
temp$var_last_given_idx1
temp$var_last_given_idx2

##

n         <- 3
rho       <- 0.15
cormat    <- matrix(c(1, rho, rho, rho, 1, rho, rho, rho, 1), nrow = n, ncol = n)
varvec    <- runif(n, min = 0.5, max = 2)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

temp      <- generate_example(covmat = covmat, noise = c(1, 1, 1))
temp$vars
temp$var_last_given_idx1
temp$var_last_given_idx2

covmat[3, 3] * ( 1 - cormat[3, 1]^2 * (covmat[1, 1]) / (covmat[1, 1] + 1) )

### case 1 ########################################################################

rm(list = ls()[!(ls() %in% c("generate_example"))])

rho       <- 0.15
n         <- 20

cormat    <- matrix(rho, nrow = n, ncol = n) + diag(1 - rho, n)

##

rm(list = ls()[!(ls() %in% c("generate_example", "rho", "n", "cormat"))])

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- rep(0.0, n)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1 ; output$var_last_given_idx2

sum(varvec[output$idx1]) ; sum(varvec[output$idx2])

sum(noise[output$idx1]) ; sum(noise[output$idx2])

##

rm(list = ls()[!(ls() %in% c("generate_example", "rho", "n", "cormat"))])

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- rep(0.1, n)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1 ; output$var_last_given_idx2

sum(varvec[output$idx1]) ; sum(varvec[output$idx2])

sum(noise[output$idx1]) ; sum(noise[output$idx2])

##

rm(list = ls()[!(ls() %in% c("generate_example", "rho", "n", "cormat"))])

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- runif(n, min = 0.1, max = 0.4)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1 ; output$var_last_given_idx2

sum(varvec[output$idx1]) ; sum(varvec[output$idx2])

sum(noise[output$idx1]) ; sum(noise[output$idx2])

### case 2 ########################################################################

rm(list = ls()[!(ls() %in% c("generate_example"))])

rho       <- 0.15
n         <- 20

cormat    <- correlationVecchia::cov_expo_iso(matrix(runif((n-1) * 2), (n-1), 2), c(1, 0.25))
cormat    <- cbind( rbind(cormat, rho ), rho) ; cormat[n, n] <- 1

##

rm(list = ls()[!(ls() %in% c("generate_example", "rho", "n", "cormat"))])

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- rep(0.0, n)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1 ; output$var_last_given_idx2

sum(varvec[output$idx1]) ; sum(varvec[output$idx2])

sum(noise[output$idx1]) ; sum(noise[output$idx2])

##

rm(list = ls()[!(ls() %in% c("generate_example", "rho", "n", "cormat"))])

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- rep(0.1, n)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1 ; output$var_last_given_idx2

sum(varvec[output$idx1]) ; sum(varvec[output$idx2])

sum(noise[output$idx1]) ; sum(noise[output$idx2])

##

rm(list = ls()[!(ls() %in% c("generate_example", "rho", "n", "cormat"))])

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- runif(n, min = 0.1, max = 0.4)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1 ; output$var_last_given_idx2

sum(varvec[output$idx1]) ; sum(varvec[output$idx2])

sum(noise[output$idx1]) ; sum(noise[output$idx2])

### case 3 ########################################################################

rm(list = ls()[!(ls() %in% c("generate_example"))])

rho       <- 0.15
n         <- 20

locs      <- matrix(runif(n * 2), n, 2)
cormat    <- correlationVecchia::cov_expo_iso(matrix(runif(n * 2), n, 2), c(1, 0.25))

##

rm(list = ls()[!(ls() %in% c("generate_example", "rho", "n", "cormat"))])

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- rep(0.0, n)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1 ; output$var_last_given_idx2

sum(varvec[output$idx1]) ; sum(varvec[output$idx2])

sum(noise[output$idx1]) ; sum(noise[output$idx2])

##

rm(list = ls()[!(ls() %in% c("generate_example", "rho", "n", "cormat"))])

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- rep(0.1, n)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1 ; output$var_last_given_idx2

sum(varvec[output$idx1]) ; sum(varvec[output$idx2])

sum(noise[output$idx1]) ; sum(noise[output$idx2])

##

rm(list = ls()[!(ls() %in% c("generate_example", "rho", "n", "cormat"))])

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- runif(n, min = 0.1, max = 0.4)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1 ; output$var_last_given_idx2

sum(varvec[output$idx1]) ; sum(varvec[output$idx2])

sum(noise[output$idx1]) ; sum(noise[output$idx2])

### case 4 ########################################################################

rm(list = ls()[!(ls() %in% c("generate_example"))])

rho       <- 0.15
n         <- 20

cormat    <- correlationVecchia::cov_expo_iso(matrix(runif((n-1) * 2), (n-1), 2), c(1, 0.25))
cormat    <- cbind( rbind(cormat, rho ), rho) ; cormat[n, n] <- 1

##

varvec1   <- runif(n, min = 0.25, max = 2)
varvec2   <- varvec1 + runif(n, min = 0.25, 2)

covmat1   <- cormat * sqrt( matrix(varvec1, n, n, byrow = TRUE) * matrix(varvec1, n, n, byrow = FALSE) )
covmat2   <- cormat * sqrt( matrix(varvec2, n, n, byrow = TRUE) * matrix(varvec2, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
idx       <- sort(sample(seq(n-1), len))
noise     <- rep(0.0, n)

output1   <- generate_example(locs = NULL, covmat = covmat1, idx1 = idx, idx2 = idx, noise = noise)
output2    <- generate_example(locs = NULL, covmat = covmat2, idx1 = idx, idx2 = idx, noise = noise)
output1$var_last_given_idx1 ; output2$var_last_given_idx1

### case 5 ########################################################################

rm(list = ls()[!(ls() %in% c("generate_example"))])

rho       <- 0.15
n         <- 20

cormat    <- correlationVecchia::cov_expo_iso(matrix(runif((n-1) * 2), (n-1), 2), c(1, 0.25))
cormat    <- cbind( rbind(cormat, rho ), rho) ; cormat[n, n] <- 1

##

varvec    <- runif(n, min = 0.25, max = 4)
covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
idx       <- sort(sample(seq(n-1), len))
noise1    <- runif(n, min = 0, max = 0.2)
noise2    <- noise1 + runif(n, min = 0, max = 0.2)

output1   <- generate_example(locs = NULL, covmat = covmat, idx1 = idx, idx2 = idx, noise = noise1)
output2    <- generate_example(locs = NULL, covmat = covmat, idx1 = idx, idx2 = idx, noise = noise2)
output1$var_last_given_idx1 ; output2$var_last_given_idx1
