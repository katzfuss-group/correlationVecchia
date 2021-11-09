####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of approximations in various settings.
###
####################################################################################

rm(list = ls())

set.seed(11082021)

##

gc()

memory.size() # memory.limit()

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

### example ########################################################################

# temp      <- generate_example()
# temp$var_last_given_idx1
# temp$var_last_given_idx2

### example ########################################################################

rho       <- 0.25
n         <- 20
varvec    <- runif(n, min = 0.25, max = 4)

cormat    <- matrix(rho, nrow = n, ncol = n) + diag(1 - rho, n)

# set.seed(11012021)
# output$var_last_given_idx1  > output$var_last_given_idx2
# sum(varvec[output$idx1])    < sum(varvec[output$idx2])
# sum(noise[output$idx1])     > sum(noise[output$idx2])

# set.seed(11042021)
# output$var_last_given_idx1  > output$var_last_given_idx2
# sum(varvec[output$idx1])    < sum(varvec[output$idx2])
# sum(noise[output$idx1])     < sum(noise[output$idx2])

# set.seed(11082021)
# output$var_last_given_idx1  > output$var_last_given_idx2
# sum(varvec[output$idx1])    > sum(varvec[output$idx2])
# sum(noise[output$idx1])     > sum(noise[output$idx2])

# set.seed(11032021)
# output$var_last_given_idx1  < output$var_last_given_idx2
# sum(varvec[output$idx1])    > sum(varvec[output$idx2])
# sum(noise[output$idx1])     < sum(noise[output$idx2])

# set.seed(11062021)
# output$var_last_given_idx1  < output$var_last_given_idx2
# sum(varvec[output$idx1])    < sum(varvec[output$idx2])
# sum(noise[output$idx1])     < sum(noise[output$idx2])

# set.seed(11072021)
# output$var_last_given_idx1  < output$var_last_given_idx2
# sum(varvec[output$idx1])    > sum(varvec[output$idx2])
# sum(noise[output$idx1])     > sum(noise[output$idx2])

# If noise = rep(0, n), then output$var_last_given_idx1 = output$var_last_given_idx2

covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- rep(0.0, n)
# noise     <- rep(0.1, n)
# noise     <- runif(n, min = 0.1, max = 1)

all(eigen(covmat)$values > 0)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1
output$var_last_given_idx2

# output$idx1
# output$idx2

sum(varvec[output$idx1])
sum(varvec[output$idx2])

sum(noise[output$idx1])
sum(noise[output$idx2])

##

idx       <- output$idx1
D         <- diag(sqrt(varvec[idx]))
U         <- diag(noise[idx])
# output$var_last_given_idx1
# covmat[n, n, drop = FALSE] - covmat[n, n, drop = FALSE] * cormat[n, idx, drop = FALSE] %*% solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)) %*% cormat[idx, n, drop = FALSE]
# covmat[n, n] - covmat[n, n] * rho^2 * sum(diag( solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)) %*% matrix(1, nrow = length(idx), ncol = length(idx)) ))

nu        <- sum( solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)))
covmat[n, n] * (1 - rho^2 * nu )

output$var_last_given_idx1

##

idx       <- output$idx2
D         <- diag(sqrt(varvec[idx]))
U         <- diag(noise[idx])
# output$var_last_given_idx1
# covmat[n, n, drop = FALSE] - covmat[n, n, drop = FALSE] * cormat[n, idx, drop = FALSE] %*% solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)) %*% cormat[idx, n, drop = FALSE]
# covmat[n, n] - covmat[n, n] * rho^2 * sum(diag( solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)) %*% matrix(1, nrow = length(idx), ncol = length(idx)) ))

nu        <- sum( solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)))
covmat[n, n] * (1 - rho^2 * nu )

output$var_last_given_idx2

### example ########################################################################

n         <- 20
rho       <- 0.1
varvec    <- runif(n, min = 0.25, max = 4)

cormat    <- correlationVecchia::cov_expo_iso(matrix(runif((n-1) * 2), (n-1), 2), c(1, 0.25))
cormat    <- cbind( rbind(cormat, rho ), rho) ; cormat[n, n] <- 1

covmat    <- cormat * sqrt( matrix(varvec, n, n, byrow = TRUE) * matrix(varvec, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
noise     <- rep(0.0, n)
# noise     <- rep(0.1, n)
# noise     <- runif(n, min = 0.1, max = 1)

output    <- generate_example(locs = matrix(runif(n * 2), n, 2), covmat = covmat, idx1 = sort(sample(seq(n-1), len)), idx2 = sort(sample(seq(n-1), len)), noise = noise)
output$var_last_given_idx1
output$var_last_given_idx2

sum(varvec[output$idx1])
sum(varvec[output$idx2])

# Even if noise = rep(0, n), situation is not trivial

##

idx       <- output$idx1
D         <- diag(sqrt(varvec[idx]))
U         <- diag(noise[idx])
# output$var_last_given_idx1
# covmat[n, n, drop = FALSE] - covmat[n, n, drop = FALSE] * cormat[n, idx, drop = FALSE] %*% solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)) %*% cormat[idx, n, drop = FALSE]
# covmat[n, n] - covmat[n, n] * rho^2 * sum(diag( solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)) %*% matrix(1, nrow = length(idx), ncol = length(idx)) ))

nu        <- sum( solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)))
covmat[n, n] * (1 - rho^2 * nu )

output$var_last_given_idx1

##

idx       <- output$idx2
D         <- diag(sqrt(varvec[idx]))
U         <- diag(noise[idx])
# output$var_last_given_idx1
# covmat[n, n, drop = FALSE] - covmat[n, n, drop = FALSE] * cormat[n, idx, drop = FALSE] %*% solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)) %*% cormat[idx, n, drop = FALSE]
# covmat[n, n] - covmat[n, n] * rho^2 * sum(diag( solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)) %*% matrix(1, nrow = length(idx), ncol = length(idx)) ))

nu        <- sum( solve(cormat[idx, idx, drop = FALSE] + solve(D) %*% U %*% solve(D)))
covmat[n, n] * (1 - rho^2 * nu )

output$var_last_given_idx2

### example ########################################################################

n         <- 20
rho       <- 0.1

varvec1   <- runif(n, min = 0.25, max = 2)
varvec2   <- varvec1 + runif(n, min = 0.25, 2)

locs      <- matrix(runif(n * 2), n, 2)
cormat    <- correlationVecchia::cov_expo_iso(matrix(runif((n-1) * 2), (n-1), 2), c(1, 0.25))
cormat    <- cbind( rbind(cormat, rho), rho) ; cormat[n, n] <- 1

covmat1   <- cormat * sqrt( matrix(varvec1, n, n, byrow = TRUE) * matrix(varvec1, n, n, byrow = FALSE) )
covmat2   <- cormat * sqrt( matrix(varvec2, n, n, byrow = TRUE) * matrix(varvec2, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
idx       <- sort(sample(seq(n-1), len))
noise     <- rep(0.0, n)

output1   <- generate_example(locs = locs, covmat = covmat1, idx1 = idx, idx2 = idx, noise = noise)
output1$var_last_given_idx1

output2    <- generate_example(locs = locs, covmat = covmat2, idx1 = idx, idx2 = idx, noise = noise)
output2$var_last_given_idx1

# For a given correlation matrix, larger noise leads to larger variability of prediction.

### example ########################################################################

n         <- 20
rho       <- 0.1

varvec1   <- runif(n, min = 0.25, max = 2)
varvec2   <- varvec1 + runif(n, min = 0.25, 2)

locs      <- matrix(runif(n * 2), n, 2)
cormat    <- correlationVecchia::cov_expo_iso(matrix(runif(n * 2), n, 2), c(1, 0.25))

covmat1   <- cormat * sqrt( matrix(varvec1, n, n, byrow = TRUE) * matrix(varvec1, n, n, byrow = FALSE) )
covmat2   <- cormat * sqrt( matrix(varvec2, n, n, byrow = TRUE) * matrix(varvec2, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
idx       <- sort(sample(seq(n-1), len))
noise     <- rep(0.0, n)

output1   <- generate_example(locs = locs, covmat = covmat1, idx1 = idx, idx2 = idx, noise = noise)
output1$var_last_given_idx1

output2    <- generate_example(locs = locs, covmat = covmat2, idx1 = idx, idx2 = idx, noise = noise)
output2$var_last_given_idx1

# For a given correlation matrix, larger variance leads to larger variability of prediction.

### example ########################################################################

n         <- 20
rho       <- 0.1

varvec1   <- runif(n, min = 0.25, max = 2)

locs      <- matrix(runif(n * 2), n, 2)
cormat    <- correlationVecchia::cov_expo_iso(matrix(runif(n * 2), n, 2), c(1, 0.25))

covmat    <- cormat * sqrt( matrix(varvec1, n, n, byrow = TRUE) * matrix(varvec1, n, n, byrow = FALSE) )

len       <- as.integer(n/2)
idx       <- sort(sample(seq(n-1), len))
noise1    <- rep(0.0, n)
noise2    <- rep(0.2, n)

output1   <- generate_example(locs = locs, covmat = covmat, idx1 = idx, idx2 = idx, noise = noise1)
output1$var_last_given_idx1

output2    <- generate_example(locs = locs, covmat = covmat, idx1 = idx, idx2 = idx, noise = noise2)
output2$var_last_given_idx1

# For a given correlation matrix, larger noise leads to larger variability of prediction.
