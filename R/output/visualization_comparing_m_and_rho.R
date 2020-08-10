####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
####################################################################################

library(correlationVecchia)
library(microbenchmark)

set.seed(08062020)
par(mfrow = c(1, 2))

##########################################################################################################################################################################

rm(list = ls())

#####

# n             <- 30^2
n             <- 20^2
locs          <- matrix(runif(n * 2, 0, 1), n, 2)
covmodel      <- "cov_expo_iso"
covparms      <- c(1, 0.1)

plot(locs, xlab = 'x1', ylab = 'x2', main = 'Randomly selected inputs from the unit square domain')

#####

# m <- 20 ; rho <- 5
m <- 10 ; rho <- 3.8

# vecchia_specify
out.m_based     <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", conditioning.method = "euclidean", covmodel = get(covmodel), covparms = covparms)
out.rho_based   <- fastcorrvecchia_specify_knownCovparms(locs = locs, rho = rho, initial.pt = out.m_based$ord[1], distype = "euclidean", covmodel = covmodel, covparms = covparms)

# sum(sizes)
sum(!is.na(out.m_based$U.prep$revNNarray)) - n
sum(!is.na(out.rho_based$U.prep$revNNarray)) - n

# orderings
sum(out.m_based$ord != out.rho_based$ord) / n # 0 = Orderings are same.

#####

locsord     <- out.rho_based$locsord
ord         <- out.rho_based$ord

condset.full <- matrix(NA, n, n)
for(i in 1:n) condset.full[i, seq(i)] <- rev(seq(i))

condset.m   <- out.m_based$U.prep$revNNarray[, rev(seq(ncol(out.m_based$U.prep$revNNarray)))]

condset.rho <- out.rho_based$U.prep$revNNarray[, rev(seq(ncol(out.rho_based$U.prep$revNNarray)))]

# dim(condset.full)
# dim(condset.m)
# dim(condset.rho)
# 
# condset.full[1:10, 1:10]
# condset.m[1:10, 1:10]
# condset.rho[1:10, 1:10]

#####

varvec.full <- c() ; varvec.m <- c() ; varvec.rho <- c()

covftn <- get(covmodel)

varvec.full[1] <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.full[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.full[i, ]), ], covparms = covparms)
  varvec.full[i]  <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

varvec.m[1] <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.m[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.m[i, ]), ], covparms = covparms)
  varvec.m[i]     <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

varvec.rho <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.rho[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.rho[i, ]), ], covparms = covparms)
  varvec.rho[i]   <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

#####

# plot(varvec.full, type = 'o', ylim = range(c(varvec.full, varvec.m, varvec.rho)), col = 'black', ylab = 'conditional variance', main = 'var(w_i | w_{g(i)})')
# points(varvec.m, type = 'o', col = 'red')
# points(varvec.rho, type = 'o', col = 'blue')

#####

ratio.m       <- varvec.m / varvec.full
ratio.rho     <- varvec.rho / varvec.full

plot(log(ratio.m), type = 'o', ylim = range(c(log(ratio.m), log(ratio.rho))), col = 'red', ylab = 'log(variance ratio)', xlab = 'index', main = 'Exponential covariance function')
points(log(ratio.rho), type = 'o', col = 'blue')
legend('topleft', legend = c('m-based Veccchia', 'rho-based Vecchia'), lty = 1, pch = 1, col = c('red', 'blue'))

#####

performance(out = out.m_based, locs = locs, covmodel = get(covmodel), covparms = covparms)
performance(out = out.rho_based, locs = locs, covmodel = get(covmodel), covparms = covparms)

0.5 * sum(log(ratio.m))
0.5 * sum(log(ratio.rho))

par(mfrow = c(1, 1))
plot(rowSums(!is.na(out.rho_based$U.prep$revNNarray)), ylab = "size of conditioning set", xlab = "index") ; abline(h = m, col = 'red')
par(mfrow = c(1, 2))

##########################################################################################################################################################################

rm(list = ls())

#####

# n             <- 30^2
n             <- 20^2
locs          <- expand.grid(x = seq(from = 0.1, to = 0.9, length.out = as.integer(sqrt(n))), y = seq(from = 0.1, to = 0.9, length.out = as.integer(sqrt(n))))
locs          <- as.matrix(locs) + 0.005 * matrix(runif(n * 2, min = -1, max = 1), n, 2)
covmodel      <- "cov_expo_iso"
covparms      <- c(1, 0.1)

plot(locs, xlab = 'x1', ylab = 'x2', main = 'Grid-ish locations (inputs) from the unit sqaure domain')

#####

# m <- 20 ; rho <- 5
m <- 10 ; rho <- 2.58

# vecchia_specify
out.m_based     <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", conditioning.method = "euclidean", covmodel = get(covmodel), covparms = covparms)
out.rho_based   <- fastcorrvecchia_specify_knownCovparms(locs = locs, rho = rho, initial.pt = out.m_based$ord[1], distype = "euclidean", covmodel = covmodel, covparms = covparms)

# sum(sizes)
sum(!is.na(out.m_based$U.prep$revNNarray)) - n
sum(!is.na(out.rho_based$U.prep$revNNarray)) - n

# orderings
sum(out.m_based$ord != out.rho_based$ord) / n # 0 = Orderings are same.

#####

locsord     <- out.rho_based$locsord
ord         <- out.rho_based$ord

condset.full <- matrix(NA, n, n)
for(i in 1:n) condset.full[i, seq(i)] <- rev(seq(i))

condset.m   <- out.m_based$U.prep$revNNarray[, rev(seq(ncol(out.m_based$U.prep$revNNarray)))]

condset.rho <- out.rho_based$U.prep$revNNarray[, rev(seq(ncol(out.rho_based$U.prep$revNNarray)))]

# dim(condset.full)
# dim(condset.m)
# dim(condset.rho)
# 
# condset.full[1:10, 1:10]
# condset.m[1:10, 1:10]
# condset.rho[1:10, 1:10]

#####

varvec.full <- c() ; varvec.m <- c() ; varvec.rho <- c()

covftn <- get(covmodel)

varvec.full[1] <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.full[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.full[i, ]), ], covparms = covparms)
  varvec.full[i]  <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

varvec.m[1] <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.m[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.m[i, ]), ], covparms = covparms)
  varvec.m[i]     <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

varvec.rho <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.rho[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.rho[i, ]), ], covparms = covparms)
  varvec.rho[i]   <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

#####

# plot(varvec.full, type = 'o', ylim = range(c(varvec.full, varvec.m, varvec.rho)), col = 'black', ylab = 'conditional variance', main = 'var(w_i | w_{g(i)})')
# points(varvec.m, type = 'o', col = 'red')
# points(varvec.rho, type = 'o', col = 'blue')

#####

ratio.m       <- varvec.m / varvec.full
ratio.rho     <- varvec.rho / varvec.full

plot(log(ratio.m), type = 'o', ylim = range(c(log(ratio.m), log(ratio.rho))), col = 'red', ylab = 'log(variance ratio)', xlab = 'index', main = 'Exponential covariance function')
points(log(ratio.rho), type = 'o', col = 'blue')
legend('topleft', legend = c('m-based Veccchia', 'rho-based Vecchia'), lty = 1, pch = 1, col = c('red', 'blue'))

#####

performance(out = out.m_based, locs = locs, covmodel = get(covmodel), covparms = covparms)
performance(out = out.rho_based, locs = locs, covmodel = get(covmodel), covparms = covparms)

0.5 * sum(log(ratio.m))
0.5 * sum(log(ratio.rho))

par(mfrow = c(1, 1))
plot(rowSums(!is.na(out.rho_based$U.prep$revNNarray)), ylab = "size of conditioning set", xlab = "index") ; abline(h = m, col = 'red')
par(mfrow = c(1, 2))

##########################################################################################################################################################################

rm(list = ls())

#####

# n             <- 30^2
n             <- 20^2
locs          <- matrix(runif(n * 2, 0, 1), n, 2)
covmodel      <- "cov_matern_2.5"
covparms      <- c(1, 0.1)

plot(locs, xlab = 'x1', ylab = 'x2', main = 'Randomly selected inputs from the unit square domain')

#####

# m <- 20 ; rho <- 5
m <- 10 ; rho <- 3.8

# vecchia_specify
out.m_based     <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", conditioning.method = "euclidean", covmodel = get(covmodel), covparms = covparms)
out.rho_based   <- fastcorrvecchia_specify_knownCovparms(locs = locs, rho = rho, initial.pt = out.m_based$ord[1], distype = "euclidean", covmodel = covmodel, covparms = covparms)

# sum(sizes)
sum(!is.na(out.m_based$U.prep$revNNarray)) - n
sum(!is.na(out.rho_based$U.prep$revNNarray)) - n

# orderings
sum(out.m_based$ord != out.rho_based$ord) / n # 0 = Orderings are same.

#####

locsord     <- out.rho_based$locsord
ord         <- out.rho_based$ord

condset.full <- matrix(NA, n, n)
for(i in 1:n) condset.full[i, seq(i)] <- rev(seq(i))

condset.m   <- out.m_based$U.prep$revNNarray[, rev(seq(ncol(out.m_based$U.prep$revNNarray)))]

condset.rho <- out.rho_based$U.prep$revNNarray[, rev(seq(ncol(out.rho_based$U.prep$revNNarray)))]

# dim(condset.full)
# dim(condset.m)
# dim(condset.rho)
# 
# condset.full[1:10, 1:10]
# condset.m[1:10, 1:10]
# condset.rho[1:10, 1:10]

#####

varvec.full <- c() ; varvec.m <- c() ; varvec.rho <- c()

covftn <- get(covmodel)

varvec.full[1] <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.full[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.full[i, ]), ], covparms = covparms)
  varvec.full[i]  <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

varvec.m[1] <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.m[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.m[i, ]), ], covparms = covparms)
  varvec.m[i]     <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

varvec.rho <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.rho[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.rho[i, ]), ], covparms = covparms)
  varvec.rho[i]   <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

#####

# plot(varvec.full, type = 'o', ylim = range(c(varvec.full, varvec.m, varvec.rho)), col = 'black', ylab = 'conditional variance', main = 'var(w_i | w_{g(i)})')
# points(varvec.m, type = 'o', col = 'red')
# points(varvec.rho, type = 'o', col = 'blue')

#####

ratio.m       <- varvec.m / varvec.full
ratio.rho     <- varvec.rho / varvec.full

plot(log(ratio.m), type = 'o', ylim = range(c(log(ratio.m), log(ratio.rho))), col = 'red', ylab = 'log(variance ratio)', xlab = 'index', main = 'Matern covariance function with nu = 2.5')
points(log(ratio.rho), type = 'o', col = 'blue')
legend('topleft', legend = c('m-based Veccchia', 'rho-based Vecchia'), lty = 1, pch = 1, col = c('red', 'blue'))

#####

performance(out = out.m_based, locs = locs, covmodel = get(covmodel), covparms = covparms)
performance(out = out.rho_based, locs = locs, covmodel = get(covmodel), covparms = covparms)

0.5 * sum(log(ratio.m))
0.5 * sum(log(ratio.rho))

par(mfrow = c(1, 1))
plot(rowSums(!is.na(out.rho_based$U.prep$revNNarray)), ylab = "size of conditioning set", xlab = "index") ; abline(h = m, col = 'red')
par(mfrow = c(1, 2))

##########################################################################################################################################################################

rm(list = ls())

#####

# n             <- 30^2
n             <- 20^2
locs          <- expand.grid(x = seq(from = 0.1, to = 0.9, length.out = as.integer(sqrt(n))), y = seq(from = 0.1, to = 0.9, length.out = as.integer(sqrt(n))))
locs          <- as.matrix(locs) + 0.005 * matrix(runif(n * 2, min = -1, max = 1), n, 2)
covmodel      <- "cov_matern_2.5"
covparms      <- c(1, 0.1)

plot(locs, xlab = 'x1', ylab = 'x2', main = 'Grid-ish locations (inputs) from the unit sqaure domain')

#####

# m <- 20 ; rho <- 5
m <- 10 ; rho <- 2.58

# vecchia_specify
out.m_based     <- corrvecchia_specify_knownCovparms_2(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", conditioning = "NN", conditioning.method = "euclidean", covmodel = get(covmodel), covparms = covparms)
out.rho_based   <- fastcorrvecchia_specify_knownCovparms(locs = locs, rho = rho, initial.pt = out.m_based$ord[1], distype = "euclidean", covmodel = covmodel, covparms = covparms)

# sum(sizes)
sum(!is.na(out.m_based$U.prep$revNNarray)) - n
sum(!is.na(out.rho_based$U.prep$revNNarray)) - n

# orderings
sum(out.m_based$ord != out.rho_based$ord) / n # 0 = Orderings are same.

#####

locsord     <- out.rho_based$locsord
ord         <- out.rho_based$ord

condset.full <- matrix(NA, n, n)
for(i in 1:n) condset.full[i, seq(i)] <- rev(seq(i))

condset.m   <- out.m_based$U.prep$revNNarray[, rev(seq(ncol(out.m_based$U.prep$revNNarray)))]

condset.rho <- out.rho_based$U.prep$revNNarray[, rev(seq(ncol(out.rho_based$U.prep$revNNarray)))]

# dim(condset.full)
# dim(condset.m)
# dim(condset.rho)
# 
# condset.full[1:10, 1:10]
# condset.m[1:10, 1:10]
# condset.rho[1:10, 1:10]

#####

varvec.full <- c() ; varvec.m <- c() ; varvec.rho <- c()

covftn <- get(covmodel)

varvec.full[1] <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.full[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.full[i, ]), ], covparms = covparms)
  varvec.full[i]  <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

varvec.m[1] <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.m[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.m[i, ]), ], covparms = covparms)
  varvec.m[i]     <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

varvec.rho <- covparms[1]
for(i in 2:n) {
  # covmat          <- covftn(locs = locs[na.omit(condset.rho[i, ]), ], covparms = covparms)
  covmat          <- covftn(locs = locsord[na.omit(condset.rho[i, ]), ], covparms = covparms)
  varvec.rho[i]   <- as.numeric(covmat[1, 1] - covmat[1, -1] %*% solve(covmat[-1, -1]) %*% covmat[-1, 1])
}

#####

# plot(varvec.full, type = 'o', ylim = range(c(varvec.full, varvec.m, varvec.rho)), col = 'black', ylab = 'conditional variance', main = 'var(w_i | w_{g(i)})')
# points(varvec.m, type = 'o', col = 'red')
# points(varvec.rho, type = 'o', col = 'blue')

#####

ratio.m       <- varvec.m / varvec.full
ratio.rho     <- varvec.rho / varvec.full

plot(log(ratio.m), type = 'o', ylim = range(c(log(ratio.m), log(ratio.rho))), col = 'red', ylab = 'log(variance ratio)', xlab = 'index', main = 'Matern covariance function with nu = 2.5')
points(log(ratio.rho), type = 'o', col = 'blue')
legend('topleft', legend = c('m-based Veccchia', 'rho-based Vecchia'), lty = 1, pch = 1, col = c('red', 'blue'))

#####

performance(out = out.m_based, locs = locs, covmodel = get(covmodel), covparms = covparms)
performance(out = out.rho_based, locs = locs, covmodel = get(covmodel), covparms = covparms)

0.5 * sum(log(ratio.m))
0.5 * sum(log(ratio.rho))

par(mfrow = c(1, 1))
plot(rowSums(!is.na(out.rho_based$U.prep$revNNarray)), ylab = "size of conditioning set", xlab = "index") ; abline(h = m, col = 'red')
par(mfrow = c(1, 2))

##########################################################################################################################################################################

par(mfrow = c(1, 1))
