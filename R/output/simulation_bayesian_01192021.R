####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: 
###
####################################################################################

rm(list = ls())
set.seed(01192020)

### load packages
library(correlationVecchia)
library(foreach)

### ncores ######################################

no_cores  <- parallel::detectCores() - 2

### cov_bivariate_expo_latDim ###################

n               <- 20^2
covparms        <- c(1, 0.1, 0.1)
process         <- generate_gp_space(nsim = 1, n = n, d = 2, p = 2, method.locs = 'random', covmodel = cov_bivariate_expo_latDim, covparms = covparms)

locs            <- process$sim$sim1$locs
locsall         <- rbind(locs$locs1, locs$locs2)

nugget          <- 0.0
z               <- process$sim$sim1$y + nugget * rnorm(2 * n)

ms              <- c(10, 30)
approxs         <- c("b1", "b2", "b3", "b4", "cc")
candid          <- expand.grid(ms, approxs) ; colnames(candid) <- c("m", "approx")

### range - computation

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.bivariate.range <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_bivariate_expo_latDim(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "range", N = 200, xlim = c(0.08, 0.13), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

### range - visualization

par(mfrow = c(5, 2))
for(i in 1:length(out.bivariate.range)) {
  temp <- out.bivariate.range[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

### distance - computation

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.bivariate.distance <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_bivariate_expo_latDim(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "distance", N = 200, xlim = c(0.01, 0.3), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

### distance - visualization

par(mfrow = c(5, 2))
for(i in 1:length(out.bivariate.distance)) {
  temp <- out.bivariate.distance[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

### save

save(n, covparms, locs, nugget, z, candid, out.bivariate.range, out.bivariate.distance, file = "pilotstudy_bayesian_multivariate.RData")

### reset

rm(list = ls())
no_cores  <- parallel::detectCores() - 2

### cov_matern_spacetime ###################

ms              <- c(10, 30)
approxs         <- c("b1", "b2", "b3", "cc")
candid          <- expand.grid(ms, approxs) ; colnames(candid) <- c("m", "approx")

### n = 400, monitoring station

n               <- 20^2
covparms        <- c(1, 0.1, 0.2, 0.5)
process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "space.random.time.grid", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs
nugget          <- 0.00

z               <- process$sim$sim1$y + nugget * rnorm(n) # fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.spatial.n400.monitoring <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.05, 0.15), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.temporal.n400.monitoring <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.01, 0.6), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

par(mfrow = c(4, 2))
for(i in 1:length(out.spatial.n400.monitoring)) {
  temp <- out.spatial.n400.monitoring[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(4, 2))
for(i in 1:length(out.temporal.n400.monitoring)) {
  temp <- out.temporal.n400.monitoring[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

### n = 900, monitoring station

n               <- 30^2
covparms        <- c(1, 0.1, 0.2, 0.5)
process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "space.random.time.grid", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs
nugget          <- 0.00

z               <- process$sim$sim1$y + nugget * rnorm(n) # fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.spatial.n900.monitoring <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.05, 0.15), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.temporal.n900.monitoring <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.01, 0.6), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

par(mfrow = c(4, 2))
for(i in 1:length(out.spatial.n900.monitoring)) {
  temp <- out.spatial.n900.monitoring[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(4, 2))
for(i in 1:length(out.temporal.n900.monitoring)) {
  temp <- out.temporal.n900.monitoring[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

### n = 400, satellite

n               <- 20^2
covparms        <- c(1, 0.1, 0.2, 0.5)
process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs
nugget          <- 0.00

z               <- process$sim$sim1$y + nugget * rnorm(n) # fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.spatial.n400.satellite <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.05, 0.15), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.temporal.n400.satellite <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.01, 0.6), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

par(mfrow = c(4, 2))
for(i in 1:length(out.spatial.n400.satellite)) {
  temp <- out.spatial.n400.satellite[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(4, 2))
for(i in 1:length(out.temporal.n400.satellite)) {
  temp <- out.temporal.n400.satellite[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

### n = 900, satellite

n               <- 30^2
covparms        <- c(1, 0.1, 0.2, 0.5)
process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs
nugget          <- 0.00

z               <- process$sim$sim1$y + nugget * rnorm(n) # fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.spatial.n900.satellite <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.05, 0.15), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

cl        <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

out.temporal.n900.satellite <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = 0, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.01, 0.6), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time()

par(mfrow = c(4, 2))
for(i in 1:length(out.spatial.n900.satellite)) {
  temp <- out.spatial.n900.satellite[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(4, 2))
for(i in 1:length(out.temporal.n900.satellite)) {
  temp <- out.temporal.n900.satellite[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density')
} ; rm(temp)
par(mfrow = c(1, 1))

### save

save(out.spatial.n400.monitoring,
     out.temporal.n400.monitoring,
     out.spatial.n900.monitoring,
     out.temporal.n900.monitoring,
     out.spatial.n400.satellite,
     out.temporal.n400.satellite,
     out.spatial.n900.satellite,
     out.temporal.n900.satellite, file = "pilotstudy_bayesian_spacetime.RData")

