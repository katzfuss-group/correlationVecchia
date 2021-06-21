####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of approximations in various settings.
###
####################################################################################

rm(list = ls())
set.seed(05072021)

##

gc()

memory.size() # memory.limit()

##

library(correlationVecchia)
library(foreach)

##### ncores #######################################################################

no_cores  <- parallel::detectCores() - 2

### setting ########################################################################

ms              <- c(5, 10, 20, 30, 40)
covparms        <- c(1, 0.1, 1.0, 0.5)
nugget          <- 0.4

n               <- 30^2
d               <- 2

### n = 900, satellite footprint #######################################

t               <- 1

## (~ 4.5 hrs)

Sys.time()

output.sptm.posterior.satellite.srange <- parSim_sptm_posterior(cand.m = ms, target = "srange", ordfix = TRUE, n = n, d = d, t = t, nuggets = nugget, method.locs = "satellite", N = 200, xlim = c(0.005, 0.155), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = -4, covparms = covparms, ncores = NULL)

Sys.time()

output.sptm.posterior.satellite.trange <- parSim_sptm_posterior(cand.m = ms, target = "trange", ordfix = TRUE, n = n, d = d, t = t, nuggets = nugget, method.locs = "satellite", N = 200, xlim = c(0.05, 3.05), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = -4, covparms = covparms, ncores = NULL)

Sys.time()

##

par(mfrow = c(1, 5))

for(i in 1:length(output.sptm.posterior.satellite.srange$simout)) matplot(output.sptm.posterior.satellite.srange$simout[[i]]$alpha, output.sptm.posterior.satellite.srange$simout[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.satellite.srange$simout.ic0)) matplot(output.sptm.posterior.satellite.srange$simout.ic0[[i]]$alpha, output.sptm.posterior.satellite.srange$simout.ic0[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.satellite.trange$simout)) matplot(output.sptm.posterior.satellite.trange$simout[[i]]$alpha, output.sptm.posterior.satellite.trange$simout[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.satellite.trange$simout.ic0)) matplot(output.sptm.posterior.satellite.trange$simout.ic0[[i]]$alpha, output.sptm.posterior.satellite.trange$simout.ic0[[i]]$post.norm, type = "l")

par(mfrow = c(1, 1))

##

save(ms, covparms, nugget, n, t, d, output.sptm.posterior.satellite.srange, output.sptm.posterior.satellite.trange, file = "simout_posterior_satellite_05242021.RData")

### n = 900, monitoring station #######################################

t               <- 36

## (~ 4.5 hrs)

Sys.time()

output.sptm.posterior.monitoring.srange <- parSim_sptm_posterior(cand.m = ms, target = "srange", ordfix = TRUE, n = n/t, d = d, t = t, nuggets = nugget, method.locs = "monitoring", N = 200, xlim = c(0.005, 0.115), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = 4, covparms = covparms, ncores = NULL)

Sys.time()

output.sptm.posterior.monitoring.trange <- parSim_sptm_posterior(cand.m = ms, target = "trange", ordfix = TRUE, n = n/t, d = d, t = t, nuggets = nugget, method.locs = "monitoring", N = 200, xlim = c(0.05, 3.05), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = 4, covparms = covparms, ncores = NULL)

Sys.time()

##

par(mfrow = c(1, 5))

for(i in 1:length(output.sptm.posterior.monitoring.srange$simout)) matplot(output.sptm.posterior.monitoring.srange$simout[[i]]$alpha, output.sptm.posterior.monitoring.srange$simout[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.monitoring.srange$simout.ic0)) matplot(output.sptm.posterior.monitoring.srange$simout.ic0[[i]]$alpha, output.sptm.posterior.monitoring.srange$simout.ic0[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.monitoring.trange$simout)) matplot(output.sptm.posterior.monitoring.trange$simout[[i]]$alpha, output.sptm.posterior.monitoring.trange$simout[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.monitoring.trange$simout.ic0)) matplot(output.sptm.posterior.monitoring.trange$simout.ic0[[i]]$alpha, output.sptm.posterior.monitoring.trange$simout.ic0[[i]]$post.norm, type = "l")

par(mfrow = c(1, 1))

##

save(ms, covparms, nugget, n, t, d, output.sptm.posterior.monitoring.srange, output.sptm.posterior.monitoring.trange, file = "simout_posterior_monitoring_05242021.RData")

### n = 900, random selection ##########################################

t               <- 1

## (~ 4.5 hrs)

Sys.time()

output.sptm.posterior.random.srange <- parSim_sptm_posterior(cand.m = ms, target = "srange", ordfix = TRUE, n = n, d = d, t = t, nuggets = nugget, method.locs = "random", N = 200, xlim = c(0.04, 0.16), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = 4, covparms = covparms, ncores = NULL)

Sys.time()

output.sptm.posterior.random.trange <- parSim_sptm_posterior(cand.m = ms, target = "trange", ordfix = TRUE, n = n, d = d, t = t, nuggets = nugget, method.locs = "random", N = 200, xlim = c(0.05, 3.05), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = 4, covparms = covparms, ncores = NULL)

Sys.time()

##

par(mfrow = c(1, 5))

for(i in 1:length(output.sptm.posterior.random.srange$simout)) matplot(output.sptm.posterior.random.srange$simout[[i]]$alpha, output.sptm.posterior.random.srange$simout[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.random.srange$simout.ic0)) matplot(output.sptm.posterior.random.srange$simout.ic0[[i]]$alpha, output.sptm.posterior.random.srange$simout.ic0[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.random.trange$simout)) matplot(output.sptm.posterior.random.trange$simout[[i]]$alpha, output.sptm.posterior.random.trange$simout[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.random.trange$simout.ic0)) matplot(output.sptm.posterior.random.trange$simout.ic0[[i]]$alpha, output.sptm.posterior.random.trange$simout.ic0[[i]]$post.norm, type = "l")

par(mfrow = c(1, 1))

##

save(ms, covparms, nugget, n, t, d, output.sptm.posterior.random.srange, output.sptm.posterior.random.trange, file = "simout_posterior_allrandom_05242021.RData")

### visualization ##################################################################

cols      <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70")
legends   <- c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP")

### random selection

output    <- output.sptm.posterior.random.srange

candids   <- output$setting$candid
ms        <- unique(candids$m)
approxs   <- unique(candids$approx)

idx.list <- list()
for(i in 1:length(ms)) idx.list[[i]] <- which(candids$m == ms[i])

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout[[idx.list[[i]][1]]]$alpha, output$simout[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout.ic0[[idx.list[[i]][1]]]$alpha, output$simout.ic0[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

#

output    <- output.sptm.posterior.random.trange

candids   <- output$setting$candid
ms        <- unique(candids$m)
approxs   <- unique(candids$approx)

idx.list <- list()
for(i in 1:length(ms)) idx.list[[i]] <- which(candids$m == ms[i])

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout[[idx.list[[i]][1]]]$alpha, output$simout[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout.ic0[[idx.list[[i]][1]]]$alpha, output$simout.ic0[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

## monitoring station

output    <- output.sptm.posterior.monitoring.srange

candids   <- output$setting$candid
ms        <- unique(candids$m)
approxs   <- unique(candids$approx)

idx.list <- list()
for(i in 1:length(ms)) idx.list[[i]] <- which(candids$m == ms[i])

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout[[idx.list[[i]][1]]]$alpha, output$simout[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout.ic0[[idx.list[[i]][1]]]$alpha, output$simout.ic0[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

#

output    <- output.sptm.posterior.monitoring.trange

candids   <- output$setting$candid
ms        <- unique(candids$m)
approxs   <- unique(candids$approx)

idx.list <- list()
for(i in 1:length(ms)) idx.list[[i]] <- which(candids$m == ms[i])

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout[[idx.list[[i]][1]]]$alpha, output$simout[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout.ic0[[idx.list[[i]][1]]]$alpha, output$simout.ic0[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

## satellite footprint

output    <- output.sptm.posterior.satellite.srange

candids   <- output$setting$candid
ms        <- unique(candids$m)
approxs   <- unique(candids$approx)

idx.list <- list()
for(i in 1:length(ms)) idx.list[[i]] <- which(candids$m == ms[i])

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout[[idx.list[[i]][1]]]$alpha, output$simout[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout.ic0[[idx.list[[i]][1]]]$alpha, output$simout.ic0[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

#

output    <- output.sptm.posterior.satellite.trange

candids   <- output$setting$candid
ms        <- unique(candids$m)
approxs   <- unique(candids$approx)

idx.list <- list()
for(i in 1:length(ms)) idx.list[[i]] <- which(candids$m == ms[i])

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout[[idx.list[[i]][1]]]$alpha, output$simout[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout[[j]]$alpha, output$simout[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(output$simout.ic0[[idx.list[[i]][1]]]$alpha, output$simout.ic0[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", ms[i]))

  if(i == length(idx.list)) legend("topright", legend = legends, col = cols, lwd = 3)

  j <- idx.list[[i]][1]
  k <- 1
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][2]
  k <- 2
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][3]
  k <- 3
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][4]
  k <- 4
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")

  j <- idx.list[[i]][5]
  k <- 5
  lines(output$simout.ic0[[j]]$alpha, output$simout.ic0[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "dotted")

}
par(mfrow = c(1, 1))

### visualization ##################################################################

library(ggplot2) ; library(dplyr) ; library(gridExtra)

vis.sptm.post.random.srange <- vis_arrange_posterior(output = output.sptm.posterior.random.srange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

vis.sptm.post.random.trange <- vis_arrange_posterior(output = output.sptm.posterior.random.trange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

ggplot2::ggsave("visout_sptm_post_random_s_05242021.pdf", vis.sptm.post.random.srange, width = 15.2, height = 5.7)

ggplot2::ggsave("visout_sptm_post_random_t_05242021.pdf", vis.sptm.post.random.trange, width = 15.2, height = 5.7)

vis.sptm.post.random.srange <- vis_arrange_posterior_ic0(output = output.sptm.posterior.random.srange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

vis.sptm.post.random.trange <- vis_arrange_posterior_ic0(output = output.sptm.posterior.random.trange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

ggplot2::ggsave("visout_sptm_post_ic0_random_s_05242021.pdf", vis.sptm.post.random.srange, width = 15.2, height = 5.7)

ggplot2::ggsave("visout_sptm_post_ic0_random_t_05242021.pdf", vis.sptm.post.random.trange, width = 15.2, height = 5.7)


vis.sptm.post.monitoring.srange <- vis_arrange_posterior(output = output.sptm.posterior.monitoring.srange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

vis.sptm.post.monitoring.trange <- vis_arrange_posterior(output = output.sptm.posterior.monitoring.trange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

ggplot2::ggsave("visout_sptm_post_monitoring_s_05242021.pdf", vis.sptm.post.monitoring.srange, width = 15.2, height = 5.7)

ggplot2::ggsave("visout_sptm_post_monitoring_t_05242021.pdf", vis.sptm.post.monitoring.trange, width = 15.2, height = 5.7)

vis.sptm.post.monitoring.srange <- vis_arrange_posterior_ic0(output = output.sptm.posterior.monitoring.srange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

vis.sptm.post.monitoring.trange <- vis_arrange_posterior_ic0(output = output.sptm.posterior.monitoring.trange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

ggplot2::ggsave("visout_sptm_post_ic0_monitoring_s_05242021.pdf", vis.sptm.post.monitoring.srange, width = 15.2, height = 5.7)

ggplot2::ggsave("visout_sptm_post_ic0_monitoring_t_05242021.pdf", vis.sptm.post.monitoring.trange, width = 15.2, height = 5.7)


vis.sptm.post.satellite.srange <- vis_arrange_posterior(output = output.sptm.posterior.satellite.srange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

vis.sptm.post.satellite.trange <- vis_arrange_posterior(output = output.sptm.posterior.satellite.trange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

ggplot2::ggsave("visout_sptm_post_satellite_s_05242021.pdf", vis.sptm.post.satellite.srange, width = 15.2, height = 5.7)

ggplot2::ggsave("visout_sptm_post_satellite_t_05242021.pdf", vis.sptm.post.satellite.trange, width = 15.2, height = 5.7)

vis.sptm.post.satellite.srange <- vis_arrange_posterior_ic0(output = output.sptm.posterior.satellite.srange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

vis.sptm.post.satellite.trange <- vis_arrange_posterior_ic0(output = output.sptm.posterior.satellite.trange, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "True"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), xlim = NULL, ylim = NULL, xlab = "x", ylab = "y", alpha = 0.8, size.line = 1, size.legend = 14, size.lab = 14, size.text = 12, size.margin = c(5.5, 20, 5.5, 5.5))

ggplot2::ggsave("visout_sptm_post_ic0_satellite_s_05242021.pdf", vis.sptm.post.satellite.srange, width = 15.2, height = 5.7)

ggplot2::ggsave("visout_sptm_post_ic0_satellite_t_05242021.pdf", vis.sptm.post.satellite.trange, width = 15.2, height = 5.7)

