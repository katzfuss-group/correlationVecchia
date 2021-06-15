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

library(vecchia)
library(foreach)

##### ncores #######################################################################

no_cores  <- parallel::detectCores() - 2

### setting ########################################################################

ms              <- c(5, 10, 20, 30, 40)
covparms        <- c(1, 0.1, 1.0, 0.5)
nugget          <- 0.00

n               <- 30^2
d               <- 2

### n = 900, satellite footprint #######################################

t               <- 1

## 

Sys.time()

output.sptm.posterior.satellite.srange <- parSim_sptm_posterior(cand.m = ms, target = "srange", ordfix = TRUE, n = n, d = d, t = t, nuggets = nugget, method.locs = "satellite", N = 200, xlim = c(0.09, 0.13), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = -4, covparms = covparms, ncores = NULL)

Sys.time()

output.sptm.posterior.satellite.trange <- parSim_sptm_posterior(cand.m = ms, target = "trange", ordfix = TRUE, n = n, d = d, t = t, nuggets = nugget, method.locs = "satellite", N = 200, xlim = c(0.05, 3.05), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = -4, covparms = covparms, ncores = NULL)

Sys.time()

##

par(mfrow = c(3, 5))

for(i in 1:length(output.sptm.posterior.satellite.srange$simout)) matplot(output.sptm.posterior.satellite.srange$simout[[i]]$alpha, output.sptm.posterior.satellite.srange$simout[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.satellite.trange$simout)) matplot(output.sptm.posterior.satellite.trange$simout[[i]]$alpha, output.sptm.posterior.satellite.trange$simout[[i]]$post.norm, type = "l")

par(mfrow = c(1, 1))

##

save(ms, covparms, nugget, n, t, d, output.sptm.posterior.satellite.srange, output.sptm.posterior.satellite.trange, file = "simout_posterior_satellite_05242021.RData")

### n = 900, monitoring station #######################################

t               <- 36

##

Sys.time()

output.sptm.posterior.monitoring.srange <- parSim_sptm_posterior(cand.m = ms, target = "srange", ordfix = TRUE, n = n/t, d = d, t = t, nuggets = nugget, method.locs = "monitoring", N = 200, xlim = c(0.005, 0.115), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = 4, covparms = covparms, ncores = NULL)

Sys.time()

output.sptm.posterior.monitoring.trange <- parSim_sptm_posterior(cand.m = ms, target = "trange", ordfix = TRUE, n = n/t, d = d, t = t, nuggets = nugget, method.locs = "monitoring", N = 200, xlim = c(0.75, 1.25), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = 4, covparms = covparms, ncores = NULL)

Sys.time()

##

par(mfrow = c(3, 5))

for(i in 1:length(output.sptm.posterior.monitoring.srange$simout)) matplot(output.sptm.posterior.monitoring.srange$simout[[i]]$alpha, output.sptm.posterior.monitoring.srange$simout[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.monitoring.trange$simout)) matplot(output.sptm.posterior.monitoring.trange$simout[[i]]$alpha, output.sptm.posterior.monitoring.trange$simout[[i]]$post.norm, type = "l")

par(mfrow = c(1, 1))

##

save(ms, covparms, nugget, n, t, d, output.sptm.posterior.monitoring.srange, output.sptm.posterior.monitoring.trange, file = "simout_posterior_monitoring_05242021.RData")

### n = 900, random selection ##########################################

t               <- 1

##

Sys.time()

output.sptm.posterior.random.srange <- parSim_sptm_posterior(cand.m = ms, target = "srange", ordfix = TRUE, n = n, d = d, t = t, nuggets = nugget, method.locs = "random", N = 200, xlim = c(0.06, 0.14), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = 4, covparms = covparms, ncores = NULL)

Sys.time()

output.sptm.posterior.random.trange <- parSim_sptm_posterior(cand.m = ms, target = "trange", ordfix = TRUE, n = n, d = d, t = t, nuggets = nugget, method.locs = "random", N = 200, xlim = c(0.05, 3.05), sdlog = 0.6, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, tol.dec = 4, covparms = covparms, ncores = NULL)

Sys.time()

##

par(mfrow = c(3, 5))

for(i in 1:length(output.sptm.posterior.random.srange$simout)) matplot(output.sptm.posterior.random.srange$simout[[i]]$alpha, output.sptm.posterior.random.srange$simout[[i]]$post.norm, type = "l")

for(i in 1:length(output.sptm.posterior.random.trange$simout)) matplot(output.sptm.posterior.random.trange$simout[[i]]$alpha, output.sptm.posterior.random.trange$simout[[i]]$post.norm, type = "l")

par(mfrow = c(1, 1))

##

save(ms, covparms, nugget, n, t, d, output.sptm.posterior.random.srange, output.sptm.posterior.random.trange, file = "simout_posterior_allrandom_05242021.RData")
