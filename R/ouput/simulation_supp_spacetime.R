####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of approximations in various settings.
###
####################################################################################

rm(list = ls())
set.seed(05242021)

##

gc()

memory.size() # memory.limit()

##

library(correlationVecchia) ; library(dplyr) ; library(tidyr) ; library(ggplot2) ; library(gridExtra)

### spacetime case 1 ################################################################# (~ 20 min)

Sys.time()

nsim    <- 10
cand.m  <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

output.sptm.gen1 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t = 1, covmodel = cov_expo_spacetime_nonsep, covparms = c(1, 0.75, 50, 25), abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

output.sptm.gen2 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 10^2, d = 2, t = 9, covmodel = cov_expo_spacetime_nonsep, covparms = c(1, 0.75, 50, 25), abs.corr = FALSE, method.locs = "monitoring", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

output.sptm.gen3 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 10^2, d = 2, t = 9, covmodel = cov_expo_spacetime_nonsep, covparms = c(1, 0.75, 50, 25), abs.corr = FALSE, method.locs = "grid", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

output.sptm.gen4 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t = 1, covmodel = cov_expo_spacetime_nonsep, covparms = c(1, 0.75, 50, 25), abs.corr = FALSE, method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## save

save(nsim, cand.m, output.sptm.gen1, output.sptm.gen2, output.sptm.gen3, output.sptm.gen4, file = "simout_sptm_supp1_05242021.RData")

rm(output.sptm.gen1, output.sptm.gen2, output.sptm.gen3, output.sptm.gen4)

### spacetime case 2 ################################################################# (~ 20 min)

output.sptm.gen1 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t = 1, covmodel = cov_expo_spacetime_nonsep, covparms = c(1, 0.75, 500, 2.5), abs.corr = FALSE, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

output.sptm.gen2 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 10^2, d = 2, t = 9, covmodel = cov_expo_spacetime_nonsep, covparms = c(1, 0.75, 500, 2.5), abs.corr = FALSE, method.locs = "monitoring", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

output.sptm.gen3 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 10^2, d = 2, t = 9, covmodel = cov_expo_spacetime_nonsep, covparms = c(1, 0.75, 500, 2.5), abs.corr = FALSE, method.locs = "grid", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

output.sptm.gen4 <- parSim_sptm_knownCovparms(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t = 1, covmodel = cov_expo_spacetime_nonsep, covparms = c(1, 0.75, 500, 2.5), abs.corr = FALSE, method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)

## save

save(nsim, cand.m, output.sptm.gen1, output.sptm.gen2, output.sptm.gen3, output.sptm.gen4, file = "simout_sptm_supp2_05242021.RData")

rm(nsim, cand.m, output.sptm.gen1, output.sptm.gen2, output.sptm.gen3, output.sptm.gen4)
