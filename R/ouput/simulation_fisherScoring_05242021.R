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

library(vecchia) ; library(dplyr) ; library(tidyr) ; library(ggplot2) ; library(gridExtra) 

### setting ########################################################################

nsim      <- 200
cand.m    <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)
covparms  <- c(1, 0.1, 1.0, 0.5, 0)

# nsim      <- 2
# cand.m    <- c(10, 30)
# covparms  <- c(1, 0.1, 0.1, 0.5, 0)

### monitoring station ############################################################# (~ 20 hours)

Sys.time()

output.monitoring <- parSim_sptm_Fisher(cand.m = cand.m, nsim = nsim, n = 25, d = 2, t = 36, covmodel = GpGp::matern_spacetime, srange.ini = covparms[2], trange.ini = covparms[2], covparms = covparms, method.locs = "monitoring", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL, tol.dec = 4)

Sys.time()

## save

save(nsim, cand.m, covparms, output.monitoring, file = "simout_fisher_monitoring_05242021.RData")

### satellite footprint ############################################################ (~ 20 hours)

Sys.time()

output.satellite <- parSim_sptm_Fisher(cand.m = cand.m, nsim = nsim, n = 30^2, d = 2, t = 1, covmodel = GpGp::matern_spacetime, srange.ini = covparms[2], trange.ini = covparms[3], covparms = covparms, method.locs = "satellite", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL, tol.dec = -4)

Sys.time()

## save

save(nsim, cand.m, covparms, output.satellite, file = "simout_fisher_satellite_05242021.RData")

### visualization ##################################################################

vdat1   <- output.monitoring$kls.average ; vdat1[, !(names(vdat1) %in% c("index", "m"))] <- log10(vdat1[, !(names(vdat1) %in% c("index", "m"))])
vdat2   <- data.frame(index = seq(length(output.monitoring$kls.average$m)), m = output.monitoring$kls.average$m, approx_1 = sqrt(output.monitoring$msd.srange[[1]]), approx_2 = sqrt(output.monitoring$msd.srange[[2]]), approx_3 = sqrt(output.monitoring$msd.srange[[3]]), approx_4 = sqrt(output.monitoring$msd.srange[[4]]))
vdat3   <- data.frame(index = seq(length(output.monitoring$kls.average$m)), m = output.monitoring$kls.average$m, approx_1 = sqrt(output.monitoring$msd.trange[[1]]), approx_2 = sqrt(output.monitoring$msd.trange[[2]]), approx_3 = sqrt(output.monitoring$msd.trange[[3]]), approx_4 = sqrt(output.monitoring$msd.trange[[4]]))

vdat1   <- vdat1 %>% filter(m != 3)
vdat2   <- vdat2 %>% filter(m != 3)
vdat3   <- vdat3 %>% filter(m != 3)
  
ylim1   <- vdat1 %>% select(approx_1, approx_2, approx_3, approx_4, approx_5) %>% range()
ylim2   <- vdat2 %>% select(approx_1, approx_2, approx_3, approx_4) %>% range()
ylim3   <- vdat3 %>% select(approx_2, approx_3, approx_4) %>% range()

vis.mon <- vis_arrange_tria(vdat1 = vdat1, vdat2 = vdat2, vdat3 = vdat3, ftn = function(x) x, xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = ylim1, yl2 = ylim2, yl3 = ylim3), xlab = c("m", "m", "m"), ylab = c("log10(KL)", "RMSD for spatial range", "RMSD for temporal range"), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), shape = c(18, 15, 17, 16, NA))

vdat1   <- output.satellite$kls.average ; vdat1[, !(names(vdat1) %in% c("index", "m"))] <- log10(vdat1[, !(names(vdat1) %in% c("index", "m"))])
vdat2   <- data.frame(index = seq(length(output.satellite$kls.average$m)), m = output.satellite$kls.average$m, approx_1 = sqrt(output.satellite$msd.srange[[1]]), approx_2 = sqrt(output.satellite$msd.srange[[2]]), approx_3 = sqrt(output.satellite$msd.srange[[3]]), approx_4 = sqrt(output.satellite$msd.srange[[4]]))
vdat3   <- data.frame(index = seq(length(output.satellite$kls.average$m)), m = output.satellite$kls.average$m, approx_1 = sqrt(output.satellite$msd.trange[[1]]), approx_2 = sqrt(output.satellite$msd.trange[[2]]), approx_3 = sqrt(output.satellite$msd.trange[[3]]), approx_4 = sqrt(output.satellite$msd.trange[[4]]))

# vdat2   <- data.frame(index = seq(length(output.satellite$kls.average$m)), m = output.satellite$kls.average$m, approx_1 = sqrt(output.satellite$mse.srange[[1]]), approx_2 = sqrt(output.satellite$mse.srange[[2]]), approx_3 = sqrt(output.satellite$mse.srange[[3]]), approx_4 = sqrt(output.satellite$mse.srange[[4]]))
# vdat3   <- data.frame(index = seq(length(output.satellite$kls.average$m)), m = output.satellite$kls.average$m, approx_1 = sqrt(output.satellite$mse.trange[[1]]), approx_2 = sqrt(output.satellite$mse.trange[[2]]), approx_3 = sqrt(output.satellite$mse.trange[[3]]), approx_4 = sqrt(output.satellite$mse.trange[[4]]))

vdat1   <- vdat1 %>% filter(m != 3)
vdat2   <- vdat2 %>% filter(m != 3)
vdat3   <- vdat3 %>% filter(m != 3)

ylim1   <- vdat1 %>% select(approx_1, approx_2, approx_3, approx_4, approx_5) %>% range()
ylim2   <- vdat2 %>% select(approx_1, approx_2, approx_3, approx_4) %>% range()
ylim3   <- vdat3 %>% select(approx_2, approx_3, approx_4) %>% range()

vis.sat <- vis_arrange_tria(vdat1 = vdat1, vdat2 = vdat2, vdat3 = vdat3, ftn = function(x) x, xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = ylim1, yl2 = ylim2, yl3 = ylim3), xlab = c("m", "m", "m"), ylab = c("log10(KL)", "RMSD for spatial range", "RMSD for temporal range"), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "Exact GP"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray30"), shape = c(18, 15, 17, 16, NA))


