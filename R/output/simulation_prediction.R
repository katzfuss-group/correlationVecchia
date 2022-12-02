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

### setting ########################################################################

nsim            <- 200
covparms        <- c(1, 0.1, 1.0, 0.5)
nugget          <- 1e-6

m               <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

### 1. Randomly selected locations ################################################# (~18 min)

n               <- 30^2
t               <- 1
n.pred          <- 10^2

output.sptm.pred.random <- parSim_sptm_prediction(cand.m = m, nsim = nsim, n = n/t, n.pred = n.pred, d = 2, t = t, nuggets = nugget, covmodel = cov_matern_spacetime, covparms = covparms, abs.corr = FALSE, method.locs = "random", method.locs.pred = "subset", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
# output.sptm.pred.random <- simulate_spacetime_prediction(nsim = nsim, n = n/t, n.pred = n.pred, d = 2, t = t, m = m, nuggets = nugget, abs.corr = FALSE, method.locs = "random", method.locs.pred = "subset", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel = cov_matern_spacetime, covparms = covparms)

output.sptm.pred.random$time.tot

save(nsim, covparms, nugget, n, t, n.pred, m, output.sptm.pred.random, file = "simout_prediction_allrandom_05242021.RData")

### 2. Monitoring-station locations ################################################ (~31 min)

n               <- 30^2
t               <- 9
n.pred          <- 5 * t

output.sptm.pred.monitoring <- parSim_sptm_prediction(cand.m = m, nsim = nsim, n = n/t, n.pred = n.pred, d = 2, t = t, nuggets = nugget, covmodel = cov_matern_spacetime, covparms = covparms, abs.corr = FALSE, method.locs = "monitoring", method.locs.pred = "subset", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
# output.sptm.pred.monitoring <- simulate_spacetime_prediction(nsim = nsim, n = n/t, n.pred = n.pred, d = 2, t = t, m = m, nuggets = nugget, abs.corr = FALSE, method.locs = "monitoring", method.locs.pred = "subset", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel = cov_matern_spacetime, covparms = covparms)

output.sptm.pred.monitoring$time.tot

save(nsim, covparms, nugget, n, t, n.pred, m, output.sptm.pred.monitoring, file = "simout_prediction_monitoring_05242021.RData")

### 3. Satellite-footprint locations ############################################### (~22 min)

n               <- 30^2
t               <- 1
n.pred          <- 10^2

output.sptm.pred.satellite <- parSim_sptm_prediction(cand.m = m, nsim = nsim, n = n/t, n.pred = n.pred, d = 2, t = t, nuggets = nugget, covmodel = cov_matern_spacetime, covparms = covparms, abs.corr = FALSE, method.locs = "satellite", method.locs.pred = "subset", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, ncores = NULL)
# output.sptm.pred.satellite <- simulate_spacetime_prediction(nsim = nsim, n = n/t, n.pred = n.pred, d = 2, t = t, m = m, nuggets = nugget, abs.corr = FALSE, method.locs = "satellite", method.locs.pred = "subset", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = TRUE, covmodel = cov_matern_spacetime, covparms = covparms)

output.sptm.pred.satellite$time.tot

save(nsim, covparms, nugget, n, t, n.pred, m, output.sptm.pred.satellite, file = "simout_prediction_satellite_05242021.RData")

### visualization ##################################################################

# vdat1 <- data.frame(index = seq(length(output.sptm.pred.random$setting$m)), m = output.sptm.pred.random$setting$m, approx_1 = output.sptm.pred.random$output$logscore[[1]], approx_2 = output.sptm.pred.random$output$logscore[[2]], approx_3 = output.sptm.pred.random$output$logscore[[3]], approx_4 = output.sptm.pred.random$output$logscore[[4]])
# vdat2 <- data.frame(index = seq(length(output.sptm.pred.monitoring$setting$m)), m = output.sptm.pred.monitoring$setting$m, approx_1 = output.sptm.pred.monitoring$output$logscore[[1]], approx_2 = output.sptm.pred.monitoring$output$logscore[[2]], approx_3 = output.sptm.pred.monitoring$output$logscore[[3]], approx_4 = output.sptm.pred.monitoring$output$logscore[[4]])
# vdat3 <- data.frame(index = seq(length(output.sptm.pred.satellite$setting$m)), m = output.sptm.pred.satellite$setting$m, approx_1 = output.sptm.pred.satellite$output$logscore[[1]], approx_2 = output.sptm.pred.satellite$output$logscore[[2]], approx_3 = output.sptm.pred.satellite$output$logscore[[3]], approx_4 = output.sptm.pred.satellite$output$logscore[[4]])
#
# vdat1   <- vdat1 %>% filter(m != 3)
# vdat2   <- vdat2 %>% filter(m != 3)
# vdat3   <- vdat3 %>% filter(m != 3)
#
# vis <- vis_arrange_tria_simple(vdat1 = vdat1, vdat2 = vdat2, vdat3 = vdat3, ftn = function(x) log10(-x), xlim = list(xl1 = NULL, xl2 = NULL, xl3 = NULL), ylim = list(yl1 = NULL, yl2 = NULL, yl3 = NULL), xlab = c("m", "m", "m"), ylab = c("log10(logscore)", "log10(logscore)", "log10(logscore)"), legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), shape = c(18, 15, 17, 16))
#
# ggplot2::ggsave("visout_prediction_05242021.pdf", vis, width = 15.2, height = 5.7)
