####################################################################
####
####  Author: Myeongjong (MJ) Kang (kmj.stat@gmail.com)
####
####  Description:
####
####################################################################

gc()
rm(list = ls())

library(GPvecchia)

library(foreach)

## To visualize results
library(tidyr) ; library(ggplot2) ; library(gridExtra) ; library(RColorBrewer)

source("2_corrvecchia/vecchia_specify_adjusted.R")
source("2_corrvecchia/corrvecchia_with_ab_corr_dist.R")
source("2_corrvecchia/kldiv.R")

set.seed(10102019)

####################################################################
#### Nonstaionary Matern covariance model (Risser MD, Calder CA (2015))
####################################################################

# covparms = c(R, period), 2 * pi / omega is a period
covparms <- c(1, 1/10)

# cov.aniso <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])

# 0 <= covparms[1] <= 3period / 2pi 
# # 3 / (2 * pi) = 0.4774648
covmodel_periodic <- function(locs, covparms) {
  h <- fields::rdist(locs)
  exp(-3 * h / covparms[1]) * cos(2 * pi * h / covparms[2])
}

covparms <- c(1 * 0.1, 1)
covparms[1] <= 3 * covparms[2] / (2 * pi) 

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')

covparms <- c(1 * 0.2, 1)
covparms[1] <= 3 * covparms[2] / (2 * pi) 

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')

covparms <- c(1 * 0.3, 1)
covparms[1] <= 3 * covparms[2] / (2 * pi) 

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')

covparms <- c(1 * 0.4, 1)
covparms[1] <= 3 * covparms[2] / (2 * pi) 

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')

covparms <- c(1 * 0.45, 1)
covparms[1] <= 3 * covparms[2] / (2 * pi) 

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')

covparms <- c(0.2 * 0.45, 0.2)
covparms[1] <= 3 * covparms[2] / (2 * pi) 

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')

covparms <- c(0.2 * 0.2, 0.2)
covparms[1] <= 3 * covparms[2] / (2 * pi) 

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')


####################################################################
#### simulation function
####################################################################