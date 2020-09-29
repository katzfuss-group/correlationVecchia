####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: .
###
####################################################################################

rm(list = ls())

set.seed(09292020)

library(correlationVecchia)
library(GpGp)
library(dplyr)
library(foreach)
library(ggplot2)
library(gridExtra)

### spacetime case - visualization #############################################################################################

out1  <- generate_gp_spacetime(nsim = 1, n = 100, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out2  <- generate_gp_spacetime(nsim = 1, n = 5^2, d = 2, t.len = 5, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out3  <- generate_gp_spacetime(nsim = 1, n = 5^2, d = 2, t.len = 5, method.locs = "all.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out4  <- generate_gp_spacetime(nsim = 1, n = 100, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))

locs <- list()
locs[[1]] <- out1$sim$sim1$locs %>% cbind(seq(100)) %>% as.data.frame()
locs[[2]] <- out2$sim$sim1$locs %>% cbind(seq(5^3)) %>% as.data.frame()
locs[[3]] <- out3$sim$sim1$locs %>% cbind(seq(5^3)) %>% as.data.frame()
locs[[4]] <- out4$sim$sim1$locs %>% cbind(seq(100)) %>% as.data.frame()

for(i in 1:4) colnames(locs[[i]]) <- c("x", "y", "t", "index")

k <- 1
p1 <- list()
for(i in 1:4) {
  
  p1[[i]] <- locs[[i]] %>% ggplot() + geom_point(aes(x = x, y = y), size = 0.6) + xlab("x") + ylab("y")
  k <- k + 1
  
  p1[[i+4]] <- locs[[i]] %>% ggplot() + geom_point(aes(x = index, y = t), size = 0.6) + xlab("index") + ylab("time")
  k <- k + 1
}

vis.all   <- arrangeGrob(grobs = p1, nrow = 2, ncol = 4)
p1        <- grid.arrange(vis.all)

### spacetime case - generation #############################################################################################

out1  <- generate_gp_spacetime(nsim = 1, n = 20^2, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out2  <- generate_gp_spacetime(nsim = 1, n = 20, d = 2, t.len = 4, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out3  <- generate_gp_spacetime(nsim = 1, n = 20, d = 2, t.len = 4, method.locs = "all.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out4  <- generate_gp_spacetime(nsim = 1, n = 20^2, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))

locs <- list()
locs[[1]] <- out1$sim$sim1$locs %>% as.matrix()
locs[[2]] <- out2$sim$sim1$locs %>% as.matrix()
locs[[3]] <- out3$sim$sim1$locs %>% as.matrix()
locs[[4]] <- out4$sim$sim1$locs %>% as.matrix()

### spacetime case - process #############################################################################################

covparms <- c(1, 0.1, 0.01, 3.5, 0) # var, range_space, range_time, nu, nugget

covmat <- list()
covmat[[1]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[1]])
covmat[[2]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[2]])
covmat[[3]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[3]])
covmat[[4]] <- GpGp::matern_spacetime(covparms = covparms, locs = locs[[4]])

y <- list()
y[[1]] <- t(chol(covmat[[1]])) %*% rnorm(n = nrow(covmat[[1]]))
y[[2]] <- t(chol(covmat[[2]])) %*% rnorm(n = nrow(covmat[[2]]))
y[[3]] <- t(chol(covmat[[3]])) %*% rnorm(n = nrow(covmat[[3]]))
y[[4]] <- t(factorize(covmat[[1]], method = "eigen-I", tol = 1e-6)$covfactor) %*% rnorm(n = nrow(covmat[[4]])) # y[[4]] <- t(chol(covmat[[4]])) %*% rnorm(n = nrow(covmat[[4]]))

### spacetime case - fit #############################################################################################

fit <- list()
fit[[1]] <- fit_corrvecchia(y = y[[1]], inputs = locs[[1]], trend = "intercept", scale = "parms", covfun = "matern_spacetime", max.it = 100)
fit[[2]] <- fit_corrvecchia(y = y[[2]], inputs = locs[[2]], trend = "intercept", scale = "parms", covfun = "matern_spacetime", max.it = 100)
fit[[3]] <- fit_corrvecchia(y = y[[3]], inputs = locs[[3]], trend = "intercept", scale = "parms", covfun = "matern_spacetime", max.it = 100)
fit[[4]] <- fit_corrvecchia(y = y[[4]], inputs = locs[[4]], trend = "intercept", scale = "parms", covfun = "matern_spacetime", max.it = 100)

fit[[1]]$covparms ; fit[[1]]$conv
fit[[2]]$covparms ; fit[[2]]$conv
fit[[3]]$covparms ; fit[[3]]$conv
fit[[4]]$covparms ; fit[[4]]$conv




