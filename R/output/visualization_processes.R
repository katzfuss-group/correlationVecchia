####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of Vecchia-like approximation in various settings.
###
####################################################################################

# use_build_ignore("R/output/visualization_processes.R", escape = TRUE)

set.seed(03222020)

library(correlationVecchia)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

### multivariate case #############################################################################################

cand.d  <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

p1        <- list()
mypal     <- brewer.pal(10, "RdYlBu")
mypal     <- mypal[seq(from = length(mypal), to = 1, by = -1)]
for(i in 1:length(cand.d)) {

  out       <- generate_gp_space(nsim = 1, n = 5^2, d = 2, p = 2, method.locs = "random", method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps, verbose = FALSE, covmodel = cov_bivariate_expo_latDim, covparms = c(1, 0.1, cand.d[i])) 
  covmat    <- out$sim$sim1$covmat

  # covmat %>% fields::image.plot()
  p1[[i]]   <- covmat %>% melt() %>% ggplot(aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) + theme_void() + theme(legend.position="none", plot.title = element_text(hjust=0.5)) + ggtitle(label = paste0("Latent coordinate distance = ", cand.d[i])) + scale_fill_gradientn(colours = mypal, limits = c(0,1))
}

vis.all   <- arrangeGrob(grobs = p1, nrow = 2, ncol = 3)
leg       <- get_legend(covmat %>% melt() %>% ggplot(aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) +  scale_fill_gradientn(colours = mypal, limits = c(0,1)))

p1        <- grid.arrange(vis.all, leg, nrow = 1, ncol = 2, widths = c(6, 1), heights = c(1))

ggplot2::ggsave("p1_multi.pdf", p1, width = 15.2, height = 5.7)


### spacetime case #############################################################################################

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

ggplot2::ggsave("p1_spacetime.pdf", p1, width = 15.2, height = 5.7)

### derivative case #############################################################################################

cand.r  <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)

locs <- matrix(runif(20), 10, 2)

p1        <- list()
mypal     <- brewer.pal(10, "RdYlBu")
mypal     <- mypal[seq(from = length(mypal), to = 1, by = -1)]
for(i in 1:length(cand.r)) {
  cormat <- corr_derivative_matern_2.5_2d(locs = locs, covparms = c(1, cand.r[i])) ; diag(cormat) <- 1
  p1[[i]]   <- cormat %>% melt() %>% ggplot(aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) + theme_void() + theme(legend.position="none", plot.title = element_text(hjust=0.5)) + ggtitle(label = paste0("Distance matrix with range parameter = ", cand.r[i])) + scale_fill_gradientn(colours = mypal, limits = c(-1,1))
}

vis.all   <- arrangeGrob(grobs = p1, nrow = 2, ncol = 4)
leg       <- get_legend(cormat %>% melt() %>% ggplot(aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) +  scale_fill_gradientn(colours = mypal, limits = c(-1,1)))

p1        <- grid.arrange(vis.all, leg, nrow = 1, ncol = 2, widths = c(9, 1), heights = c(1))

ggplot2::ggsave("p1_deriv.pdf", p1, width = 15.2, height = 5.7)


p2        <- list()
mypal     <- brewer.pal(10, "RdYlBu")
mypal     <- mypal[seq(from = length(mypal), to = 1, by = -1)]
for(i in 1:length(cand.r)) {
  cormat <- corr_derivative_matern_2.5_2d(locs = locs, covparms = c(1, cand.r[i])) ; diag(cormat) <- 1
  cormat <- abs(cormat)
  p2[[i]]   <- cormat %>% melt() %>% ggplot(aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) + theme_void() + theme(legend.position="none", plot.title = element_text(hjust=0.5)) + ggtitle(label = paste0("Distance matrix with range parameter = ", cand.r[i])) + scale_fill_gradientn(colours = mypal, limits = c(-1,1))
}

vis.all   <- arrangeGrob(grobs = p2, nrow = 2, ncol = 4)
leg       <- get_legend(cormat %>% melt() %>% ggplot(aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) +  scale_fill_gradientn(colours = mypal, limits = c(-1,1)))

p2        <- grid.arrange(vis.all, leg, nrow = 1, ncol = 2, widths = c(9, 1), heights = c(1))

ggplot2::ggsave("p2_deriv.pdf", p2, width = 15.2, height = 5.7)


par(mfrow = c(1, 2))

locs <- seq(from = 0, to = 0.6, by = 0.01)
locs <- matrix(locs, nrow = length(locs), ncol = 1)
locs <- cbind(locs, 0)

covparms = c(1, 0.1)

covvec.gp <- correlationVecchia::cov_matern_2.5(locs = locs, covparms = covparms)[1, ]
covvec.dgp <- correlationVecchia:::.partial_cov_matern_2.5(locs = locs, covparms = covparms, coord = 1)[1, ]
covvec.d2gp <- correlationVecchia:::.double_partial_cov_matern_2.5(locs = locs, covparms = covparms, coord = 1)[1, ]

plot(locs[, 1], covvec.gp, ylim = c(-40, 170), type = 'l', col = 'black', xlab = 'distance', ylab = 'cov', lwd = 2)
lines(locs[, 1], covvec.dgp, col = 'red', lwd = 2)
lines(locs[, 1], covvec.d2gp, col = 'blue', lwd = 2)
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = c("cov(GP, GP)", "cov(GP', GP)", "cov(GP', GP')"), lty = 1, lwd = 2, col = c("black", "red", "blue"))

corvec.gp <- covvec.gp / covvec.gp[1]
corvec.dgp <- covvec.dgp / sqrt(covvec.gp[1] * covvec.d2gp[1])
corvec.d2gp <- covvec.d2gp / covvec.d2gp[1]

plot(locs[, 1], corvec.gp, ylim = c(-0.5, 1), type = 'l', col = 'black', xlab = 'distance', ylab = 'corr', lwd = 2)
lines(locs[, 1], corvec.dgp, col = 'red', lwd = 2)
lines(locs[, 1], corvec.d2gp, col = 'blue', lwd = 2)
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = c("corr(GP, GP)", "corr(GP', GP)", "corr(GP', GP')"), lty = 1, lwd = 2, col = c("black", "red", "blue"))

par(mfrow = c(1, 1))


par(mfrow = c(1, 2))

locs <- seq(from = 0, to = 3, by = 0.1)
locs <- matrix(locs, nrow = length(locs), ncol = 1)
locs <- cbind(locs, 0)

covparms = c(1, 1)

covvec.gp <- correlationVecchia::cov_matern_2.5(locs = locs, covparms = covparms)[1, ]
covvec.dgp <- correlationVecchia:::.partial_cov_matern_2.5(locs = locs, covparms = covparms, coord = 1)[1, ]
covvec.d2gp <- correlationVecchia:::.double_partial_cov_matern_2.5(locs = locs, covparms = covparms, coord = 1)[1, ]

plot(locs[, 1], covvec.gp, ylim = c(-0.5, 1.7), type = 'l', col = 'black', xlab = 'distance', ylab = 'cov', lwd = 2)
lines(locs[, 1], covvec.dgp, col = 'red', lwd = 2)
lines(locs[, 1], covvec.d2gp, col = 'blue', lwd = 2)
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = c("cov(GP, GP)", "cov(GP', GP)", "cov(GP', GP')"), lty = 1, lwd = 2, col = c("black", "red", "blue"))

corvec.gp <- covvec.gp / covvec.gp[1]
corvec.dgp <- covvec.dgp / sqrt(covvec.gp[1] * covvec.d2gp[1])
corvec.d2gp <- covvec.d2gp / covvec.d2gp[1]

plot(locs[, 1], corvec.gp, ylim = c(-0.5, 1), type = 'l', col = 'black', xlab = 'distance', ylab = 'corr', lwd = 2)
lines(locs[, 1], corvec.dgp, col = 'red', lwd = 2)
lines(locs[, 1], corvec.d2gp, col = 'blue', lwd = 2)
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = c("corr(GP, GP)", "corr(GP', GP)", "corr(GP', GP')"), lty = 1, lwd = 2, col = c("black", "red", "blue"))

par(mfrow = c(1, 1))
