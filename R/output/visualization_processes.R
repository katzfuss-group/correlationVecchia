####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of Vecchia-like approximation in various settings.
###
####################################################################################

# use_build_ignore("R/output/visualization_processes.R", escape = TRUE)

set.seed(05042020)

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

######

n <- 30

locs <- matrix(sort(runif(n)), nrow = n, ncol = 1)
locs <- list(locs1 = locs, locs2 = locs)
locsall <- rbind(locs$locs1, locs$locs2)

covmat1 <- cov_bivariate_expo_latDim(locs = locs, covparms = c(1, 0.1, 0.0))
covmat2 <- cov_bivariate_expo_latDim(locs = locs, covparms = c(1, 0.1, 0.3))
covmat3 <- cov_bivariate_expo_latDim(locs = locs, covparms = c(1, 0.1, 0.6))

covfac1 <- factorize(covmat = covmat1, method = "eigen-I")$covfactor
covfac2 <- factorize(covmat = covmat2, method = "eigen-I")$covfactor
covfac3 <- factorize(covmat = covmat3, method = "eigen-I")$covfactor

y1 <- as.numeric(t(covfac1) %*% rnorm(n * 2)) ; y11 <- y1[1:n] ; y12 <- y1[seq(from = n+1, to = n * 2, by = 1)]
y2 <- as.numeric(t(covfac2) %*% rnorm(n * 2)) ; y21 <- y2[1:n] ; y22 <- y2[seq(from = n+1, to = n * 2, by = 1)]
y3 <- as.numeric(t(covfac3) %*% rnorm(n * 2)) ; y31 <- y3[1:n] ; y32 <- y3[seq(from = n+1, to = n * 2, by = 1)]

x <- locs$locs1[, 1]

plot(x, y11, type = 'o', ylim = c(-3, 27), col = 'red')
points(x, y21, type = 'o', col = 'blue')
points(x, y31, type = 'o', col = 'green')

points(x, y12 + 8, type = 'o', col = 'red')
points(x, y22 + 16, type = 'o', col = 'blue')
points(x, y32 + 24, type = 'o', col = 'green')

for(i in 1:4) abline(h = 8 * (i-1))
for(i in 1:n) abline(v = x[i], lty = "dotted")

diff <- 8
df <- data.frame(x = x, y11 = y11, y12 = y12 + diff, y21 = y21, y22 = y22 + diff * 2, y31 = y31, y32 = y32 + diff * 3)
df <- df %>% tidyr::gather(key = "type", value = "processes", -x)
df$latent.d <- as.character(c(rep(0, n*2), rep(0.3, n*2), rep(0.6, n*2)))

p <- ggplot(data = df, aes(x = x, y = processes, group = type)) + geom_line(aes(color = latent.d)) + geom_point(aes(color = latent.d)) + geom_hline(yintercept = c(0, diff, diff * 2, diff * 3)) + geom_vline(xintercept = x, linetype = "dotted") + theme(axis.ticks = element_blank()) + scale_y_continuous(breaks=c(0, diff, diff*2, diff*3), labels = c("original", "latent.d = 0.0", "latent.d = 0.3", "latent.d = 0.6")) + scale_x_continuous(breaks=NULL)

ggplot2::ggsave("processes.pdf", p, width = 15.2, height = 5.7)


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

### derivative case 2 #############################################################################################

locs <- seq(from = 0, to = 0.5, by = 0.01)
locs <- matrix(locs, nrow = length(locs), ncol = 1)
locs <- cbind(locs, 0)

covparms = c(1, 0.1)

covvec.mat25 <- correlationVecchia::cov_matern_2.5(locs = locs, covparms = covparms)[1, ]
covvec.mat45 <- correlationVecchia::cov_matern_4.5(locs = locs, covparms = covparms)[1, ]
covvec.sqexp <- correlationVecchia::cov_squared_expo(locs = locs, covparms = covparms)[1, ]

corvec.mat25 <- covvec.mat25 / covvec.mat25[1]
corvec.mat45 <- covvec.mat45 / covvec.mat45[1]
corvec.sqexp <- covvec.sqexp / covvec.sqexp[1]

covvec.dgp.mat25 <- correlationVecchia:::.partial_cov_matern_2.5(locs = locs, covparms = covparms, coord = c(1))[1, ]
covvec.dgp.mat45 <- correlationVecchia:::.partial_cov_matern_4.5(locs = locs, covparms = covparms, coord = c(1))[1, ]
covvec.dgp.sqexp <- correlationVecchia:::.partial_cov_squared_expo(locs = locs, covparms = covparms, coord = c(1))[1, ]

covvec.d2gp.mat25 <- correlationVecchia:::.double_partial_cov_matern_2.5(locs = locs, covparms = covparms, coord = c(1))[1, ]
covvec.d2gp.mat45 <- correlationVecchia:::.double_partial_cov_matern_4.5(locs = locs, covparms = covparms, coord = c(1))[1, ]
covvec.d2gp.sqexp <- correlationVecchia:::.double_partial_cov_squared_expo(locs = locs, covparms = covparms, coord = c(1))[1, ]

corvec.dgp.mat25 <- covvec.dgp.mat25 / sqrt(covvec.mat25[1]) / sqrt(covvec.d2gp.mat25[1])
corvec.dgp.mat45 <- covvec.dgp.mat45 / sqrt(covvec.mat45[1]) / sqrt(covvec.d2gp.mat45[1])
corvec.dgp.sqexp <- covvec.dgp.sqexp / sqrt(covvec.sqexp[1]) / sqrt(covvec.d2gp.sqexp[1])

corvec.d2gp.mat25 <- covvec.d2gp.mat25 / covvec.d2gp.mat25[1]
corvec.d2gp.mat45 <- covvec.d2gp.mat45 / covvec.d2gp.mat45[1]
corvec.d2gp.sqexp <- covvec.d2gp.sqexp / covvec.d2gp.sqexp[1]

legvec <- c("cov(GP, GP)    ", "cov(GP', GP)    ", "cov(GP', GP')    ")

par(mfrow = c(3, 2))

plot(locs[, 1], covvec.mat25, type = 'l', lwd = 2, col = 'black', xlab = 'distance', ylab = 'covariance', main = 'Matern covariance with nu = 2.5', ylim = range(c(covvec.mat25, covvec.dgp.mat25, covvec.d2gp.mat25)))
lines(locs[, 1], covvec.dgp.mat25, lwd = 2, col = 'red')
lines(locs[, 1], covvec.d2gp.mat25, lwd = 2, col = 'blue')
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = legvec, lty = 1, lwd = 1, col = c("black", "red", "blue"))

plot(locs[, 1], corvec.mat25, type = 'l', lwd = 2, col = 'black', xlab = 'distance', ylab = 'correlation', main = 'Matern correlation with nu = 2.5', ylim = range(c(corvec.mat25, corvec.dgp.mat25, corvec.d2gp.mat25)))
lines(locs[, 1], corvec.dgp.mat25, lwd = 2, col = 'red')
lines(locs[, 1], corvec.d2gp.mat25, lwd = 2, col = 'blue')
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = legvec, lty = 1, lwd = 1, col = c("black", "red", "blue"))

plot(locs[, 1], covvec.mat45, type = 'l', lwd = 2, col = 'black', xlab = 'distance', ylab = 'covariance', main = 'Matern covariance with nu = 4.5', ylim = range(c(covvec.mat45, covvec.dgp.mat45, covvec.d2gp.mat45)))
lines(locs[, 1], covvec.dgp.mat45, lwd = 2, col = 'red')
lines(locs[, 1], covvec.d2gp.mat45, lwd = 2, col = 'blue')
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = legvec, lty = 1, lwd = 1, col = c("black", "red", "blue"))

plot(locs[, 1], corvec.mat45, type = 'l', lwd = 2, col = 'black', xlab = 'distance', ylab = 'correlation', main = 'Matern correlation with nu = 4.5', ylim = range(c(corvec.mat45, corvec.dgp.mat45, corvec.d2gp.mat45)))
lines(locs[, 1], corvec.dgp.mat45, lwd = 2, col = 'red')
lines(locs[, 1], corvec.d2gp.mat45, lwd = 2, col = 'blue')
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = legvec, lty = 1, lwd = 1, col = c("black", "red", "blue"))

plot(locs[, 1], covvec.sqexp, type = 'l', lwd = 2, col = 'black', xlab = 'distance', ylab = 'covariance', main = 'Squared exponential covariance', ylim = range(c(covvec.sqexp, covvec.dgp.sqexp, covvec.d2gp.sqexp)))
lines(locs[, 1], covvec.dgp.sqexp, lwd = 2, col = 'red')
lines(locs[, 1], covvec.d2gp.sqexp, lwd = 2, col = 'blue')
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = legvec, lty = 1, lwd = 1, col = c("black", "red", "blue"))

plot(locs[, 1], corvec.sqexp, type = 'l', lwd = 2, col = 'black', xlab = 'distance', ylab = 'correlation', main = 'Squared exponential correlation', ylim = range(c(corvec.sqexp, corvec.dgp.sqexp, corvec.d2gp.sqexp)))
lines(locs[, 1], corvec.dgp.sqexp, lwd = 2, col = 'red')
lines(locs[, 1], corvec.d2gp.sqexp, lwd = 2, col = 'blue')
abline(h = 0, v = 0, col = 'gray', lwd = 2)
legend("topright", legend = legvec, lty = 1, lwd = 1, col = c("black", "red", "blue"))

par(mfrow = c(1, 1))

df1 <- data.frame(distance = locs[, 1], covgp = covvec.sqexp, covdgp = covvec.dgp.sqexp, cov.d2gp = covvec.d2gp.sqexp)
df1 <- df1 %>% tidyr::gather(key = "type", value = "covariance", -distance)


df2 <- data.frame(distance = locs[, 1], corgp = corvec.sqexp, cordgp = corvec.dgp.sqexp, cor.d2gp = corvec.d2gp.sqexp)
df2 <- df2 %>% tidyr::gather(key = "type", value = "correlation", -distance)

p <- list()
p[[1]] <- ggplot(data = df1, aes(x = distance, y = covariance, group = type)) + geom_line(aes(color = type), size = 1.1) + scale_color_manual(values=c("#4DAF4A", "#377EB8", "#E41A1C"), labels = c("cov(y', y')", "cov(y', y)", "cov(y, y)")) + xlim(0, 0.5) + ylim(-50, 100) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme(legend.title = element_blank())
p[[2]] <- ggplot(data = df2, aes(x = distance, y = correlation, group = type)) + geom_line(aes(color = type), size = 1.1) + scale_color_manual(values=c("#4DAF4A", "#377EB8", "#E41A1C"), labels = c("corr(y', y')", "corr(y', y)", "corr(y, y)")) + xlim(0, 0.5) + ylim(-0.5, 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme(legend.title = element_blank())

p <- ggpubr::ggarrange(p[[1]], p[[2]], nrow = 2, ncol = 1)

ggplot2::ggsave("sqexp.pdf", p, width = 15.2, height = 8)

### derivative case 3 #############################################################################################

n <- 20
d <- 2

locs <- matrix(runif(n * d), nrow = n, ncol = d)
locsall <- rbind(locs, locs, locs)
locs <- list(locs1 = locs, locs2 = locs, locs3 = locs)

covmat <- corr_derivative_squared_expo_2d(locs = locs, covparms = c(1, 0.1))
fields::image.plot(abs(covmat))

m <- 6

out.b1 <- baseline_1_multivariate_specify(locs = locs, m = m)
ord.b1 <- out.b1$ord
con.b1 <- out.b1$U.prep$revNNarray[, seq(from = m, to = 1, by = -1)]

out.b2 <- baseline_2_multivariate_specify(locs = locs, m = m)
ord.b2 <- out.b2$ord
con.b2 <- out.b2$U.prep$revNNarray[, seq(from = m, to = 1, by = -1)]

out.b3 <- baseline_3_multivariate_specify(locs = locs, m = m)
ord.b3 <- out.b3$ord
con.b3 <- out.b3$U.prep$revNNarray[, seq(from = m, to = 1, by = -1)]

out.b4 <- baseline_4_multivariate_specify(locs = locs, m = m, abs.corr = TRUE, covmodel = cov_derivative_squared_expo_2d, covparms = c(1, 0.1))
ord.b4 <- out.b4$ord
con.b4 <- out.b4$U.prep$revNNarray[, seq(from = m, to = 1, by = -1)]

out.cc <- corrvecchia_specify_knownCovparms(locs = locsall, m = m, ordering = "maxmin", ordering.method = "correlation", abs.corr = TRUE, conditioning = "NN", conditioning.method = "correlation", covmodel = covmat, covparms = c(1))
ord.cc <- out.b4$ord
con.cc <- out.b4$U.prep$revNNarray[, seq(from = m, to = 1, by = -1)]

################################################################################

ind <- 2 * n + 1

################################################################################

par(mfrow = c(1, 2))

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Ordering (Baseline 1)')

locsord <- locsall[ord.b1, ]

for(i in 1:n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 4)
for(i in 1:n+n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 3)
for(i in 1:n+2*n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 2)

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Conditioning (Baseline 1)')
points(locsord[ind, 1], locsord[ind, 2], pch = 15, cex = 2)

for(i in 1:m) {
  k <- con.b1[ind, i]
  points(locsord[k, 1], locsord[k, 2], pch = 13)
  text(locsord[k, 1], locsord[k, 2], label = i, pos = 4)
}

par(mfrow = c(1, 1))

################################################################################

par(mfrow = c(1, 2))

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Ordering (Baseline 2)')

locsord <- locsall[ord.b2, ]

for(i in 1:n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 4)
for(i in 1:n+n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 3)
for(i in 1:n+2*n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 2)

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Conditioning (Baseline 2)')
points(locsord[ind, 1], locsord[ind, 2], pch = 15, cex = 2)

for(i in 1:m) {
  k <- con.b2[ind, i]
  points(locsord[k, 1], locsord[k, 2], pch = 13)
  text(locsord[k, 1], locsord[k, 2], label = i, pos = 4)
}

par(mfrow = c(1, 1))

###############################################################################

par(mfrow = c(1, 2))

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Ordering (Baseline 3)')

locsord <- locsall[ord.b3, ]

for(i in 1:n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 4)
for(i in 1:n+n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 3)
for(i in 1:n+2*n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 2)

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Conditioning (Baseline 3)')
points(locsord[ind, 1], locsord[ind, 2], pch = 15, cex = 2)

for(i in 1:m) {
  k <- con.b3[ind, i]
  points(locsord[k, 1], locsord[k, 2], pch = 13)
  text(locsord[k, 1], locsord[k, 2], label = i, pos = 4)
}

par(mfrow = c(1, 1))

#############################################################################

par(mfrow = c(1, 2))

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Ordering (Baseline 4)')

locsord <- locsall[ord.b4, ]

for(i in 1:n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 4)
for(i in 1:n+n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 3)
for(i in 1:n+2*n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 2)

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Conditioning (Baseline 4)')
points(locsord[ind, 1], locsord[ind, 2], pch = 15, cex = 2)

for(i in 1:m) {
  k <- con.b4[ind, i]
  points(locsord[k, 1], locsord[k, 2], pch = 13)
  text(locsord[k, 1], locsord[k, 2], label = i, pos = 4)
}

par(mfrow = c(1, 1))

############################################################################

par(mfrow = c(1, 2))

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Ordering (CorrVecchia)')

locsord <- locsall[ord.cc, ]

for(i in 1:n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 4)
for(i in 1:n+n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 3)
for(i in 1:n+2*n) text(locsord[i, 1], locsord[i, 2], label = i, pos = 2)

plot(locs[[1]], xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2', main = 'Conditioning (CorrVecchia)')
points(locsord[ind, 1], locsord[ind, 2], pch = 15, cex = 2)

for(i in 1:m) {
  k <- con.cc[ind, i]
  points(locsord[k, 1], locsord[k, 2], pch = 13)
  text(locsord[k, 1], locsord[k, 2], label = i, pos = 4)
}

par(mfrow = c(1, 1))

### periodic case #############################################################################################

locs <- cbind(seq(from = 0, to = sqrt(2), by = 0.01), 0)

covmat <- cov_wave(locs = locs, covparms = c(1, 0.2), method = "Dampedsine")
dampedsine <- covmat[, 1]
df1 <- data.frame(locs = locs[, 1], covmat = dampedsine, type = "Damped Sine")

covmat <- cov_wave(locs = locs, covparms = c(1, 0.1, 0.2), method = "Dampedcosine")
dampedcosine <- covmat[, 1]
df2 <- data.frame(locs = locs[, 1], covmat = dampedcosine, type = "Damped Cosine")

covmat <- cov_wave(locs = locs, covparms = c(1, 0, 0.2), method = "BesselJ")
besselJ <- covmat[, 1]
df3 <- data.frame(locs = locs[, 1], covmat = besselJ, type = "Bessel J")

df <- rbind(df1, df2, df3) ; rm(dampedsine, dampedcosine, besselJ)

ggplot(data = df, aes(x = locs, y = covmat, group = type)) + geom_line(aes(color = type)) + scale_color_manual(values=c("#984EA3", "#4DAF4A", "#377EB8"))

df3 <- df3 %>% mutate(distance = sqrt(1-abs(covmat))) %>% select(-type) %>% rename(edist = locs, "periodic correlation" = covmat, "correlation-based distance" = distance)
df3 <- df3 %>% tidyr::gather(key = "type", value = "value", -edist)

p <- ggplot(data = df3, aes(x = edist, y = value, group = type)) + geom_line(aes(color= type), size = 1.1) + scale_color_manual(values=c("#377EB8", "#E41A1C")) + scale_linetype_manual(values=c("dotted", "solid")) + xlab("Euclidean-based distance") + geom_hline(yintercept = c(0, 1)) + geom_vline(xintercept = 0)

ggplot2::ggsave("besselJ.pdf", p, width = 15.2, height = 5.7)
