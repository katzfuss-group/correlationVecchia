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

cand.r  <- c(1, 2, 3, 4, 5, 6)

locs <- matrix(runif(20), 10, 2)

p1        <- list()
mypal     <- brewer.pal(10, "RdYlBu")
mypal     <- mypal[seq(from = length(mypal), to = 1, by = -1)]
for(i in 1:length(cand.r)) {
  covmat <- cov_derivative_matern_2.5_2d(locs = locs, covparms = c(1, cand.r[i]))
  covmat <- 1 - (covmat - min(covmat))/(max(covmat) - min(covmat))
  p1[[i]]   <- covmat %>% melt() %>% ggplot(aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) + theme_void() + theme(legend.position="none", plot.title = element_text(hjust=0.5)) + ggtitle(label = paste0("Distance matrix with range parameter = ", cand.r[i])) + scale_fill_gradientn(colours = mypal, limits = c(0,1))
}

vis.all   <- arrangeGrob(grobs = p1, nrow = 2, ncol = 3)
leg       <- get_legend(covmat %>% melt() %>% ggplot(aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) +  scale_fill_gradientn(colours = mypal, limits = c(0,1)))

p1        <- grid.arrange(vis.all, leg, nrow = 1, ncol = 2, widths = c(6, 1), heights = c(1))

ggplot2::ggsave("p1_deriv.pdf", p1, width = 15.2, height = 5.7)

