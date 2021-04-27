####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of Vecchia-like approximation in various settings.
###
####################################################################################

# use_build_ignore("R/output/visualization_processes.R", escape = TRUE)

set.seed(04272021)

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



### spacetime case #############################################################################################

out1  <- generate_gp_spacetime(nsim = 1, n = 30^2, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out2  <- generate_gp_spacetime(nsim = 1, n = 5^2, d = 2, t.len = 36, method.locs = "space.random.time.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out3  <- generate_gp_spacetime(nsim = 1, n = 5^2, d = 2, t.len = 36, method.locs = "all.grid", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))
out4  <- generate_gp_spacetime(nsim = 1, n = 30^2, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_spacetime_expo, covparms = c(1, 0.75, 50, 25))

locs <- list()
locs[[1]] <- out1$sim$sim1$locs %>% cbind(seq(900)) %>% as.data.frame()
locs[[2]] <- out2$sim$sim1$locs %>% cbind(seq(900)) %>% as.data.frame()
locs[[3]] <- out3$sim$sim1$locs %>% cbind(seq(900)) %>% as.data.frame()
locs[[4]] <- out4$sim$sim1$locs %>% cbind(seq(900)) %>% as.data.frame()

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

ggplot2::ggsave("simulation_spacetime_locations.pdf", p1, width = 15.2, height = 5.7)

### derivative case #############################################################################################


