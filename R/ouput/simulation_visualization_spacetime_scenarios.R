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

####################################################################################

alpha       <- 0.4
size.point  <- 0.4
size.line   <- 1
size.lab    <- 18
size.text   <- 16
size.margin <- c(5.5, 20, 5.5, 5.5)

scens <- c("random", "station", "gridded", "satellite")

####################################################################################

out1  <- generate_gp_spacetime(nsim = 1, n = 30^2, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0))
out2  <- generate_gp_spacetime(nsim = 1, n = 10^2, d = 2, t.len = 9, method.locs = "space.random.time.grid", covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0))
out3  <- generate_gp_spacetime(nsim = 1, n = 10^2, d = 2, t.len = 9, method.locs = "all.grid", covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0))
out4  <- generate_gp_spacetime(nsim = 1, n = 30^2, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_expo_spacetime, covparms = c(1, 0.1, 1.0))

locs <- list()
locs[[1]] <- out1$sim$sim1$locs %>% cbind(seq(900)) %>% as.data.frame()
locs[[2]] <- out2$sim$sim1$locs %>% cbind(seq(900)) %>% as.data.frame()
locs[[3]] <- out3$sim$sim1$locs %>% cbind(seq(900)) %>% as.data.frame()
locs[[4]] <- out4$sim$sim1$locs %>% cbind(seq(900)) %>% as.data.frame()

for(i in 1:4) colnames(locs[[i]]) <- c("x", "y", "t", "index")

k <- 1
p1 <- list()
for(i in 1:4) {

  locs[[i]]$scen <- scens[i]

  p1[[i]] <- locs[[i]] %>% ggplot() +
    geom_point(aes(x = x, y = y), size = size.point, alpha = alpha) + xlab("x1") + ylab("x2") + # labs(title = "Space domain") +
    theme(axis.title.x = element_text(size = size.lab),
          axis.text.x = element_text(size = size.text),
          axis.title.y = element_text(size = size.lab),
          axis.text.y = element_text(size = size.text),
          strip.text = element_text(size = size.lab),
          plot.margin = unit(size.margin, "pt")) +
    facet_wrap(.~scen, scales = "free")
  k <- k + 1

  p1[[i+4]] <- locs[[i]] %>% ggplot() +
    geom_point(aes(x = index, y = t), size = size.point, alpha = alpha) + xlab("index") + ylab("time") + # labs(title = "Time") +
    theme(axis.title.x = element_text(size = size.lab),
          axis.text.x = element_text(size = size.text),
          axis.title.y = element_text(size = size.lab),
          axis.text.y = element_text(size = size.text),
          strip.text = element_text(size = size.lab),
          plot.margin = unit(size.margin, "pt"))
  facet_wrap(.~scen, scales = "free")
  k <- k + 1
}

# vis.all   <- arrangeGrob(grobs = p1, nrow = 2, ncol = 4)
# p1        <- grid.arrange(vis.all)

p1        <- list(plot1 = p1[[1]], plot2 = p1[[2]], plot3 = p1[[3]], plot4 = p1[[4]])
vis.all   <- arrangeGrob(grobs = p1, nrow = 1, ncol = 4)
p1        <- grid.arrange(vis.all)

ggplot2::ggsave("simulation_spacetime_scenarios.pdf", p1, width = 15.2, height = 3.7)



