##########################################################################
#####
##### Author: Myeongjong Kang (kmj.stat@gmail.com)
#####
##### Description: 
#####
##########################################################################

rm(list = ls()) ; gc()

## To employ the Vecchia approximation 
library(GPvecchia) ; library(fields)
## To visualize results
library(tidyr) ; library(ggplot2) ; library(RColorBrewer)

## Call the vecchia_specify_adjusted function
source("1_Pilot_Study/2_vecchia_specify_adjusted.R")
## Call the function to compute multivariate normal KL divergence between true distribution N(mu0,covmat0) and approx distr N(mu1,covmat1)
source("1_Pilot_Study/3_kldiv.R") 
## Call the simulation function
source("1_Pilot_Study/4_simulation.R") 

set.seed(20190814)


##########################################################################
##### Simulation 1: Illustraion of the Vecchia approximation with grid locations
##########################################################################

sim1 <- simulation(ordering = "maxmin", which.coord = NULL, conditioning = "NN", 
                   m = 10, scale.aniso = 10, n = 15^2, 
                   range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

### Visualize the process
# y <- as.numeric(t(chol(sim1$Sigma)) %*% rnorm(sim1$n))
# quilt.plot(sim1$locs[,1], sim1$locs[,2], y)

### Visualize the true covariance matrix and its estimated matrices
# par(mfrow=c(1,3))
# image.plot(sim1$Sigma)
# for(i in 1:2) image.plot(sim1$Sigma.hat[[i]])
# par(mfrow=c(1,1))

### illustration of ordering and conditioning sets
n       <- sim1$n
approx  <- sim1$approx

par(mfrow=c(1,2))
for(i in 1:2){ # i=1: i=2: 
  ## Conditioning set of the n-th observation
  quilt.plot(approx[[i]]$locsord, floor(log(1:n)), 
             col=heat.colors(99), xaxs="r", yaxs='r', 
             main = "Cond. set of the n-th obs.")
  points(approx[[i]]$locsord[n,1], 
         approx[[i]]$locsord[n,2], cex=2, pch='x')
  condlocs <- approx[[i]]$locsord[approx[[i]]$U.prep$revNNarray[n,], ]
  points(condlocs[, 1], condlocs[, 2], cex=1.5)
  
  ## Conditioning set of the (n/2)-th observation
  quilt.plot(approx[[i]]$locsord, floor(log(1:n)), 
             col=heat.colors(99), xaxs="r", yaxs='r', 
             main = "Cond. set of the (n/2)-th obs.")
  points(approx[[i]]$locsord[as.integer(n/2),1], 
         approx[[i]]$locsord[as.integer(n/2),2], cex=2, pch='x')
  condlocs <- approx[[i]]$locsord[approx[[i]]$U.prep$revNNarray[as.integer(n/2), ],]
  points(condlocs[, 1], condlocs[, 2], cex=1.5)
} ; par(mfrow=c(1,1))


##########################################################################
##### Simulation 2: Comparison using KL divergence with grid locations
##########################################################################

cand.m      <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.scale  <- c(1, 5, 10, 15, 20, 25)

sim.maxmin  <- repeated_simulation(ordering = "maxmin", which.coord = NULL, conditioning = "NN", 
                                   ms = cand.m, scales = cand.scale, 
                                   n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

sim.xcoord  <- repeated_simulation(ordering = "coord", which.coord = 1, conditioning = "NN", 
                                   ms = cand.m, scales = cand.scale, 
                                   n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

sim.ycoord  <- repeated_simulation(ordering = "coord", which.coord = 2, conditioning = "NN", 
                                   ms = cand.m, scales = cand.scale, 
                                   n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

# sim.xycoord <- repeated_simulation(ordering = "coord", which.coord = c(1, 2), conditioning = "NN", ms = cand.m, scales = cand.scale, n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

##### KL divergence against m #####

kls.maxmin.euclidean    <- sim.maxmin$kls[[1]][,3]
kls.maxmin.corr         <- sim.maxmin$kls[[2]][,3]
kls.xcoord.euclidean    <- sim.xcoord$kls[[1]][,3]
# kls.xcoord.corr         <- sim.xcoord$kls[[2]][,3] # non-sense
kls.ycoord.euclidean    <- sim.ycoord$kls[[1]][,3]
# kls.ycoord.corr         <- sim.ycoord$kls[[2]][,3] # non-sense

vis.list  <- c(kls.maxmin.euclidean, kls.maxmin.corr, 
               kls.xcoord.euclidean, kls.ycoord.euclidean)
min.kls   <- min(vis.list)
max.kls   <- max(vis.list)

plot(cand.m, log10(kls.maxmin.euclidean), type = "l", col = 1, lty = 1, lwd = 3,
     ylim = c(log10(min.kls), log10(max.kls)), 
     xlab = "m", ylab = "log10(KL)", main = NULL, cex.lab=1.5)
lines(cand.m, log10(kls.maxmin.corr), type = "l", col = 2, lty = 2, lwd = 3)
lines(cand.m, log10(kls.xcoord.euclidean), type = "l", col = 3, lty = 3, lwd = 3)
lines(cand.m, log10(kls.ycoord.euclidean), type = "l", col = 4, lty = 4, lwd = 3)
legend("bottomleft", 
       legend=c("maxmin + Euclidean", "maxmin + correlation",
                "x-coord + Euclidean", "y-coord + Euclidean"),
       col=1:5, lty=1:5, lwd = 3, cex=1)

##### KL divergence against scale.aniso when m = 10 #####

kls.maxmin.euclidean    <- sim.maxmin$kls[[1]][3,]
kls.maxmin.corr         <- sim.maxmin$kls[[2]][3,]
kls.xcoord.euclidean    <- sim.xcoord$kls[[1]][3,]
# kls.xcoord.corr         <- sim.xcoord$kls[[2]][3,] # non-sense
kls.ycoord.euclidean    <- sim.ycoord$kls[[1]][3,]
# kls.ycoord.corr         <- sim.ycoord$kls[[2]][3,] # non-sense

vis.list  <- c(kls.maxmin.euclidean, kls.maxmin.corr, 
               kls.xcoord.euclidean, kls.ycoord.euclidean)
min.kls   <- min(vis.list)
max.kls   <- max(vis.list)

plot(cand.scale, log10(kls.maxmin.euclidean), type = "l", col = 1, lty = 1, lwd = 3,
     ylim = c(log10(min.kls), log10(max.kls)), 
     xlab = "m", ylab = "log10(KL)", main = NULL, cex.lab=1.5)
lines(cand.scale, log10(kls.maxmin.corr), type = "l", col = 2, lty = 2, lwd = 3)
lines(cand.scale, log10(kls.xcoord.euclidean), type = "l", col = 3, lty = 3, lwd = 3)
lines(cand.scale, log10(kls.ycoord.euclidean), type = "l", col = 4, lty = 4, lwd = 3)
legend("bottomleft", 
       legend=c("maxmin + Euclidean", "maxmin + correlation",
                "x-coord + Euclidean", "y-coord + Euclidean"),
       col=1:5, lty=1:5, lwd = 3, cex=1)

##### KL divergence against scale.aniso when m = 30 #####

kls.maxmin.euclidean    <- sim.maxmin$kls[[1]][7,]
kls.maxmin.corr         <- sim.maxmin$kls[[2]][7,]
kls.xcoord.euclidean    <- sim.xcoord$kls[[1]][7,]
# kls.xcoord.corr         <- sim.xcoord$kls[[2]][7,] # non-sense
kls.ycoord.euclidean    <- sim.ycoord$kls[[1]][7,]
# kls.ycoord.corr         <- sim.ycoord$kls[[2]][7,] # non-sense

vis.list  <- c(kls.maxmin.euclidean, kls.maxmin.corr, 
               kls.xcoord.euclidean, kls.ycoord.euclidean)
min.kls   <- min(vis.list)
max.kls   <- max(vis.list)

plot(cand.scale, log10(kls.maxmin.euclidean), type = "l", col = 1, lty = 1, lwd = 3,
     ylim = c(log10(min.kls), log10(max.kls)), 
     xlab = "m", ylab = "log10(KL)", main = NULL, cex.lab=1.5)
lines(cand.scale, log10(kls.maxmin.corr), type = "l", col = 2, lty = 2, lwd = 3)
lines(cand.scale, log10(kls.xcoord.euclidean), type = "l", col = 3, lty = 3, lwd = 3)
lines(cand.scale, log10(kls.ycoord.euclidean), type = "l", col = 4, lty = 4, lwd = 3)
legend("bottomleft", 
       legend=c("maxmin + Euclidean", "maxmin + correlation",
                "x-coord + Euclidean", "y-coord + Euclidean"),
       col=1:5, lty=1:5, lwd = 3, cex=1)


##########################################################################
##### Simulation 3: Illustraion of the Vecchia approximation with non-grid locations
##########################################################################

sim3  <- simulation_nongrid(ordering = "maxmin", which.coord = NULL, conditioning = "NN", 
                            m = 10, scale.aniso = 10, 
                            n = 15^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

### Visualize the process
# y <- as.numeric(t(chol(sim3$Sigma)) %*% rnorm(sim3$n))
# quilt.plot(sim3$locs[,1], sim3$locs[,2], y)

### Visualize the true covariance matrix and its estimated matrices
# par(mfrow=c(1,3))
# image.plot(sim3$Sigma)
# for(i in 1:2) image.plot(sim3$Sigma.hat[[i]])
# par(mfrow=c(1,1))

## illustration of ordering and conditioning sets
n       <- sim3$n
approx  <- sim3$approx

par(mfrow=c(1,2))
for(i in 1:2){ # i=1: i=2: 
  ## Conditioning set of the n-th observation
  quilt.plot(approx[[i]]$locsord, floor(log(1:n)), 
             col=heat.colors(99), xaxs="r", yaxs='r', 
             main = "Cond. set of the n-th obs.")
  points(approx[[i]]$locsord[n,1], 
         approx[[i]]$locsord[n,2], cex=2, pch='x')
  condlocs <- approx[[i]]$locsord[approx[[i]]$U.prep$revNNarray[n,], ]
  points(condlocs[, 1], condlocs[, 2], cex=1.5)
  
  ## Conditioning set of the (n/2)-th observation
  quilt.plot(approx[[i]]$locsord, floor(log(1:n)), 
             col=heat.colors(99), xaxs="r", yaxs='r', 
             main = "Cond. set of the (n/2)-th obs.")
  points(approx[[i]]$locsord[as.integer(n/2),1], 
         approx[[i]]$locsord[as.integer(n/2),2], cex=2, pch='x')
  condlocs <- approx[[i]]$locsord[approx[[i]]$U.prep$revNNarray[as.integer(n/2), ],]
  points(condlocs[, 1], condlocs[, 2], cex=1.5)
} ; par(mfrow=c(1,1))


##########################################################################
##### Simulation 4: Comparison using KL divergence with non-grid locations
##########################################################################

cand.m      <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)
cand.scale  <- c(1, 5, 10, 15, 20, 25)

sim.maxmin.nongrid  <- repeated_simulation_nongrid(ordering = "maxmin", which.coord = NULL, conditioning = "NN", 
                                                   ms = cand.m, scales = cand.scale, 
                                                   n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

sim.xcoord.nongrid  <- repeated_simulation_nongrid(ordering = "coord", which.coord = 1, conditioning = "NN", 
                                                   ms = cand.m, scales = cand.scale, 
                                                   n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

sim.ycoord.nongrid  <- repeated_simulation_nongrid(ordering = "coord", which.coord = 2, conditioning = "NN", 
                                                   ms = cand.m, scales = cand.scale, 
                                                   n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

# sim.xycoord.nongrid <- repeated_simulation_nongrid(ordering = "coord", which.coord = c(1, 2), conditioning = "NN", ms = cand.m, scales = cand.scale, n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

##### KL divergence against m #####

kls.maxmin.euclidean.nongrid    <- sim.maxmin.nongrid$kls[[1]][,3]
kls.maxmin.corr.nongrid         <- sim.maxmin.nongrid$kls[[2]][,3]
kls.xcoord.euclidean.nongrid    <- sim.xcoord.nongrid$kls[[1]][,3]
# kls.xcoord.corr.nongrid         <- sim.xcoord.nongrid$kls[[2]][,3] # non-sense
kls.ycoord.euclidean.nongrid    <- sim.ycoord.nongrid$kls[[1]][,3]
# kls.ycoord.corr.nongrid         <- sim.ycoord.nongrid$kls[[2]][,3] # non-sense

vis.dat <- data.frame(kls.maxmin.euclidean.nongrid, kls.maxmin.corr.nongrid, kls.xcoord.euclidean.nongrid, kls.ycoord.euclidean.nongrid)
vis.dat <- vis.dat[, order(colnames(vis.dat))]
vis.dat <- cbind(rep(cand.m, times = 4), gather(vis.dat))
colnames(vis.dat) <- c("m", "method", "KL")
head(vis.dat)

unique(vis.dat[,"method"])
# unique(vis.dat[,"method"])[order(unique(vis.dat[,"method"]))]

kls.legend <- c("Correlation + Maxmin         ", "Euclidean + Maxmin", "Euclidean + x-coord", "Euclidean + y-coord")

size <- 28
p <- ggplot(vis.dat, aes(x=m, y = log10(KL), col = method)) + 
     geom_point(aes(shape = method), size = 5) + geom_line(size = 2, alpha = 0.7) +
     ylab('log10(KL)') + scale_x_discrete(name = 'm', limits=cand.m, labels=as.character(cand.m)) +
     scale_color_manual(values = brewer.pal(4, "Set1"), labels = kls.legend) +
     scale_shape_manual(values = c(15, 16, 17, 18), labels = kls.legend) +
     theme(axis.title.x = element_text(size = size), 
           axis.text.x = element_text(size = size-4),
           axis.title.y = element_text(size = size), 
           axis.text.y = element_text(size = size-4),
           legend.position = 'top',
           legend.title = element_blank(),
           legend.direction = 'horizontal',
           legend.text=element_text(size = size)) +
     guides(color=guide_legend(nrow=2,byrow=FALSE))
print(p)

##### KL divergence against scale.aniso when m = 30 #####

kls.maxmin.euclidean.nongrid    <- sim.maxmin.nongrid$kls[[1]][7,]
kls.maxmin.corr.nongrid         <- sim.maxmin.nongrid$kls[[2]][7,]
kls.xcoord.euclidean.nongrid    <- sim.xcoord.nongrid$kls[[1]][7,]
# kls.xcoord.corr.nongrid         <- sim.xcoord.nongrid$kls[[2]][7,] # non-sense
kls.ycoord.euclidean.nongrid    <- sim.ycoord.nongrid$kls[[1]][7,]
# kls.ycoord.corr.nongrid         <- sim.ycoord.nongrid$kls[[2]][7,] # non-sense
kls.xycoord.euclidean.nongrid   <- sim.xycoord.nongrid$kls[[1]][7,]
# kls.xycoord.corr.nongrid        <- sim.xycoord.nongrid$kls[[2]][7,] # non-sense

vis.dat <- data.frame(kls.maxmin.euclidean.nongrid, kls.maxmin.corr.nongrid, kls.xcoord.euclidean.nongrid, kls.ycoord.euclidean.nongrid)
vis.dat <- vis.dat[, order(colnames(vis.dat))]
vis.dat <- cbind(rep(cand.scale, times = 4), gather(vis.dat))
colnames(vis.dat) <- c("scale", "method", "KL")
head(vis.dat)

kls.legend <- c("Correlation + Maxmin         ", "Euclidean + Maxmin", "Euclidean + x-coord", "Euclidean + y-coord")

size <- 28
p <- ggplot(vis.dat, aes(x=scale, y = log10(KL), col = method)) + 
  geom_point(aes(shape = method), size = 5) + geom_line(size = 2, alpha = 0.7) +
  ylab('log10(KL)') + scale_x_discrete(name = 'a', limits=cand.scale, labels=as.character(cand.scale)) +
  scale_color_manual(values = brewer.pal(4, "Set1"), labels = kls.legend) +
  scale_shape_manual(values = c(15, 16, 17, 18), labels = kls.legend) +
  theme(axis.title.x = element_text(size = size), 
        axis.text.x = element_text(size = size-4),
        axis.title.y = element_text(size = size), 
        axis.text.y = element_text(size = size-4),
        legend.position = 'top',
        legend.title = element_blank(),
        legend.direction = 'horizontal',
        legend.text=element_text(size = size)) +
  guides(color=guide_legend(nrow=2,byrow=FALSE))
print(p)
