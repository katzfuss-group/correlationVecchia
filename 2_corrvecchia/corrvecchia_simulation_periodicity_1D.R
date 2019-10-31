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
source("2_corrvecchia/corrvecchia.R")
source("2_corrvecchia/kldiv.R")

set.seed(10102019)


####################################################################
#### Wave covariance model (Risser MD, Calder CA (2015))
####################################################################

# cov.aniso <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[,2])) / covparms[2])

# covparms = c(sigma, period)
covmodel_dampedsine <- function(locs, covparms) {
  h       <- fields::rdist(locs)
  
  ind     <- which(h < 1e-8, arr.ind = T)
  c       <- covparms[1] * sin(h / covparms[2])* (covparms[2] / h)
  c[ind]  <- covparms[1]
  
  c
}

# covparms = c(sigma, range, period)
covmodel_dampedcosine <- function(locs, covparms) {
  h <- fields::rdist(locs)
  
  covparms[1] * exp(- h / covparms[2]) * cos(h / covparms[3])
}

# covparms = c(sigma, nu, period)
covmodel_besselJ <- function(locs, covparms) {
  h <- fields::rdist(locs)
  c <- gamma(covparms[2] + 1) * (2*covparms[3]/c(h))^covparms[2] * Bessel::BesselJ(c(h) / covparms[3], covparms[2])
  c[which(is.nan(c))] <- 1
  
  covparms[1] * matrix(c, nrow = nrow(h), ncol = ncol(h))
}

covmodel_wave <- function(locs, covparms, covm = "Dampedsine") {
  if(covm == "Dampedsine") {
    return(covmodel_dampedsine(locs, covparms))
  } else if(covm == "Dampedcosine") {
    return(covmodel_dampedcosine(locs, covparms))
  } else if(covm == "BesselJ") {
    return(covmodel_besselJ(locs, covparms))
  } else {
    stop("Select one of the following covariance models: Dampedsine, Dampedcosine, or besselJ.")
  }
}

# n     <- 15^2
# locs  <- matrix(runif(n, 0, 1), n, 1)
# h     <- fields::rdist(locs)
# ind   <- order(unlist(h))
# hord  <- unlist(h)[ind]
# 
# covparms  <- c(1, 1/10)
# covh      <- covmodel_wave(locs = locs, covparms = covparms, covm = "Dampedsine")
# covhord   <- unlist(covh)[ind]
# plot(hord, covhord, type = 'l', ylim = c(-0.8, 1), col = 1)
# abline(h = 0, col = 'gray')
# 
# covparms  <- c(1, 1, 1/10)
# covh      <- covmodel_wave(locs = locs, covparms = covparms, covm = "Dampedcosine")
# covhord   <- unlist(covh)[ind]
# lines(hord, covhord, type = 'l', col = 2)
# 
# covparms  <- c(1, 1, 1/20)
# covh      <- covmodel_wave(locs = locs, covparms = covparms, covm = "Dampedcosine")
# covhord   <- unlist(covh)[ind]
# lines(hord, covhord, type = 'l', col = 3)
# 
# covparms  <- c(1, 1/10, 1/10)
# covh      <- covmodel_wave(locs = locs, covparms = covparms, covm = "BesselJ")
# covhord   <- unlist(covh)[ind]
# lines(hord, covhord, type = 'l', col = 4)
# 
# legend("topright", legend = c("Damped Sine with (1, 0.1)", "Damped Cosine with (1, 1, 0.1)", "Damped Cosine with (1, 1, 0.05)","Bessel J with (1, 0.1, 0.1)"), col = 1:4, lty = 1)


####################################################################
#### simulation function
####################################################################

positive_def <- function(Sigma, tol){
  eig.decomp  <- eigen(Sigma)
  diagvec     <- ifelse(eig.decomp$values < tol, tol, eig.decomp$values)
  
  Sigma.modified <- eig.decomp$vectors %*% diag(diagvec) %*% t(eig.decomp$vectors)
  return(Sigma.modified)
}

simulation <- function(n = 15^2, m = 10, covparms = c(1, 1/10), covm = "Dampedsine", def.dist = NULL, tol = 1e-6) {
  
  locs      <- matrix(runif(n, 0, 1), n, 1)
  Sigma     <- covmodel_wave(locs = locs, covparms = covparms, covm = covm)
  
  Sigma.modified <- positive_def(Sigma, tol)
  
  y         <- as.numeric(t(chol(Sigma.modified)) %*% rnorm(n))
  
  ### specify vecchia approximations
  approx <- list()
  
  # standard vecchia with maxmin ordering
  approx[[1]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
  # standard vecchia with x coord ordering
  approx[[2]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "coord", which.coord = 1, cond.yz='y', conditioning = "NN")
  # standard vecchia with y coord ordering
  approx[[3]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "coord", which.coord = 2, cond.yz='y', conditioning = "NN")
  # euclidean-based ordering + euclidean-based NN conditioning
  approx[[4]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = def.dist, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma.modified, covparms = covparms)
  # euclidean-based ordering + correlation-based NN conditioning
  approx[[5]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = def.dist, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = covparms)
  # correlation-based ordering + euclidean-based NN conditioning
  approx[[6]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = def.dist, ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma.modified, covparms = covparms)
  # correlation-based ordering + correlation-based NN conditioning
  approx[[7]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = def.dist, ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = covparms)
  # x-coordinate-based ordering + correlation-based NN conditioning
  approx[[8]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(1), def.dist = def.dist, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = covparms)
  # y-coordinate-based ordering + correlation-based NN conditioning
  approx[[9]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(2), def.dist = def.dist, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = covparms)
  # (x+y)-coordinate-based ordering + correlation-based NN conditioning
  approx[[10]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(1, 2), def.dist = def.dist, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = covparms)
  # (x+y)-coordinate-based ordering + euclidean-based NN conditioning
  approx[[11]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(1, 2), def.dist = def.dist, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma.modified, covparms = covparms)
  # standard vecchia with maxmin ordering
  approx[[12]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "coord", which.coord = c(1, 2), cond.yz='y', conditioning = "NN")
  
  ### compute approximate covariance matrices
  n.approx    <- length(approx)
  Sigma.hat   <- list()
  kls         <- rep(NA, n.approx)
  for(i in 1:n.approx){
    
    Sigma.ord       <- covmodel_wave(locs = approx[[i]]$locsord, covparms = covparms, covm = covm) # true cov in appropriate ordering
    
    Sigma.ord.modified <- positive_def(Sigma.ord, tol)
    
    U               <- createU(approx[[i]], covparms, 0, covmodel = Sigma.ord.modified)$U
    revord          <- order(approx[[i]]$ord)
    Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    kls[i]          <- kldiv(Sigma.modified, Sigma.hat[[i]])
  }
  
  result                  <- list()
  result$n                <- n
  result$m                <- m
  result$covparms         <- covparms
  result$covm             <- covm
  result$locs             <- locs
  result$Sigma            <- Sigma
  result$Sigma.modified   <- Sigma.modified
  result$approx           <- approx
  result$Sigma.hat        <- Sigma.hat
  result$kls              <- kls
  
  return(result)  
}


####################################################################
#### visualization
####################################################################

vis_arrange <- function(vdat1, vdat2, combined.legend, color.pal = brewer.pal(6, "Set1"), shape.pal = c(8, 13, 15, 16, 17, 18), alpha.value = 0.7, size.legend = 28, size.lab = 28, size.text = 18){
  
  xlabel1 <- sort(unique(vdat1$m))
  plot1   <- ggplot(vdat1, aes(x=m, y = log10(KL), col = method)) + 
    geom_point(aes(shape = method), size = 3) + 
    geom_line(size = 1, alpha = alpha.value) +
    ylab('log10(KL)') + 
    scale_x_continuous(name = 'm', limits=range(xlabel1), breaks=xlabel1) +
    scale_color_manual(values = color.pal, labels = combined.legend) +
    scale_shape_manual(values = shape.pal, labels = combined.legend) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text),
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text),
          legend.title = element_blank(),
          legend.text=element_text(size = size.legend),
          legend.direction = 'horizontal',
          plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt")) # t, r, b, l
  
  xlabel2 <- sort(unique(vdat2$period))
  plot2   <- ggplot(vdat2, aes(x=period, y = log10(KL), col = method)) + 
    geom_point(aes(shape = method), size = 2) + 
    geom_line(size = 1, alpha = alpha.value) +
    ylab('log10(KL)') + 
    scale_x_continuous(name = 'period', limits=range(xlabel2), breaks=xlabel2) +
    scale_color_manual(values = color.pal, labels = combined.legend) +
    scale_shape_manual(values = shape.pal, labels = combined.legend) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text),
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text),
          plot.margin = unit(c(5.5, 5.5, 5.5, 20), "pt"))
  
  tmp <- ggplot_gtable(ggplot_build(plot1))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  mylegend <- tmp$grobs[[leg]]
  
  grid.arrange(mylegend, 
               arrangeGrob(plot1 + theme(legend.position="none"), 
                           plot2 + theme(legend.position="none"),
                           nrow=1), 
               nrow=2,heights=c(1, 10))
}


####################################################################
#### simulation 1: Dampedcosine + '1-rho' distance 
####################################################################

cand.m            <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45) ; n.cand.m <- length(cand.m)
cand.period       <- c(0.05, 0.07, 0.09, 0.1) ; n.cand.period <- length(cand.period)
sim1              <- list()

# # small case
# cand.m            <- c(10, 20, 30) ; n.cand.m <- length(cand.m)
# cand.period        <- c(0.05, 0.1) ; n.cand.period <- length(cand.period)
# sim1              <- list()

cand.all            <- expand.grid(cand.m, cand.period)
cand.all            <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all)  <- c('index', 'm', 'period')
n.cand.all          <- n.cand.m * n.cand.period

no_cores            <- parallel::detectCores() - 2
cl                  <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)
sim1 <- foreach(m = cand.all$m, p = cand.all$period, .export = c("order_coordinate", "covmodel_dampedsine", "covmodel_dampedcosine", "covmodel_besselJ", "covmodel_wave", "positive_def", "conditioning_nn", "correlation", "corrvecchia_knownCovparms", "distance_correlation", "kldiv", "order_maxmin_correlation", "order_maxmin_correlation_old", "order_maxmin_euclidean", "simulation", "vecchia_specify_adjusted"), .packages='GPvecchia') %dopar% simulation(n = 30^2, m = m, covparms = c(1, 1, p), covm = "Dampedcosine", def.dist = NULL, tol = 1e-4)
parallel::stopCluster(cl)

kls.maxmin.eucord.euccond.ref   <- rep(NA, n.cand.all)
kls.xcoord.eucord.euccond       <- rep(NA, n.cand.all)
kls.ycoord.eucord.euccond       <- rep(NA, n.cand.all)
kls.maxmin.eucord.euccond       <- rep(NA, n.cand.all)
kls.maxmin.eucord.corcond       <- rep(NA, n.cand.all)
kls.maxmin.corord.euccond       <- rep(NA, n.cand.all)
kls.maxmin.corord.corcond       <- rep(NA, n.cand.all)
kls.xcoord.eucord.corcond       <- rep(NA, n.cand.all)
kls.ycoord.eucord.corcond       <- rep(NA, n.cand.all)
kls.xycoord.eucord.corcond      <- rep(NA, n.cand.all)
kls.xycoord.eucord.euccond      <- rep(NA, n.cand.all)
kls.xycoord.eucord.euccond.ref  <- rep(NA, n.cand.all)
for(i in 1:n.cand.all) {
  kls.maxmin.eucord.euccond.ref[i]  <- sim1[[i]]$kls[1]
  kls.xcoord.eucord.euccond[i]      <- sim1[[i]]$kls[2]
  kls.ycoord.eucord.euccond[i]      <- sim1[[i]]$kls[3]
  kls.maxmin.eucord.euccond[i]      <- sim1[[i]]$kls[4]
  kls.maxmin.eucord.corcond[i]      <- sim1[[i]]$kls[5]
  kls.maxmin.corord.euccond[i]      <- sim1[[i]]$kls[6]
  kls.maxmin.corord.corcond[i]      <- sim1[[i]]$kls[7]
  kls.xcoord.eucord.corcond[i]      <- sim1[[i]]$kls[8]
  kls.ycoord.eucord.corcond[i]      <- sim1[[i]]$kls[9]
  kls.xycoord.eucord.corcond[i]     <- sim1[[i]]$kls[10]
  kls.xycoord.eucord.euccond[i]     <- sim1[[i]]$kls[11]
  kls.xycoord.eucord.euccond.ref[i] <- sim1[[i]]$kls[12]
}

sqrt(sum((kls.maxmin.eucord.euccond.ref - kls.maxmin.eucord.euccond)^2))
sqrt(sum((kls.xycoord.eucord.euccond.ref - kls.xycoord.eucord.euccond)^2))

err.modifying1 <- c()
for(i in 1:length(sim1)) err.modifying1[i] <- sqrt(sum((sim1[[i]]$Sigma - sim1[[i]]$Sigma.modified))^2)
max(err.modifying1)

set.period  <- 0.05
ind         <- cand.all$period == set.period
vis.dat1    <- data.frame(kls.xcoord.eucord.euccond[ind], kls.ycoord.eucord.euccond[ind], kls.maxmin.eucord.euccond[ind], kls.maxmin.eucord.corcond[ind], kls.maxmin.corord.euccond[ind], kls.maxmin.corord.corcond[ind])
vis.dat1    <- vis.dat1[, order(colnames(vis.dat1))]
vis.dat1    <- cbind(rep(cand.m, times = ncol(vis.dat1)), tidyr::gather(vis.dat1))
colnames(vis.dat1) <- c("m", "method", "KL")
head(vis.dat1)

set.m       <- 30
ind         <- cand.all$m == set.m
vis.dat2    <- data.frame(kls.xcoord.eucord.euccond[ind], kls.ycoord.eucord.euccond[ind], kls.maxmin.eucord.euccond[ind], kls.maxmin.eucord.corcond[ind], kls.maxmin.corord.euccond[ind], kls.maxmin.corord.corcond[ind])
vis.dat2    <- vis.dat2[, order(colnames(vis.dat2))]
vis.dat2    <- cbind(rep(cand.period, times = ncol(vis.dat2)), tidyr::gather(vis.dat2))
colnames(vis.dat2) <- c("period", "method", "KL")
head(vis.dat2)

kls.legend <- c("C-Maxmin + C-NN", "C-Maxmin + E-NN", "E-Maxmin + C-NN", "E-Maxmin + E-NN", "X-Coord + E-NN", "Y-Coord + E-NN")
vis_arrange(vdat1 = vis.dat1, vdat2 = vis.dat2, combined.legend = kls.legend, color.pal = brewer.pal(6, "Set1"), shape.pal = c(16, 17, 15, 18, 8, 13), alpha.value = 0.7, size.legend = 14, size.lab = 14, size.text = 12)

# save(sim1, cand.all, set.period, set.m, vis.dat1, vis.dat2, kls.legend, err.modifying1, vis_arrange, file='2_corrvecchia/sim_periodicity_1D_1.RData')
# rm(sim1, cand.all, vis.dat1, vis.dat2, kls.legend, err.modifying1)
# load(file='2_corrvecchia/sim_periodicity_1D_1.RData')


####################################################################
#### simulation 2: Dampedcosine + '1-|rho|' distance 
####################################################################

cand.m            <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45) ; n.cand.m <- length(cand.m)
cand.period       <- c(0.05, 0.07, 0.09, 0.1) ; n.cand.period <- length(cand.period)
sim2              <- list()

# # small case
# cand.m            <- c(10, 20, 30) ; n.cand.m <- length(cand.m)
# cand.period        <- c(0.05, 0.1) ; n.cand.period <- length(cand.period)
# sim2              <- list()

cand.all            <- expand.grid(cand.m, cand.period)
cand.all            <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all)  <- c('index', 'm', 'period')
n.cand.all          <- n.cand.m * n.cand.period

no_cores            <- parallel::detectCores() - 2
cl                  <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)
sim2 <- foreach(m = cand.all$m, p = cand.all$period, .export = c("order_coordinate", "covmodel_dampedsine", "covmodel_dampedcosine", "covmodel_besselJ", "covmodel_wave", "positive_def", "conditioning_nn", "correlation", "corrvecchia_knownCovparms", "distance_correlation", "kldiv", "order_maxmin_correlation", "order_maxmin_correlation_old", "order_maxmin_euclidean", "simulation", "vecchia_specify_adjusted"), .packages='GPvecchia') %dopar% simulation(n = 30^2, m = m, covparms = c(1, 1, p), covm = "Dampedcosine", def.dist = "abs", tol = 1e-4)
parallel::stopCluster(cl)

kls.maxmin.eucord.euccond.ref   <- rep(NA, n.cand.all)
kls.xcoord.eucord.euccond       <- rep(NA, n.cand.all)
kls.ycoord.eucord.euccond       <- rep(NA, n.cand.all)
kls.maxmin.eucord.euccond       <- rep(NA, n.cand.all)
kls.maxmin.eucord.corcond       <- rep(NA, n.cand.all)
kls.maxmin.corord.euccond       <- rep(NA, n.cand.all)
kls.maxmin.corord.corcond       <- rep(NA, n.cand.all)
kls.xcoord.eucord.corcond       <- rep(NA, n.cand.all)
kls.ycoord.eucord.corcond       <- rep(NA, n.cand.all)
kls.xycoord.eucord.corcond      <- rep(NA, n.cand.all)
kls.xycoord.eucord.euccond      <- rep(NA, n.cand.all)
kls.xycoord.eucord.euccond.ref  <- rep(NA, n.cand.all)
for(i in 1:n.cand.all) {
  kls.maxmin.eucord.euccond.ref[i]  <- sim2[[i]]$kls[1]
  kls.xcoord.eucord.euccond[i]      <- sim2[[i]]$kls[2]
  kls.ycoord.eucord.euccond[i]      <- sim2[[i]]$kls[3]
  kls.maxmin.eucord.euccond[i]      <- sim2[[i]]$kls[4]
  kls.maxmin.eucord.corcond[i]      <- sim2[[i]]$kls[5]
  kls.maxmin.corord.euccond[i]      <- sim2[[i]]$kls[6]
  kls.maxmin.corord.corcond[i]      <- sim2[[i]]$kls[7]
  kls.xcoord.eucord.corcond[i]      <- sim2[[i]]$kls[8]
  kls.ycoord.eucord.corcond[i]      <- sim2[[i]]$kls[9]
  kls.xycoord.eucord.corcond[i]     <- sim2[[i]]$kls[10]
  kls.xycoord.eucord.euccond[i]     <- sim2[[i]]$kls[11]
  kls.xycoord.eucord.euccond.ref[i] <- sim2[[i]]$kls[12]
}

sqrt(sum((kls.maxmin.eucord.euccond.ref - kls.maxmin.eucord.euccond)^2))
sqrt(sum((kls.xycoord.eucord.euccond.ref - kls.xycoord.eucord.euccond)^2))

err.modifying2 <- c()
for(i in 1:length(sim2)) err.modifying2[i] <- sqrt(sum((sim2[[i]]$Sigma - sim2[[i]]$Sigma.modified))^2)
max(err.modifying2)

set.period  <- 0.05
ind         <- cand.all$period == set.period
vis.dat1    <- data.frame(kls.xcoord.eucord.euccond[ind], kls.ycoord.eucord.euccond[ind], kls.maxmin.eucord.euccond[ind], kls.maxmin.eucord.corcond[ind], kls.maxmin.corord.euccond[ind], kls.maxmin.corord.corcond[ind])
vis.dat1    <- vis.dat1[, order(colnames(vis.dat1))]
vis.dat1    <- cbind(rep(cand.m, times = ncol(vis.dat1)), tidyr::gather(vis.dat1))
colnames(vis.dat1) <- c("m", "method", "KL")
head(vis.dat1)

set.m       <- 30
ind         <- cand.all$m == set.m
vis.dat2    <- data.frame(kls.xcoord.eucord.euccond[ind], kls.ycoord.eucord.euccond[ind], kls.maxmin.eucord.euccond[ind], kls.maxmin.eucord.corcond[ind], kls.maxmin.corord.euccond[ind], kls.maxmin.corord.corcond[ind])
vis.dat2    <- vis.dat2[, order(colnames(vis.dat2))]
vis.dat2    <- cbind(rep(cand.period, times = ncol(vis.dat2)), tidyr::gather(vis.dat2))
colnames(vis.dat2) <- c("period", "method", "KL")
head(vis.dat2)

kls.legend <- c("C-Maxmin + C-NN", "C-Maxmin + E-NN", "E-Maxmin + C-NN", "E-Maxmin + E-NN", "X-Coord + E-NN", "Y-Coord + E-NN")
vis_arrange(vdat1 = vis.dat1, vdat2 = vis.dat2, combined.legend = kls.legend, color.pal = brewer.pal(6, "Set1"), shape.pal = c(16, 17, 15, 18, 8, 13), alpha.value = 0.7, size.legend = 14, size.lab = 14, size.text = 12)

# save(sim2, cand.all, set.period, set.m, vis.dat1, vis.dat2, kls.legend, err.modifying2, vis_arrange, file='2_corrvecchia/sim_periodicity_1D_2.RData')
# rm(sim2, cand.all, vis.dat1, vis.dat2, kls.legend, err.modifying2)
# load(file='2_corrvecchia/sim_periodicity_1D_2.RData')