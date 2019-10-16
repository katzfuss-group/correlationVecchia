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

# covparms = c(omega, nu), omega > 0, nu > (d-2)*d/2
covparms <- c(1, 1)

covmodel_periodic <- function(locs, covparms) {
  h <- fields::rdist(locs)
  c <- gamma(covparms[2] + 1) * (2/covparms[1]/c(h))^covparms[2] * Bessel::BesselJ(covparms[1] * c(h), covparms[2])
  c[which(is.nan(c))] <- 1
  matrix(c, nrow = nrow(h), ncol = ncol(h))
}

covparms <- c(2 , 0)

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')

covparms <- c(4 , 0)

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')

covparms <- c(6 , 0)

n <- 15^2
locs <- matrix(runif(n * 2, 0, 1), n, 2)
h <- fields::rdist(locs) ; covh <- covmodel_periodic(locs, covparms)
min(eigen(covh)$value)

ind <- order(unlist(h)) ; hord <- unlist(h)[ind]
covhord <- unlist(covh)[ind]
plot(hord, covhord, type = 'l')
abline(h = 0, col = 'red')

covparms <- c(8 , 0)

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

positive_def <- function(Sigma, tol){
  eig.decomp  <- eigen(Sigma)
  diagvec     <- ifelse(eig.decomp$values < tol, tol, eig.decomp$values)
  
  Sigma.modified <- eig.decomp$vectors %*% diag(diagvec) %*% t(eig.decomp$vectors)
  return(Sigma.modified)
}

simulation <- function(n = 15^2, m = 10, covparms = c(1, 1/12), tol = 1e-6) {
  
  locs      <- matrix(runif(n * 2, 0, 1), n, 2)
  Sigma     <- covmodel_periodic(locs, covparms)
  
  Sigma.modified <- positive_def(Sigma, tol)
  
  y         <- as.numeric(t(chol(Sigma.modified)) %*% rnorm(n))
  
  ### specify vecchia approximations
  approx <- list()
  
  # standard vecchia with maxmin ordering
  approx[[1]]           <- vecchia_specify_adjusted(locs, m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
  # standard vecchia with x coord ordering
  approx[[2]]           <- vecchia_specify_adjusted(locs, m, ordering = "coord", which.coord = 1, cond.yz='y', conditioning = "NN")
  # standard vecchia with y coord ordering
  approx[[3]]           <- vecchia_specify_adjusted(locs, m, ordering = "coord", which.coord = 2, cond.yz='y', conditioning = "NN")
  # correlation-based vecchia with the corrvecchia function
  approx[[4]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", covmodel = Sigma.modified, covparms = covparms)
  
  ### compute approximate covariance matrices
  Sigma.hat   <- list()
  kls         <- c()
  for(i in 1:4){
    
    Sigma.ord       <- covmodel_periodic(approx[[i]]$locsord, covparms) # true cov in appropriate ordering
    
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
  result$locs             <- locs
  result$approx           <- approx
  result$kls              <- kls
  result$Sigma            <- Sigma
  result$Sigma.modified   <- Sigma.modified
  result$Sigma.hat        <- Sigma.hat
  
  return(result)  
}


####################################################################
#### visualization
####################################################################

vis_arrange <- function(vdat1, vdat2, combined.legend, color.pal = brewer.pal(4, "Set1"), alpha.value = 0.7, size.legend = 28, size.lab = 28, size.text = 18){
  
  xlabel1 <- sort(unique(vdat1$m))
  plot1   <- ggplot(vdat1, aes(x=m, y = log10(KL), col = method)) + 
    geom_point(aes(shape = method), size = 3) + 
    geom_line(size = 1, alpha = alpha.value) +
    ylab('log10(KL)') + 
    scale_x_discrete(name = 'm', limits=xlabel1, labels=as.character(xlabel1)) +
    scale_color_manual(values = color.pal, labels = combined.legend) +
    scale_shape_manual(values = c(15, 16, 17, 18), labels = combined.legend) +
    theme(axis.title.x = element_text(size = size.lab), 
          axis.text.x = element_text(size = size.text),
          axis.title.y = element_text(size = size.lab), 
          axis.text.y = element_text(size = size.text),
          legend.title = element_blank(),
          legend.text=element_text(size = size.legend),
          legend.direction = 'horizontal',
          plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt")) # t, r, b, l
  
  xlabel2 <- sort(unique(vdat2$omega))
  plot2   <- ggplot(vdat2, aes(x=omega, y = log10(KL), col = method)) +
    geom_point(aes(shape = method), size = 2) + 
    geom_line(size = 1, alpha = alpha.value) +
    ylab('log10(KL)') + 
    scale_x_continuous(name = 'omega', limits=range(xlabel2), breaks=xlabel2) +
    scale_color_manual(values = color.pal, labels = combined.legend) +
    scale_shape_manual(values = c(15, 16, 17, 18), labels = combined.legend) +
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
#### simulation 1: 
####################################################################

cand.m            <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45) ; n.cand.m <- length(cand.m)
cand.omega        <- c(2, 4, 6, 8) ; n.cand.omega <- length(cand.omega)
sim1              <- list()

cand.all            <- expand.grid(cand.m, cand.omega )
cand.all            <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all)  <- c('index', 'm', 'omega ')
n.cand.all          <- n.cand.m * n.cand.omega 

no_cores            <- parallel::detectCores() - 2
cl                  <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)
sim1 <- foreach(m = cand.all$m, o = cand.all$omega, .export = c("positive_def", "conditioning_nn", "correlation", "corrvecchia_knownCovparms", "distance_correlation", "kldiv", "order_maxmin_correlation", "order_maxmin_correlation_old", "order_maxmin_euclidean", "simulation", "vecchia_specify_adjusted"), .packages='GPvecchia') %dopar% simulation(30^2, m = m, covparms = c(o, 0), tol = 1e-3)
parallel::stopCluster(cl)

kls.maxmin.euclidean    <- rep(NA, n.cand.all)
kls.maxmin.corr         <- rep(NA, n.cand.all)
kls.xcoord.euclidean    <- rep(NA, n.cand.all)
kls.ycoord.euclidean    <- rep(NA, n.cand.all)
for(i in 1:n.cand.all) {
  kls.maxmin.euclidean[i]    <- sim1[[i]]$kls[1]
  kls.maxmin.corr[i]         <- sim1[[i]]$kls[4]
  kls.xcoord.euclidean[i]    <- sim1[[i]]$kls[2]
  kls.ycoord.euclidean[i]    <- sim1[[i]]$kls[3]
}

set.omega   <- 6
ind         <- cand.all$omega == set.omega
vis.dat1    <- data.frame(kls.maxmin.euclidean[ind], kls.maxmin.corr[ind], kls.xcoord.euclidean[ind], kls.ycoord.euclidean[ind])
vis.dat1    <- vis.dat1[, order(colnames(vis.dat1))]
vis.dat1    <- cbind(rep(cand.m, times = 4), tidyr::gather(vis.dat1))
colnames(vis.dat1) <- c("m", "method", "KL")
head(vis.dat1)

set.m       <- 30
ind         <- cand.all$m == set.m
vis.dat2    <- data.frame(kls.maxmin.euclidean[ind], kls.maxmin.corr[ind], kls.xcoord.euclidean[ind], kls.ycoord.euclidean[ind])
vis.dat2    <- vis.dat2[, order(colnames(vis.dat2))]
vis.dat2    <- cbind(rep(cand.omega, times = 4), tidyr::gather(vis.dat2))
colnames(vis.dat2) <- c("omega", "method", "KL")
head(vis.dat2)

kls.legend <- c("Correlation + Maxmin     ", "Euclidean + Maxmin     ", "Euclidean + x-coord     ", "Euclidean + y-coord")
vis_arrange(vdat1 = vis.dat1, vdat2 = vis.dat2, combined.legend = kls.legend, color.pal = brewer.pal(4, "Set1"), alpha.value = 0.7, size.legend = 16, size.lab = 16, size.text = 12)

err.modifying <- c()
for(i in 1:length(sim1)) err.modifying[i] <- sqrt(sum((sim1[[i]]$Sigma - sim1[[i]]$Sigma.modified))^2)
max(err.modifying)

# save(sim1, cand.all, vis.dat1, vis.dat2, kls.legend, err.modifying, file='2_corrvecchia/sim_bessel_with_ab_corr_dist_1.RData')
# load(file='2_corrvecchia/sim_bessel_with_ab_corr_dist_1.RData')

