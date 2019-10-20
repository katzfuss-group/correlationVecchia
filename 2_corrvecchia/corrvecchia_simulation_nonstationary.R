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
#### Nonstaionary Matern covariance model (Risser MD, Calder CA (2015))
####################################################################

covparms <- c(1)

a <- function(loc) 0.47 * loc[1] + 0.03
b <- function(loc) 1
angle <- function(loc) 0

# spatially-varying standard deviation
sigma <- function(loc) determinant(aniso_mat(loc), logarithm = F)[[1]][1]^0.25
# spatially-varying local anisotropy (controlling both the range and direction of dependence)
aniso_mat<- function(loc) {
  
  eta <- angle(loc)
  rot.mat <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
  
  range <- c(a(loc)^(-2), b(loc)^(-2)) * 0.01
  diag.mat <- diag(range, nrow = length(loc))
  
  aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
  
  return(aniso.mat)
}
# matern's smoothness
smoothness <- function(loc) 0.2 * exp(loc[1])

# Katzfuss M. 2013. Bayesian nonstationary spatial modeling for very large datasets. Environmetrics 24(3):189–200.
# Stein ML. 2005. Nonstationary spatial covariance functions, Technical Report, University of Chicago, Department of Statistics.

matern_ns <- function(locs1, locs2 = NULL) {
  
  if(is.null(locs2)) locs2 = locs1
  
  n1 <- nrow(locs1) ; n2 <- nrow(locs2) ; d <- ncol(locs1)
  
  mat.cov <- matrix(NA, nrow = n1, ncol = n2)
  
  for(i in 1:n1){
    for(j in 1:n2){
      sigma.ij    <- sigma(locs1[i, ]) * sigma(locs2[j, ])
      kernel.ij   <- ( aniso_mat(locs1[i, ]) + aniso_mat(locs2[j, ]) ) / 2 
      smooth.ij   <- ( smoothness(locs1[i, ]) + smoothness(locs2[j, ]) ) / 2
      
      q.ij <- as.numeric(crossprod( locs1[i, ] - locs2[j, ], solve(kernel.ij, locs1[i, ] - locs2[j, ]) ))
      
      mat.cov[i,j] <- sigma.ij * fields::Matern( sqrt(q.ij), nu = smooth.ij) / sqrt( determinant(kernel.ij, logarithm = F)[[1]][1] )
    }
  }
  
  return(mat.cov)
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
  
  xlabel2 <- sort(unique(vdat2$m))
  plot2   <- ggplot(vdat2, aes(x=m, y = log10(KL), col = method)) + 
    geom_point(aes(shape = method), size = 2) + 
    geom_line(size = 1, alpha = alpha.value) +
    ylab('log10(KL)') + 
    scale_x_continuous(name = 'm', limits=range(xlabel2), breaks=xlabel2) +
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
#### simulation function
####################################################################

positive_def <- function(Sigma, tol){
  eig.decomp  <- eigen(Sigma)
  diagvec     <- ifelse(eig.decomp$values < tol, tol, eig.decomp$values)
  
  Sigma.modified <- eig.decomp$vectors %*% diag(diagvec) %*% t(eig.decomp$vectors)
  return(Sigma.modified)
}

simulation <- function(n = 15^2, m = 10, covparms = c(1), tol = 1e-6) {

  locs      <- matrix(runif(n * 2, 0, 1), n, 2)
  Sigma     <- matern_ns(locs)
  
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
  # euclidean-based ordering + euclidean-based conditioning
  approx[[4]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma.modified, covparms = covparms)
  # euclidean-based ordering + correlation-based conditioning
  approx[[5]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = covparms)
  # correlation-based ordering + euclidean-based conditioning
  approx[[6]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma.modified, covparms = covparms)
  # correlation-based ordering + correlation-based conditioning
  approx[[7]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma.modified, covparms = covparms)
  
  ### compute approximate covariance matrices
  n.approx    <- length(approx)
  Sigma.hat   <- list()
  kls         <- rep(NA, n.approx)
  for(i in 1:n.approx){
    
    Sigma.ord       <- matern_ns(locs1 = approx[[i]]$locsord, locs2 = NULL) # true cov in appropriate ordering
    
    Sigma.ord.modified <- positive_def(Sigma.ord, tol)
    
    U               <- createU(approx[[i]], c(1, 1, 1), 0, covmodel = Sigma.ord.modified)$U
    revord          <- order(approx[[i]]$ord)
    Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    kls[i]          <- kldiv(Sigma.modified, Sigma.hat[[i]])
  }
  
  result                  <- list()
  result$n                <- n
  result$m                <- m
  result$covparms         <- covparms
  result$locs             <- locs
  result$Sigma            <- Sigma
  result$Sigma.modified   <- Sigma.modified
  result$approx           <- approx
  result$Sigma.hat        <- Sigma.hat
  result$kls              <- kls

  return(result)  
}


####################################################################
#### simulation 1: smoothness
####################################################################

cand.m    <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45) ; n.cand.m <- length(cand.m)
sim1      <- list()

# # small case
# cand.m    <- c(1, 10, 20, 30, 40) ; n.cand.m <- length(cand.m)
# sim1      <- list()

a <- function(loc) 1
b <- function(loc) 1
angle <- function(loc) 0
smoothness <- function(loc) 0.2 + 1.3 * loc[1]

# plot(seq(0, 5, 0.01), fields::Matern(seq(0, 5, 0.01), smoothness = 1), type = "l", col = 1)
# lines(seq(0, 5, 0.01), fields::Matern(seq(0, 5, 0.01), smoothness = 0.2), type = "l", col = 2)
# lines(seq(0, 5, 0.01), fields::Matern(seq(0, 5, 0.01), smoothness = 1.5), type = "l", col = 3)
# lines(seq(0, 5, 0.01), fields::Matern(seq(0, 5, 0.01), smoothness = 2), type = "l", col = 4)

no_cores            <- parallel::detectCores() - 2
cl                  <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)
sim1 <- foreach(m = cand.m, .export = c("a", "b", "aniso_mat", "conditioning_nn", "correlation", "corrvecchia_knownCovparms", "distance_correlation", "kldiv", "matern_ns", "order_maxmin_correlation", "order_maxmin_correlation_old", "order_maxmin_euclidean", "simulation", "smoothness", "vecchia_specify_adjusted"), .packages='GPvecchia') %dopar% simulation(30^2, m = m, covparms = c(1), tol = 1e-4)
parallel::stopCluster(cl)

# for(i in 1:n.cand.m) sim[[i]] <- simulation(30^2, m = cand.m[i], covparms = c(1))

kls.maxmin.eucord.euccond.ref <- rep(NA, n.cand.m)
kls.xcoord.eucord.euccond     <- rep(NA, n.cand.m)
kls.ycoord.eucord.euccond     <- rep(NA, n.cand.m)
kls.maxmin.eucord.euccond     <- rep(NA, n.cand.m)
kls.maxmin.eucord.corcond     <- rep(NA, n.cand.m)
kls.maxmin.corord.euccond     <- rep(NA, n.cand.m)
kls.maxmin.corord.corcond     <- rep(NA, n.cand.m)
for(i in 1:n.cand.m) {
  kls.maxmin.eucord.euccond.ref[i]  <- sim1[[i]]$kls[1]
  kls.xcoord.eucord.euccond[i]      <- sim1[[i]]$kls[2]
  kls.ycoord.eucord.euccond[i]      <- sim1[[i]]$kls[3]
  kls.maxmin.eucord.euccond[i]      <- sim1[[i]]$kls[4]
  kls.maxmin.eucord.corcond[i]      <- sim1[[i]]$kls[5]
  kls.maxmin.corord.euccond[i]      <- sim1[[i]]$kls[6]
  kls.maxmin.corord.corcond[i]      <- sim1[[i]]$kls[7]
}

sqrt(sum((kls.maxmin.eucord.euccond.ref - kls.maxmin.eucord.euccond)^2))

vis.dat1 <- data.frame(kls.xcoord.eucord.euccond, kls.ycoord.eucord.euccond, kls.maxmin.eucord.euccond, kls.maxmin.eucord.corcond, kls.maxmin.corord.euccond, kls.maxmin.corord.corcond)
vis.dat1 <- vis.dat1[, order(colnames(vis.dat1))]
head(vis.dat1)

err.modifying1 <- c()
for(i in 1:length(sim1)) err.modifying1[i] <- sqrt(sum((sim1[[i]]$Sigma - sim1[[i]]$Sigma.modified))^2)
max(err.modifying1)

plot(cand.m, log10(vis.dat1$kls.maxmin.eucord.euccond), type = "o", col = 1, lty = 1, lwd = 3, ylim = c(min(log10(vis.dat1)), max(log10(vis.dat1))), xlab = "m", ylab = "log10(KL)", main = NULL)
lines(cand.m, log10(vis.dat1$kls.maxmin.corord.corcond), type = "o", col = 2, lty = 2, lwd = 3)
lines(cand.m, log10(vis.dat1$kls.maxmin.eucord.corcond), type = "o", col = 3, lty = 3, lwd = 3)
lines(cand.m, log10(vis.dat1$kls.maxmin.corord.euccond), type = "o", col = 4, lty = 4, lwd = 3)
lines(cand.m, log10(vis.dat1$kls.xcoord.eucord.euccond), type = "o", col = 5, lty = 5, lwd = 3)
lines(cand.m, log10(vis.dat1$kls.ycoord.eucord.euccond), type = "o", col = 6, lty = 6, lwd = 3)
legend("topright", legend=c("Maxmin + E.ord + E.cond", "Maxmin + C.ord + C.cond", "Maxmin + E.ord + C.cond", "Maxmin + C.ord + E.cond", "X-coord + E.ord + E.cond", "Y-coord + E.ord + E.cond"), col=1:6, lty=1:6, lwd = 3, cex=1)

# save(sim1, cand.m, vis.dat1, err.modifying1, file='2_corrvecchia/sim_nonstationarity_1.RData')
# rm(sim1, cand.m)
# load(file='2_corrvecchia/sim_nonstationarity_1.RData')


####################################################################
#### simulation 2: range
####################################################################

cand.m    <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45) ; n.cand.m <- length(cand.m)
sim2      <- list()

# # small case
# cand.m    <- c(1, 10, 20, 30, 40) ; n.cand.m <- length(cand.m)
# sim2      <- list()

a <- function(loc) 1/(0.01 + 0.99 * loc[1])/10
b <- function(loc) 1/(0.01 + 0.99 * loc[1])/10
angle <- function(loc) 0
smoothness <- function(loc) 0.5

no_cores            <- parallel::detectCores() - 2
cl                  <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)
sim2 <- foreach(m = cand.m, .export = c("a", "b", "aniso_mat", "conditioning_nn", "correlation", "corrvecchia_knownCovparms", "distance_correlation", "kldiv", "matern_ns", "order_maxmin_correlation", "order_maxmin_correlation_old", "order_maxmin_euclidean", "simulation", "smoothness", "vecchia_specify_adjusted"), .packages='GPvecchia') %dopar% simulation(30^2, m = m, covparms = c(1), tol = 1e-4)
parallel::stopCluster(cl)

kls.maxmin.eucord.euccond.ref <- rep(NA, n.cand.m)
kls.xcoord.eucord.euccond     <- rep(NA, n.cand.m)
kls.ycoord.eucord.euccond     <- rep(NA, n.cand.m)
kls.maxmin.eucord.euccond     <- rep(NA, n.cand.m)
kls.maxmin.eucord.corcond     <- rep(NA, n.cand.m)
kls.maxmin.corord.euccond     <- rep(NA, n.cand.m)
kls.maxmin.corord.corcond     <- rep(NA, n.cand.m)
for(i in 1:n.cand.m) {
  kls.maxmin.eucord.euccond.ref[i]  <- sim2[[i]]$kls[1]
  kls.xcoord.eucord.euccond[i]      <- sim2[[i]]$kls[2]
  kls.ycoord.eucord.euccond[i]      <- sim2[[i]]$kls[3]
  kls.maxmin.eucord.euccond[i]      <- sim2[[i]]$kls[4]
  kls.maxmin.eucord.corcond[i]      <- sim2[[i]]$kls[5]
  kls.maxmin.corord.euccond[i]      <- sim2[[i]]$kls[6]
  kls.maxmin.corord.corcond[i]      <- sim2[[i]]$kls[7]
}

sqrt(sum((kls.maxmin.eucord.euccond.ref - kls.maxmin.eucord.euccond)^2))

vis.dat2 <- data.frame(kls.xcoord.eucord.euccond, kls.ycoord.eucord.euccond, kls.maxmin.eucord.euccond, kls.maxmin.eucord.corcond, kls.maxmin.corord.euccond, kls.maxmin.corord.corcond)
vis.dat2 <- vis.dat2[, order(colnames(vis.dat2))]
head(vis.dat2)

err.modifying2 <- c()
for(i in 1:length(sim2)) err.modifying2[i] <- sqrt(sum((sim2[[i]]$Sigma - sim2[[i]]$Sigma.modified))^2)
max(err.modifying2)

plot(cand.m, log10(vis.dat2$kls.maxmin.eucord.euccond), type = "o", col = 1, lty = 1, lwd = 3, ylim = c(min(log10(vis.dat2)), max(log10(vis.dat2))), xlab = "m", ylab = "log10(KL)", main = NULL)
lines(cand.m, log10(vis.dat2$kls.maxmin.corord.corcond), type = "o", col = 2, lty = 2, lwd = 3)
lines(cand.m, log10(vis.dat2$kls.maxmin.eucord.corcond), type = "o", col = 3, lty = 3, lwd = 3)
lines(cand.m, log10(vis.dat2$kls.maxmin.corord.euccond), type = "o", col = 4, lty = 4, lwd = 3)
lines(cand.m, log10(vis.dat2$kls.xcoord.eucord.euccond), type = "o", col = 5, lty = 5, lwd = 3)
lines(cand.m, log10(vis.dat2$kls.ycoord.eucord.euccond), type = "o", col = 6, lty = 6, lwd = 3)
legend("topright", legend=c("Maxmin + E.ord + E.cond", "Maxmin + C.ord + C.cond", "Maxmin + E.ord + C.cond", "Maxmin + C.ord + E.cond", "X-coord + E.ord + E.cond", "Y-coord + E.ord + E.cond"), col=1:6, lty=1:6, lwd = 3, cex=1)

# save(sim2, cand.m, vis.dat2, err.modifying2, file='2_corrvecchia/sim_nonstationarity_2.RData')
# rm(sim2, cand.m, kls.xcoord.eucord.euccond, kls.ycoord.eucord.euccond, kls.maxmin.eucord.euccond, kls.maxmin.eucord.corcond, kls.maxmin.corord.euccond, kls.maxmin.corord.corcond)
# load(file='2_corrvecchia/sim_nonstationarity_2.RData')


####################################################################
#### visualization step
####################################################################

vis.dat1    <- vis.dat1[, order(colnames(vis.dat1))]
vis.dat1    <- cbind(rep(cand.m, times = ncol(vis.dat1)), tidyr::gather(vis.dat1))
colnames(vis.dat1) <- c("m", "method", "KL")
head(vis.dat1)

vis.dat2    <- vis.dat2[, order(colnames(vis.dat2))]
vis.dat2    <- cbind(rep(cand.m, times = ncol(vis.dat2)), tidyr::gather(vis.dat2))
colnames(vis.dat2) <- c("m", "method", "KL")
head(vis.dat2)

kls.legend <- c("Maxmin + C.ord + C.cond", "Maxmin + C.ord + E.cond", "Maxmin + E.ord + C.cond", "Maxmin + E.ord + E.cond", "X-coord + E.ord + E.cond", "Y-coord + E.ord + E.cond")
vis_arrange(vdat1 = vis.dat1, vdat2 = vis.dat2, combined.legend = kls.legend, color.pal = brewer.pal(6, "Set1"), shape.pal = c(16, 17, 15, 18, 8, 13), alpha.value = 0.7, size.legend = 14, size.lab = 14, size.text = 12)

