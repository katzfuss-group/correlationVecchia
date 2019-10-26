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

# spatially-varying standard deviation
sigma <- function(loc, deg.aniso) determinant(aniso_mat(loc, deg.aniso), logarithm = F)[[1]][1]^0.25
# spatially-varying local anisotropy (controlling both the range and direction of dependence)
aniso_mat <- function(loc, deg.aniso) {
  
  eta       <- 0
  rot.mat   <- matrix(c(cos(eta), sin(eta), -sin(eta), cos(eta)), nrow = length(loc), ncol = length(loc), byrow = T)
  
  a         <- function(loc) deg.aniso
  
  range     <- c(a(loc)^(-2), 1) * 0.01
  diag.mat  <- diag(range, nrow = length(loc))
  
  aniso.mat <- t(rot.mat) %*% diag.mat %*% rot.mat
  
  return(aniso.mat)
}
# matern's smoothness
smoothness <- function(loc) 0.5

# Katzfuss M. 2013. Bayesian nonstationary spatial modeling for very large datasets. Environmetrics 24(3):189–200.
# Stein ML. 2005. Nonstationary spatial covariance functions, Technical Report, University of Chicago, Department of Statistics.

matern_ns <- function(locs1, locs2 = NULL, deg.aniso) {
  
  if(is.null(locs2)) locs2 <- locs1
  
  n1 <- nrow(locs1) ; n2 <- nrow(locs2) ; d <- ncol(locs1)
  
  mat.cov <- matrix(NA, nrow = n1, ncol = n2)
  
  for(i in 1:n1){
    for(j in 1:n2){
      sigma.ij    <- sigma(locs1[i, ], deg.aniso) * sigma(locs2[j, ], deg.aniso)
      kernel.ij   <- ( aniso_mat(locs1[i, ], deg.aniso) + aniso_mat(locs2[j, ], deg.aniso) ) / 2 
      smooth.ij   <- ( smoothness(locs1[i, ]) + smoothness(locs2[j, ]) ) / 2
      
      q.ij <- as.numeric(crossprod( locs1[i, ] - locs2[j, ], solve(kernel.ij, locs1[i, ] - locs2[j, ]) ))
      
      mat.cov[i,j] <- sigma.ij * fields::Matern( sqrt(q.ij), nu = smooth.ij) / sqrt( determinant(kernel.ij, logarithm = F)[[1]][1] )
    }
  }
  
  return(mat.cov)
}

# n         <- 15^2
# locs      <- matrix(runif(n * 2, 0, 1), n, 2)
# 
# Sigma.new <- matern_ns(locs1 = locs, locs2 = NULL, deg.aniso = 10)
# 
# cov.iso       <- function(locs, covparms) covparms[1] * exp(-fields::rdist(locs) / covparms[2])
# cov.aniso     <- function(locs, covparms) covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[ ,2])) / covparms[2])
# covparms      <- c(1, 0.1, 10)
# 
# Sigma.old <- cov.aniso(locs, covparms)
# 
# sum((Sigma.old - Sigma.new)^2)

####################################################################
#### simulation function
####################################################################

simulation <- function(n = 30^2, m = 30, deg.aniso = 10, covparms = c(1)) {

  locs      <- matrix(runif(n * 2, 0, 1), n, 2)
  Sigma     <- matern_ns(locs1 = locs, locs2 = NULL, deg.aniso = deg.aniso)
  y         <- as.numeric(t(chol(Sigma)) %*% rnorm(n))
  
  ### specify vecchia approximations
  approx <- list()
  
  # standard vecchia with maxmin ordering
  approx[[1]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "maxmin", which.coord = NULL, cond.yz='y', conditioning = "NN")
  # standard vecchia with x coord ordering
  approx[[2]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "coord", which.coord = 1, cond.yz='y', conditioning = "NN")
  # standard vecchia with y coord ordering
  approx[[3]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "coord", which.coord = 2, cond.yz='y', conditioning = "NN")
  # euclidean-based ordering + euclidean-based NN conditioning
  approx[[4]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma, covparms = covparms)
  # euclidean-based ordering + correlation-based NN conditioning
  approx[[5]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma, covparms = covparms)
  # correlation-based ordering + euclidean-based NN conditioning
  approx[[6]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma, covparms = covparms)
  # correlation-based ordering + correlation-based NN conditioning
  approx[[7]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma, covparms = covparms)
  # x-coordinate-based ordering + correlation-based NN conditioning
  approx[[8]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(1), def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma, covparms = covparms)
  # y-coordinate-based ordering + correlation-based NN conditioning
  approx[[9]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(2), def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma, covparms = covparms)
  # (x+y)-coordinate-based ordering + correlation-based NN conditioning
  approx[[10]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(1, 2), def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel = Sigma, covparms = covparms)
  # (x+y)-coordinate-based ordering + euclidean-based NN conditioning
  approx[[11]]           <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "coord", coordinate = c(1, 2), def.dist = NULL, ordering.method = "euclidean", initial.pt = NULL, conditioning = "NN", conditioning.method = "euclidean", covmodel = Sigma, covparms = covparms)
  # standard vecchia with maxmin ordering
  approx[[12]]           <- vecchia_specify_adjusted(locs = locs, m = m, ordering = "coord", which.coord = c(1, 2), cond.yz='y', conditioning = "NN")
  
  ### compute approximate covariance matrices
  n.approx    <- length(approx)
  Sigma.hat   <- list()
  kls         <- rep(NA, n.approx)
  for(i in 1:n.approx){
    
    Sigma.ord       <- matern_ns(locs1 = approx[[i]]$locsord, locs2 = NULL, deg.aniso = deg.aniso) # true cov in appropriate ordering
    
    U               <- createU(approx[[i]], c(1, 1, 1), 0, covmodel = Sigma.ord)$U
    revord          <- order(approx[[i]]$ord)
    Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    
    kls[i]          <- kldiv(Sigma, Sigma.hat[[i]])
  }
  
  result                  <- list()
  result$n                <- n
  result$m                <- m
  result$covparms         <- covparms
  result$locs             <- locs
  result$Sigma            <- Sigma
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
  
  xlabel2 <- sort(unique(vdat2$scale))
  plot2   <- ggplot(vdat2, aes(x=scale, y = log10(KL), col = method)) + 
    geom_point(aes(shape = method), size = 2) + 
    geom_line(size = 1, alpha = alpha.value) +
    ylab('log10(KL)') + 
    scale_x_continuous(name = 'a', limits=range(xlabel2), breaks=xlabel2) +
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
#### simulation 1
####################################################################

cand.m            <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45) ; n.cand.m <- length(cand.m)
cand.scale        <- c(1, 5, 10, 15, 20, 25) ; n.cand.scale <- length(cand.scale)
sim1              <- list()

# # small case
# cand.m            <- c(10, 20, 30) ; n.cand.m <- length(cand.m)
# cand.scale        <- c(5, 10) ; n.cand.scale <- length(cand.scale)
# sim1              <- list()

cand.all            <- expand.grid(cand.m, cand.scale)
cand.all            <- cbind(seq(nrow(cand.all)), cand.all)
colnames(cand.all)  <- c('index', 'm', 'scale')
n.cand.all          <- n.cand.m * n.cand.scale

no_cores            <- parallel::detectCores() - 2
cl                  <- parallel::makeCluster(no_cores)

doParallel::registerDoParallel(cl)
sim1 <- foreach(m = cand.all$m, a = cand.all$scale, .export = c("order_coordinate", "aniso_mat", "conditioning_nn", "correlation", "corrvecchia_knownCovparms", "distance_correlation", "kldiv", "matern_ns", "order_maxmin_correlation", "order_maxmin_correlation_old", "order_maxmin_euclidean", "simulation", "smoothness", "vecchia_specify_adjusted"), .packages='GPvecchia') %dopar% simulation(n = 30^2, m = m, deg.aniso = a, covparms = c(1))
parallel::stopCluster(cl)

# for(i in cand.all$index){
# 
#   start.time      <- proc.time()
#   sim1[[i]]       <- simulation(30^2, m = cand.all[i, 2], deg.aniso = cand.all[i, 3], covparms = c(1))
#   end.time        <- proc.time()
# 
#   proctime        <- end.time - start.time
# 
#   print(paste0("simulation ", i, " is done. [ ", proctime[3], "s ]"))
# }

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

set.scale   <- 10
ind         <- cand.all$scale == set.scale
vis.dat1    <- data.frame(kls.xcoord.eucord.euccond[ind], kls.ycoord.eucord.euccond[ind], kls.maxmin.eucord.euccond[ind], kls.xycoord.eucord.corcond[ind], kls.maxmin.corord.euccond[ind], kls.maxmin.corord.corcond[ind])
vis.dat1    <- vis.dat1[, order(colnames(vis.dat1))]
vis.dat1    <- cbind(rep(cand.m, times = ncol(vis.dat1)), tidyr::gather(vis.dat1))
colnames(vis.dat1) <- c("m", "method", "KL")
head(vis.dat1)

set.m       <- 30
ind         <- cand.all$m == set.m
vis.dat2    <- data.frame(kls.xcoord.eucord.euccond[ind], kls.ycoord.eucord.euccond[ind], kls.maxmin.eucord.euccond[ind], kls.xycoord.eucord.corcond[ind], kls.maxmin.corord.euccond[ind], kls.maxmin.corord.corcond[ind])
vis.dat2    <- vis.dat2[, order(colnames(vis.dat2))]
vis.dat2    <- cbind(rep(cand.scale, times = ncol(vis.dat2)), tidyr::gather(vis.dat2))
colnames(vis.dat2) <- c("scale", "method", "KL")
head(vis.dat2)

kls.legend <- c("C-Maxmin + C-NN", "C-Maxmin + E-NN", "E-Maxmin + E-NN", "X-Coord + E-NN", "(X+Y)-Coord + C-NN", "Y-Coord + E-NN")
vis_arrange(vdat1 = vis.dat1, vdat2 = vis.dat2, combined.legend = kls.legend, color.pal = brewer.pal(6, "Set1"), shape.pal = c(16, 17, 15, 18, 8, 13), alpha.value = 0.7, size.legend = 14, size.lab = 14, size.text = 12)

# save(sim1, cand.all, vis.dat1, vis.dat2, kls.legend, file='2_corrvecchia/sim_anisotropy_1.RData')
# load(file='2_corrvecchia/sim_anisotropy_1.RData')


