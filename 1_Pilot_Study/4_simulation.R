##########################################################################
#####
##### Author: Myeongjong Kang (kmj.stat@gmail.com)
#####
##### Description: 
#####
##########################################################################

# library(fields)

cov.iso       <- function(locs, range) exp(-rdist(locs)/range)
cov.aniso     <- function(locs, range, scale.aniso) exp(-rdist(cbind(locs[,1]*scale.aniso,locs[,2]))/range)

simulation    <- function(ordering = "maxmin", which.coord = NULL, conditioning = "NN", 
                          m = 10, scale.aniso = 10,
                          n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0) {
  
  ##### grid of locations: n, and #####
  grid.oneside    <- seq(0,1,length=sqrt(n))
  locs            <- as.matrix(expand.grid(grid.oneside,grid.oneside))
  locs.trans      <- cbind(locs[,1]*scale.aniso,locs[,2])
  
  #####  true.cov #####
  Sigma           <- cov.aniso(locs, range, scale.aniso)
  
  ##### compute KL divergence #####
  
  # specify approximations
  approx                <- list()
  # "standard" vecchia
  approx[[1]]           <- vecchia_specify_adjusted(locs, m, ordering = ordering, which.coord = which.coord, cond.yz='y', conditioning = conditioning)
  # correlation-based vecchia
  approx[[2]]           <- vecchia_specify_adjusted(locs.trans, m, ordering = ordering, which.coord = which.coord, cond.yz='y', conditioning = conditioning)
  approx[[2]]$locsord   <- locs[approx[[2]]$ord,]
  
  kls         <- rep(-9999, 2)
  Sigma.hat   <- list()
  for(i in 1:2){
    Sigma.ord       <- cov.aniso(approx[[i]]$locsord, range, scale.aniso) # true cov in appropriate ordering
    U               <- createU(approx[[i]], covparms, nugget, covmodel = Sigma.ord)$U
    revord          <- order(approx[[i]]$ord)
    Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    kls[i]          <- kldiv(Sigma,Sigma.hat[[i]])
  }

  ##### Return #####
  
  result                  <- list()
  result$ordering         <- ordering
  result$which.coord      <- which.coord
  result$conditioning     <- conditioning
  result$m                <- m
  result$n                <- n
  result$range            <- range
  result$scale.aniso      <- scale.aniso
  result$covparms         <- covparms
  result$nugget           <- nugget
  result$locs             <- locs
  result$locs.trans       <- locs.trans
  result$approx           <- approx
  result$kls              <- kls
  result$Sigma            <- Sigma
  result$Sigma.hat        <- Sigma.hat

  return(result)  
}

# sim.maxmin  <- simulation(ordering = "maxmin", which.coord = NULL, conditioning = "NN", m = 10, scale.aniso = 10, n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)
# sim.xcoord  <- simulation(ordering = "coord", which.coord = 1, conditioning = "NN", m = 10, scale.aniso = 10, n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)
# sim.ycoord  <- simulation(ordering = "coord", which.coord = 2, conditioning = "NN", m = 10, scale.aniso = 10, n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0)

repeated_simulation <- function(ordering = "maxmin", which.coord = NULL, conditioning = "NN", 
                                ms = c(1, 2, 5, 7, 10, 12, 15, 20, 25, 30, 35, 40), scales = c(1, 3, 7, 9, 11, 13),
                                n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0) {
  
  len.ms      <- length(ms)
  len.scales  <- length(scales)
  
  kls         <- list()
  kls[[1]]    <- matrix(-9999, nrow = len.ms, ncol = len.scales)
  kls[[2]]    <- matrix(-9999, nrow = len.ms, ncol = len.scales)
  
  for(i.m in 1:len.ms ){
    for(i.s in 1:len.scales){
      sim <- simulation(ordering = ordering, which.coord = which.coord, conditioning = conditioning, m = ms[i.m], scale.aniso = scales[i.s], n = n, range = range, covparms = covparms, nugget = nugget)
      kls[[1]][i.m,i.s] <- sim$kls[1]
      kls[[2]][i.m,i.s] <- sim$kls[2]
    }
  }
  
  result <- list()
  result$ordering         <- ordering
  result$which.coord      <- which.coord
  result$conditioning     <- conditioning
  result$ms               <- ms
  result$scales           <- scales
  result$n                <- n
  result$range            <- range
  result$covparms         <- covparms
  result$nugget           <- nugget
  result$kls              <- kls
  
  return(result)
}

simulation_nongrid <- function(ordering = "maxmin", which.coord = NULL, conditioning = "NN", 
                               m = 10, scale.aniso = 10,
                               n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0) {
  
  ##### grid of locations: n, and #####
  locs            <- matrix(runif(2 * n, 0, 1), n, 2)
  colnames(locs)  <- c("Var1", "Var2")
  locs.trans      <- cbind(locs[,1]*scale.aniso,locs[,2])
  
  #####  true.cov #####
  Sigma           <- cov.aniso(locs, range, scale.aniso)
  
  ##### compute KL divergence #####
  
  # specify approximations
  approx                <- list()
  # "standard" vecchia
  approx[[1]]           <- vecchia_specify_adjusted(locs, m, ordering = ordering, which.coord = which.coord, cond.yz='y', conditioning = conditioning)
  # correlation-based vecchia
  approx[[2]]           <- vecchia_specify_adjusted(locs.trans, m, ordering = ordering, which.coord = which.coord, cond.yz='y', conditioning = conditioning)
  approx[[2]]$locsord   <- locs[approx[[2]]$ord,]
  
  kls         <- rep(-9999, 2)
  Sigma.hat   <- list()
  for(i in 1:2){
    Sigma.ord       <- cov.aniso(approx[[i]]$locsord, range, scale.aniso) # true cov in appropriate ordering
    U               <- createU(approx[[i]], covparms, nugget, covmodel = Sigma.ord)$U
    revord          <- order(approx[[i]]$ord)
    Sigma.hat[[i]]  <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]
    kls[i]          <- kldiv(Sigma,Sigma.hat[[i]])
  }
  
  ##### Return #####
  
  result                  <- list()
  result$ordering         <- ordering
  result$which.coord      <- which.coord
  result$conditioning     <- conditioning
  result$m                <- m
  result$n                <- n
  result$range            <- range
  result$scale.aniso      <- scale.aniso
  result$covparms         <- covparms
  result$nugget           <- nugget
  result$locs             <- locs
  result$locs.trans       <- locs.trans
  result$approx           <- approx
  result$kls              <- kls
  result$Sigma            <- Sigma
  result$Sigma.hat        <- Sigma.hat
  
  return(result)  
}

repeated_simulation_nongrid <- function(ordering = "maxmin", which.coord = NULL, conditioning = "NN", 
                                        ms = c(1, 2, 5, 7, 10, 12, 15, 20, 25, 30, 35, 40), scales = c(1, 3, 7, 9, 11, 13),
                                        n = 30^2, range = 0.1, covparms = c(1, 0.1, 0.5), nugget = 0) {
  
  len.ms      <- length(ms)
  len.scales  <- length(scales)
  
  kls         <- list()
  kls[[1]]    <- matrix(-9999, nrow = len.ms, ncol = len.scales)
  kls[[2]]    <- matrix(-9999, nrow = len.ms, ncol = len.scales)
  
  for(i.m in 1:len.ms ){
    for(i.s in 1:len.scales){
      sim <- simulation_nongrid(ordering = ordering, which.coord = which.coord, conditioning = conditioning, m = ms[i.m], scale.aniso = scales[i.s], n = n, range = range, covparms = covparms, nugget = nugget)
      kls[[1]][i.m,i.s] <- sim$kls[1]
      kls[[2]][i.m,i.s] <- sim$kls[2]
    }
  }
  
  result <- list()
  result$ordering         <- ordering
  result$which.coord      <- which.coord
  result$conditioning     <- conditioning
  result$ms               <- ms
  result$scales           <- scales
  result$n                <- n
  result$range            <- range
  result$covparms         <- covparms
  result$nugget           <- nugget
  result$kls              <- kls
  
  return(result)
}