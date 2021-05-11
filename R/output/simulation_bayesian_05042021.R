####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: 
###
####################################################################################

rm(list = ls())
# set.seed(05042021)
set.seed(05072021)

### load packages
library(correlationVecchia)
library(foreach)

##### ncores #######################################################################

no_cores  <- parallel::detectCores() - 2

##### cov_matern_spacetime #########################################################

ms              <- c(5, 10, 20, 30, 40)
approxs         <- c("b1", "b2", "b3", "cc")
candid          <- expand.grid(ms, approxs) ; colnames(candid) <- c("m", "approx")

n               <- 30^2
covparms        <- c(1, 0.1, 1.0, 0.5)
nugget          <- 0.00

### n = 900, monitoring station #######################################

process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "space.random.time.grid", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs

z               <- process$sim$sim1$y + nugget * rnorm(n) ; fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.spatial.n900.monitoring <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.05, 0.12), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.temporal.n900.monitoring <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.05, 3.05), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

### n = 900, satellite ################################################

process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs

z               <- process$sim$sim1$y + nugget * rnorm(n) # fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.spatial.n900.satellite <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.08, 0.12), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.temporal.n900.satellite <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.05, 3.05), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

### n = 900, all random locs ##########################################

process         <- generate_gp_spacetime(nsim = 2, n = n, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)
locs            <- process$sim$sim1$locs

z               <- process$sim$sim1$y + nugget * rnorm(n) ; fields::quilt.plot(locs[, 1], locs[, 2], z, main = "process")
# covmat          <- process$sim$sim1$covmat + nugget * diag(n)

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.spatial.n900.allrandom <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "spatialRange", N = 200, xlim = c(0.08, 0.16), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

cl <- parallel::makeCluster(no_cores) #######
doParallel::registerDoParallel(cl)

out.temporal.n900.allrandom <- foreach::foreach(i = seq(nrow(candid)), .packages = c("correlationVecchia", "GPvecchia", "GpGp", "mvtnorm")) %dopar% posterior_for_cov_matern_spacetime(z = z, locs = locs, covparms = covparms, nugget = nugget, m = candid$m[i], approx = candid$approx[i], target = "temporalRange", N = 200, xlim = c(0.05, 2.25), sdlog = 0.6)

parallel::stopCluster(cl)
Sys.time() ##################################

##### Save #########################################################################

save(candid, n, covparms, nugget, out.spatial.n900.allrandom, out.temporal.n900.allrandom, out.spatial.n900.monitoring, out.temporal.n900.monitoring, out.spatial.n900.satellite, out.temporal.n900.satellite, file = "bayesian_05072021.RData")

##### Visualization  1 #############################################################

par(mfrow = c(4, 5))
for(i in 1:length(out.spatial.n900.allrandom)) {
  temp <- out.spatial.n900.allrandom[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(4, 5))
for(i in 1:length(out.temporal.n900.allrandom)) {
  temp <- out.temporal.n900.allrandom[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(4, 5))
for(i in 1:length(out.spatial.n900.monitoring)) {
  temp <- out.spatial.n900.monitoring[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(4, 5))
for(i in 1:length(out.temporal.n900.monitoring)) {
  temp <- out.temporal.n900.monitoring[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(4, 5))
for(i in 1:length(out.spatial.n900.monitoring)) {
  temp <- out.spatial.n900.satellite[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

par(mfrow = c(4, 5))
for(i in 1:length(out.temporal.n900.monitoring)) {
  temp <- out.temporal.n900.satellite[[i]]
  matplot(temp$alpha, temp$post.norm, type = 'l', xlim = range(temp$alpha), lty = 1, lwd = 2, xlab = 'range', ylab = 'posterior density', main = paste0("m=", candid[i, 1], ", approx=", candid[i, 2]))
} ; rm(temp)
par(mfrow = c(1, 1))

##### Visualization  2 #############################################################

# m <- unique(candid$m)
# cols <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")
# 
# idx.list <- list()
# for(i in 1:length(m)) idx.list[[i]] <- which(candid$m == m[i])
# 
# par(mfrow = c(1, 5))
# for(i in 1:length(idx.list)) {
#   k <- 1
#   plot(out.spatial.n900.monitoring[[idx.list[[i]][1]]]$alpha, out.spatial.n900.monitoring[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
#   for(j in idx.list[[i]]) {
#     lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = min(k + 1, 4), lty = "dashed")
#     k <- k + 1
#   }
# }
# par(mfrow = c(1, 1))
# 
# par(mfrow = c(1, 5))
# for(i in 1:length(idx.list)) {
#   k <- 1
#   plot(out.temporal.n900.monitoring[[idx.list[[i]][1]]]$alpha, out.temporal.n900.monitoring[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
#   for(j in idx.list[[i]]) {
#     lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = min(k + 1, 4), lty = "dashed")
#     k <- k + 1
#   }
# }
# par(mfrow = c(1, 1))
# 
# par(mfrow = c(1, 5))
# for(i in 1:length(idx.list)) {
#   k <- 1
#   plot(out.spatial.n900.satellite[[idx.list[[i]][1]]]$alpha, out.spatial.n900.satellite[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
#   for(j in idx.list[[i]]) {
#     lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 2)
#     k <- k + 1
#   }
# }
# par(mfrow = c(1, 1))
# 
# par(mfrow = c(1, 5))
# for(i in 1:length(idx.list)) {
#   k <- 1
#   plot(out.temporal.n900.satellite[[idx.list[[i]][1]]]$alpha, out.temporal.n900.satellite[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
#   for(j in idx.list[[i]]) {
#     lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 2)
#     k <- k + 1
#   }
# }
# par(mfrow = c(1, 1))

##### Visualization  3 #############################################################

m <- unique(candid$m)
cols <- c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C")

idx.list <- list()
for(i in 1:length(m)) idx.list[[i]] <- which(candid$m == m[i])


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.spatial.n900.allrandom[[idx.list[[i]][1]]]$alpha, out.spatial.n900.allrandom[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topleft", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.spatial.n900.allrandom[[j]]$alpha, out.spatial.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
}
par(mfrow = c(1, 1))


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.temporal.n900.allrandom[[idx.list[[i]][1]]]$alpha, out.temporal.n900.allrandom[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.temporal.n900.allrandom[[j]]$alpha, out.temporal.n900.allrandom[[j]]$post.norm[, 2], col = cols[k], lwd = 2, lty = "solid")
  
}
par(mfrow = c(1, 1))


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {

  plot(out.spatial.n900.monitoring[[idx.list[[i]][1]]]$alpha, out.spatial.n900.monitoring[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topleft", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "twodash")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "longdash")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotdash")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.spatial.n900.monitoring[[j]]$alpha, out.spatial.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dashed")
  
}
par(mfrow = c(1, 1))


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.temporal.n900.monitoring[[idx.list[[i]][1]]]$alpha, out.temporal.n900.monitoring[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lwd = 3)
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "twodash")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "longdash")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotdash")
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.temporal.n900.monitoring[[j]]$alpha, out.temporal.n900.monitoring[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dashed")
  
}
par(mfrow = c(1, 1))

par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.spatial.n900.satellite[[idx.list[[i]][1]]]$alpha, out.spatial.n900.satellite[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lwd = 3)
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 1, lty = "solid")
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "twodash")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "longdash")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.spatial.n900.satellite[[j]]$alpha, out.spatial.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotdash")

}
par(mfrow = c(1, 1))


par(mfrow = c(1, 5))
for(i in 1:length(idx.list)) {
  
  plot(out.temporal.n900.satellite[[idx.list[[i]][1]]]$alpha, out.temporal.n900.satellite[[idx.list[[i]][1]]]$post.norm[, 1], type = "l", lwd = 3, xlab = "range", ylab = "density", main = paste0("m = ", m[i]))
  
  if(i == length(idx.list)) legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C"), lwd = 3)
  
  j <- idx.list[[i]][4]
  k <- 4
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 1, lty = "solid")
  
  j <- idx.list[[i]][1]
  k <- 1
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "twodash")
  
  j <- idx.list[[i]][2]
  k <- 2
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "longdash")
  
  j <- idx.list[[i]][3]
  k <- 3
  lines(out.temporal.n900.satellite[[j]]$alpha, out.temporal.n900.satellite[[j]]$post.norm[, 2], col = cols[k], lwd = 3, lty = "dotdash")
  
}
par(mfrow = c(1, 1))


