####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: 
###
####################################################################################

rm(list = ls())

set.seed(05042021)

library(correlationVecchia) ; library(ggplot2) ; library(dplyr) ; library(RColorBrewer) ; library(gridExtra)

##### 1. Randomly selected locations ###############################################

nsim            <- 100
covparms        <- c(1, 0.1, 1.0, 0.5) # covparms        <- c(1, 0.1, 0.2, 0.5)
nugget          <- 0.0001

n               <- 30^2
n.pred          <- 10^2
m               <- c(5, 10, 20, 30, 40)

#####

# process         <- generate_gp_spacetime(nsim = nsim, n = n, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)
# process.pred    <- generate_gp_spacetime(nsim = nsim, n = n.pred, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)
# 
# locs <- list()
# for(i in 1:nsim) locs[[i]] <- process$sim[[i]]$locs
# 
# locs.pred <- list()
# for(i in 1:nsim) locs.pred[[i]] <- process.pred$sim[[i]]$locs

#####

process         <- generate_gp_spacetime(nsim = nsim, n = n, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)

locs            <- list()
locs.pred       <- list()
for(i in 1:nsim) {
  
  locs[[i]]       <- process$sim[[i]]$locs
  
  ind.unobs       <- sort(sample(seq(n), n.pred))
  
  locs.pred[[i]]  <- locs[[i]][ind.unobs, , drop = FALSE]
  locs[[i]]       <- locs[[i]][-ind.unobs, , drop = FALSE]
}

plot(locs[[1]][, 1:2], col = "red", main = "Random locations")
points(locs.pred[[1]][, 1:2], col = "blue")

#####

z <- list() ; z.pred <- list()
for(i in 1:nsim) {
  
  covmat          <- cov_matern_spacetime(locs = rbind(locs[[i]], locs.pred[[i]]), covparms = covparms) + nugget * diag(nrow(locs[[i]]) + nrow(locs.pred[[i]]))
  
  z.all           <- as.numeric(t(chol(x = covmat)) %*% as.matrix(rnorm(nrow(locs[[i]]) + nrow(locs.pred[[i]]))))
  z[[i]]          <- z.all[seq(nrow(locs[[i]]))]
  z.pred[[i]]     <- z.all[seq(nrow(locs.pred[[i]])) + nrow(locs[[i]])]
}

Sys.time() ######

result.allrandom       <- list()
result.allrandom[[1]]  <- list(n = n, n.pred = n.pred, m = m, covparms = covparms, nugget = nugget, locs = locs, locs.pred = locs.pred, z = z, z.pred = z.pred)

simout <- list()
for(i in 1:5) simout[[i]] <- list(mspe = matrix(NA, nrow = nsim, ncol = length(m)), logscore = matrix(NA, nrow = nsim, ncol = length(m)))

names(simout) <- c("b1", "b2", "b3", "cc", "eu")

for(i in 1:nsim) {
  
  for(j in 1:length(m)) {
    
    out.baseline1   <- prediction_baseline_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covparms = covparms, nuggets = nugget, method = "b1", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.baseline2   <- prediction_baseline_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covparms = covparms, nuggets = nugget, method = "b2", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.baseline3   <- prediction_baseline_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covparms = covparms, nuggets = nugget, method = "b3", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.euclidean   <- prediction_corrvecchia_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], initial.pt = NULL, covparms = covparms, nuggets = nugget, method = "euclidean", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.corrvecchia <- prediction_corrvecchia_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], initial.pt = NULL, covparms = covparms, nuggets = nugget, method = "correlation", predcond.method = "general", var.exact = TRUE, return.values = "all")
    
    simout$b1$mspe[i, j] <- mean((out.baseline1$predict$mu.pred - z.pred[[i]])^2)
    simout$b2$mspe[i, j] <- mean((out.baseline2$predict$mu.pred - z.pred[[i]])^2)
    simout$b3$mspe[i, j] <- mean((out.baseline3$predict$mu.pred - z.pred[[i]])^2)
    simout$cc$mspe[i, j] <- mean((out.corrvecchia$predict$mu.pred - z.pred[[i]])^2)
    simout$eu$mspe[i, j] <- mean((out.euclidean$predict$mu.pred - z.pred[[i]])^2)
    
    simout$b1$logscore[i, j] <- mean(logscore(out.baseline1$predict$mu.pred, out.baseline1$predict$var.pred, z.pred[[i]]))
    simout$b2$logscore[i, j] <- mean(logscore(out.baseline2$predict$mu.pred, out.baseline2$predict$var.pred, z.pred[[i]]))
    simout$b3$logscore[i, j] <- mean(logscore(out.baseline3$predict$mu.pred, out.baseline3$predict$var.pred, z.pred[[i]]))
    simout$cc$logscore[i, j] <- mean(logscore(out.corrvecchia$predict$mu.pred, out.corrvecchia$predict$var.pred, z.pred[[i]]))
    simout$eu$logscore[i, j] <- mean(logscore(out.euclidean$predict$mu.pred, out.euclidean$predict$var.pred, z.pred[[i]]))
  }
  
  print(i) ; print(Sys.time())
}

result.allrandom[[2]]  <- simout

mspe.b1 <- colMeans(result.allrandom[[2]]$b1$mspe)
mspe.b2 <- colMeans(result.allrandom[[2]]$b2$mspe)
mspe.b3 <- colMeans(result.allrandom[[2]]$b3$mspe)
mspe.cc <- colMeans(result.allrandom[[2]]$cc$mspe)
mspe.eu <- colMeans(result.allrandom[[2]]$eu$mspe)

logs.b1 <- colMeans(result.allrandom[[2]]$b1$logscore)
logs.b2 <- colMeans(result.allrandom[[2]]$b2$logscore)
logs.b3 <- colMeans(result.allrandom[[2]]$b3$logscore)
logs.cc <- colMeans(result.allrandom[[2]]$cc$logscore)
logs.eu <- colMeans(result.allrandom[[2]]$eu$logscore)

result.allrandom[[3]] <- list(mspe = list(b1 = mspe.b1, b2 = mspe.b2, b3 = mspe.b3, cc = mspe.cc, eu = mspe.eu), logscore = list(b1 = logs.b1, b2 = logs.b2, b3 = logs.b3, cc = logs.cc, eu = logs.eu))

names(result.allrandom) <- c("setting", "simout", "output")
save(result.allrandom, file = "sim_prediction_allrandom_100outof900_05042021.RData")

##### 2. Monitoring-station locations ##############################################

nsim            <- 100
covparms        <- c(1, 0.1, 1.0, 0.5) # covparms        <- c(1, 0.1, 0.2, 0.5)
nugget          <- 0.0001

n               <- 30^2
t               <- 36
n.pred          <- 5 * t
m               <- c(5, 10, 20, 30, 40)

#####

# process         <- generate_gp_spacetime(nsim = nsim, n = n/t, d = 2, t.len = t, method.locs = "space.random.time.grid", covmodel = cov_matern_spacetime, covparms = covparms)
# process.pred    <- generate_gp_spacetime(nsim = nsim, n = n.pred, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)
# 
# locs <- list()
# for(i in 1:nsim) locs[[i]] <- process$sim[[i]]$locs
# 
# locs.pred <- list()
# for(i in 1:nsim) locs.pred[[i]] <- process.pred$sim[[i]]$locs

#####

process         <- generate_gp_spacetime(nsim = nsim, n = n/t, d = 2, t.len = t, method.locs = "space.random.time.grid", covmodel = cov_matern_spacetime, covparms = covparms)

locs            <- list()
locs.pred       <- list()
for(i in 1:nsim) {
  
  locs[[i]]       <- as.data.frame(process$sim[[i]]$locs)
  
  locs.spatial    <- dplyr::distinct(locs[[i]][, 1:2])
  locs.spatial    <- cbind(locs.spatial, seq(nrow(locs.spatial)))
  
  locs[[i]]       <- dplyr::left_join(locs[[i]], locs.spatial, by = c("V1", "V2"))
  colnames(locs[[i]]) <- c("coord1", "coord2", "time", "index")
  
  ind.unobs       <- sort(sample(seq(n/t), n.pred/t))
  bool.unobs      <- locs[[i]]$index %in% ind.unobs
  
  locs.pred[[i]]  <- locs[[i]][bool.unobs, , drop = FALSE]
  locs[[i]]       <- locs[[i]][!bool.unobs, , drop = FALSE]
  
  locs.pred[[i]]  <- as.matrix(locs.pred[[i]][, c("coord1", "coord2", "time")])
  locs[[i]]       <- as.matrix(locs[[i]][, c("coord1", "coord2", "time")])
}

plot(locs[[1]][, 1:2], col = "red", main = "Monitoring stations")
points(locs.pred[[1]][, 1:2], col = "blue")

#####

z <- list() ; z.pred <- list()
for(i in 1:nsim) {
  
  covmat          <- cov_matern_spacetime(locs = rbind(locs[[i]], locs.pred[[i]]), covparms = covparms) + nugget * diag(nrow(locs[[i]]) + nrow(locs.pred[[i]]))
  
  z.all           <- as.numeric(t(chol(x = covmat)) %*% as.matrix(rnorm(nrow(locs[[i]]) + nrow(locs.pred[[i]]))))
  z[[i]]          <- z.all[seq(nrow(locs[[i]]))]
  z.pred[[i]]     <- z.all[seq(nrow(locs.pred[[i]])) + nrow(locs[[i]])]
}

Sys.time() ######

result.monitoring       <- list()
result.monitoring[[1]]  <- list(n = n, n.pred = n.pred, m = m, covparms = covparms, nugget = nugget, locs = locs, locs.pred = locs.pred, z = z, z.pred = z.pred)

simout <- list()
for(i in 1:5) simout[[i]] <- list(mspe = matrix(NA, nrow = nsim, ncol = length(m)), logscore = matrix(NA, nrow = nsim, ncol = length(m)))

names(simout) <- c("b1", "b2", "b3", "cc", "eu")

for(i in 1:nsim) {
  
  for(j in 1:length(m)) {
    
    out.baseline1   <- prediction_baseline_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covparms = covparms, nuggets = nugget, method = "b1", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.baseline2   <- prediction_baseline_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covparms = covparms, nuggets = nugget, method = "b2", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.baseline3   <- prediction_baseline_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covparms = covparms, nuggets = nugget, method = "b3", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.euclidean   <- prediction_corrvecchia_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], initial.pt = NULL, covparms = covparms, nuggets = nugget, method = "euclidean", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.corrvecchia <- prediction_corrvecchia_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], initial.pt = NULL, covparms = covparms, nuggets = nugget, method = "correlation", predcond.method = "general", var.exact = TRUE, return.values = "all")
    
    simout$b1$mspe[i, j] <- mean((out.baseline1$predict$mu.pred - z.pred[[i]])^2)
    simout$b2$mspe[i, j] <- mean((out.baseline2$predict$mu.pred - z.pred[[i]])^2)
    simout$b3$mspe[i, j] <- mean((out.baseline3$predict$mu.pred - z.pred[[i]])^2)
    simout$cc$mspe[i, j] <- mean((out.corrvecchia$predict$mu.pred - z.pred[[i]])^2)
    simout$eu$mspe[i, j] <- mean((out.euclidean$predict$mu.pred - z.pred[[i]])^2)
    
    simout$b1$logscore[i, j] <- mean(logscore(out.baseline1$predict$mu.pred, out.baseline1$predict$var.pred, z.pred[[i]]))
    simout$b2$logscore[i, j] <- mean(logscore(out.baseline2$predict$mu.pred, out.baseline2$predict$var.pred, z.pred[[i]]))
    simout$b3$logscore[i, j] <- mean(logscore(out.baseline3$predict$mu.pred, out.baseline3$predict$var.pred, z.pred[[i]]))
    simout$cc$logscore[i, j] <- mean(logscore(out.corrvecchia$predict$mu.pred, out.corrvecchia$predict$var.pred, z.pred[[i]]))
    simout$eu$logscore[i, j] <- mean(logscore(out.euclidean$predict$mu.pred, out.euclidean$predict$var.pred, z.pred[[i]]))
  }
  
  print(i) ; print(Sys.time())
}

result.monitoring[[2]]  <- simout

mspe.b1 <- colMeans(result.monitoring[[2]]$b1$mspe)
mspe.b2 <- colMeans(result.monitoring[[2]]$b2$mspe)
mspe.b3 <- colMeans(result.monitoring[[2]]$b3$mspe)
mspe.cc <- colMeans(result.monitoring[[2]]$cc$mspe)
mspe.eu <- colMeans(result.monitoring[[2]]$eu$mspe)

logs.b1 <- colMeans(result.monitoring[[2]]$b1$logscore)
logs.b2 <- colMeans(result.monitoring[[2]]$b2$logscore)
logs.b3 <- colMeans(result.monitoring[[2]]$b3$logscore)
logs.cc <- colMeans(result.monitoring[[2]]$cc$logscore)
logs.eu <- colMeans(result.monitoring[[2]]$eu$logscore)

result.monitoring[[3]] <- list(mspe = list(b1 = mspe.b1, b2 = mspe.b2, b3 = mspe.b3, cc = mspe.cc, eu = mspe.eu), logscore = list(b1 = logs.b1, b2 = logs.b2, b3 = logs.b3, cc = logs.cc, eu = logs.eu))

names(result.monitoring) <- c("setting", "simout", "output")
save(result.monitoring, file = "sim_prediction_monitoring_100outof900_05042021.RData")

##### 3. Satellite-footprint locations #############################################

nsim            <- 100
covparms        <- c(1, 0.1, 1.0, 0.5) # covparms        <- c(1, 0.1, 0.2, 0.5)
nugget          <- 0.0001

n               <- 30^2
n.pred          <- 10^2
m               <- c(5, 10, 20, 30, 40)

#####

# process         <- generate_gp_spacetime(nsim = nsim, n = n, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_matern_spacetime, covparms = covparms)
# process.pred    <- generate_gp_spacetime(nsim = nsim, n = n.pred, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)
# 
# locs <- list()
# for(i in 1:nsim) locs[[i]] <- process$sim[[i]]$locs
# 
# locs.pred <- list()
# for(i in 1:nsim) locs.pred[[i]] <- process.pred$sim[[i]]$locs

#####

process         <- generate_gp_spacetime(nsim = nsim, n = n, d = 2, t.len = 1, method.locs = "satellite", covmodel = cov_matern_spacetime, covparms = covparms)
process.pred    <- generate_gp_spacetime(nsim = nsim, n = n.pred, d = 2, t.len = 1, method.locs = "all.random", covmodel = cov_matern_spacetime, covparms = covparms)

locs            <- list()
locs.pred       <- list()
for(i in 1:nsim) {
  
  locs[[i]]       <- process$sim[[i]]$locs
  
  ind.unobs       <- sample(seq(n), 1)
  ind.unobs       <- seq(from = ind.unobs, by = 1, length.out = n.pred)
  
  ind.unobs[ind.unobs > n] <- ind.unobs[ind.unobs > n] - n
  ind.unobs       <- sort(ind.unobs)
  
  locs.pred[[i]]  <- locs[[i]][ind.unobs, , drop = FALSE]
  locs[[i]]       <- locs[[i]][-ind.unobs, , drop = FALSE]
}

plot(locs[[1]][, 1:2], col = "red", main = "Remote sensing")
points(locs.pred[[1]][, 1:2], col = "blue")

#####

z <- list() ; z.pred <- list()
for(i in 1:nsim) {
  
  covmat          <- cov_matern_spacetime(locs = rbind(locs[[i]], locs.pred[[i]]), covparms = covparms) + nugget * diag(nrow(locs[[i]]) + nrow(locs.pred[[i]]))
  
  z.all           <- as.numeric(t(chol(x = covmat)) %*% as.matrix(rnorm(nrow(locs[[i]]) + nrow(locs.pred[[i]]))))
  z[[i]]          <- z.all[seq(nrow(locs[[i]]))]
  z.pred[[i]]     <- z.all[seq(nrow(locs.pred[[i]])) + nrow(locs[[i]])]
}

Sys.time() ######

result.satellite       <- list()
result.satellite[[1]]  <- list(n = n, n.pred = n.pred, m = m, covparms = covparms, nugget = nugget, locs = locs, locs.pred = locs.pred, z = z, z.pred = z.pred)

simout <- list()
for(i in 1:5) simout[[i]] <- list(mspe = matrix(NA, nrow = nsim, ncol = length(m)), logscore = matrix(NA, nrow = nsim, ncol = length(m)))

names(simout) <- c("b1", "b2", "b3", "cc", "eu")

for(i in 1:nsim) {
  
  for(j in 1:length(m)) {
    
    out.baseline1   <- prediction_baseline_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covparms = covparms, nuggets = nugget, method = "b1", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.baseline2   <- prediction_baseline_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covparms = covparms, nuggets = nugget, method = "b2", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.baseline3   <- prediction_baseline_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], coordinate = NULL, covparms = covparms, nuggets = nugget, method = "b3", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.euclidean   <- prediction_corrvecchia_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], initial.pt = NULL, covparms = covparms, nuggets = nugget, method = "euclidean", predcond.method = "general", var.exact = TRUE, return.values = "all")
    out.corrvecchia <- prediction_corrvecchia_for_spacetime(z = z[[i]], locs = locs[[i]], locs.pred = locs.pred[[i]], m = m[j], initial.pt = NULL, covparms = covparms, nuggets = nugget, method = "correlation", predcond.method = "general", var.exact = TRUE, return.values = "all")
    
    simout$b1$mspe[i, j] <- mean((out.baseline1$predict$mu.pred - z.pred[[i]])^2)
    simout$b2$mspe[i, j] <- mean((out.baseline2$predict$mu.pred - z.pred[[i]])^2)
    simout$b3$mspe[i, j] <- mean((out.baseline3$predict$mu.pred - z.pred[[i]])^2)
    simout$cc$mspe[i, j] <- mean((out.corrvecchia$predict$mu.pred - z.pred[[i]])^2)
    simout$eu$mspe[i, j] <- mean((out.euclidean$predict$mu.pred - z.pred[[i]])^2)
    
    simout$b1$logscore[i, j] <- mean(logscore(out.baseline1$predict$mu.pred, out.baseline1$predict$var.pred, z.pred[[i]]))
    simout$b2$logscore[i, j] <- mean(logscore(out.baseline2$predict$mu.pred, out.baseline2$predict$var.pred, z.pred[[i]]))
    simout$b3$logscore[i, j] <- mean(logscore(out.baseline3$predict$mu.pred, out.baseline3$predict$var.pred, z.pred[[i]]))
    simout$cc$logscore[i, j] <- mean(logscore(out.corrvecchia$predict$mu.pred, out.corrvecchia$predict$var.pred, z.pred[[i]]))
    simout$eu$logscore[i, j] <- mean(logscore(out.euclidean$predict$mu.pred, out.euclidean$predict$var.pred, z.pred[[i]]))
  }
  
  print(i) ; print(Sys.time())
}

result.satellite[[2]]  <- simout

mspe.b1 <- colMeans(result.satellite[[2]]$b1$mspe)
mspe.b2 <- colMeans(result.satellite[[2]]$b2$mspe)
mspe.b3 <- colMeans(result.satellite[[2]]$b3$mspe)
mspe.cc <- colMeans(result.satellite[[2]]$cc$mspe)
mspe.eu <- colMeans(result.satellite[[2]]$eu$mspe)

logs.b1 <- colMeans(result.satellite[[2]]$b1$logscore)
logs.b2 <- colMeans(result.satellite[[2]]$b2$logscore)
logs.b3 <- colMeans(result.satellite[[2]]$b3$logscore)
logs.cc <- colMeans(result.satellite[[2]]$cc$logscore)
logs.eu <- colMeans(result.satellite[[2]]$eu$logscore)

result.satellite[[3]] <- list(mspe = list(b1 = mspe.b1, b2 = mspe.b2, b3 = mspe.b3, cc = mspe.cc, eu = mspe.eu), logscore = list(b1 = logs.b1, b2 = logs.b2, b3 = logs.b3, cc = logs.cc, eu = logs.eu))

names(result.satellite) <- c("setting", "simout", "output")
save(result.satellite, file = "sim_prediction_satellite_100outof900_05042021.RData")

##### 4. rm ########################################################################

rm(list = ls())

##### 5. Visualization #############################################################

load("sim_prediction_allrandom_100outof900_05042021.RData")
# load("sim_prediction_allrandom_100and900_05042021.RData")

b1 <- log10(-result.allrandom$output$logscore$b1)
b2 <- log10(-result.allrandom$output$logscore$b2)
b3 <- log10(-result.allrandom$output$logscore$b3)
cc <- log10(-result.allrandom$output$logscore$cc)
eu <- log10(-result.allrandom$output$logscore$eu)

vdat1 <- data.frame(index = seq(length(b1)), m = result.allrandom$setting$m, approx_1 = b1, approx_2 = b2, approx_3 = b3, approx_4 = cc, approx_5 = eu)

par(mar = c(5.1, 4.1 + 0.9, 4.1, 2.1))

plot(c(5, 10, 20, 30, 40), b1, type = 'o', ylim = range(c(b1, b2, b3, cc, eu)) + c(0, 0.1), col = "#984EA3", lwd = 3, xlab = "m", ylab = "log10(logscore)", main = "")
points(c(5, 10, 20, 30, 40), b2, type = 'o', col = "#4DAF4A", lwd = 3)
points(c(5, 10, 20, 30, 40), b3, type = 'o', col = "#377EB8", lwd = 5)
points(c(5, 10, 20, 30, 40), cc, type = 'o', col = "#E41A1C", lwd = 3)
points(c(5, 10, 20, 30, 40), eu, type = 'o', col = "gray70", lwd = 3)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "E-MM + E-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lwd = 3)

#####

load("sim_prediction_monitoring_100outof900_05042021.RData")
# load("sim_prediction_monitoring_100and900_05042021.RData")

b1 <- log10(-result.monitoring$output$logscore$b1)
b2 <- log10(-result.monitoring$output$logscore$b2)
b3 <- log10(-result.monitoring$output$logscore$b3)
cc <- log10(-result.monitoring$output$logscore$cc)
eu <- log10(-result.monitoring$output$logscore$eu)

vdat2 <- data.frame(index = seq(length(b1)), m = result.allrandom$setting$m, approx_1 = b1, approx_2 = b2, approx_3 = b3, approx_4 = cc, approx_5 = eu)

par(mar = c(5.1, 4.1 + 0.9, 4.1, 2.1))

plot(c(5, 10, 20, 30, 40), b1, type = 'o', ylim = range(c(b1, b2, b3, cc, eu)) + c(0, 0.025), col = "#984EA3", lwd = 3, xlab = "m", ylab = "log10(logscore)", main = "")
points(c(5, 10, 20, 30, 40), b2, type = 'o', col = "#4DAF4A", lwd = 3)
points(c(5, 10, 20, 30, 40), b3, type = 'o', col = "#377EB8", lwd = 3)
points(c(5, 10, 20, 30, 40), cc, type = 'o', col = "#E41A1C", lwd = 3)
points(c(5, 10, 20, 30, 40), eu, type = 'o', col = "gray70", lwd = 3)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "E-MM + E-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lwd = 3)

#####

load("sim_prediction_satellite_100outof900_05042021.RData")
# load("sim_prediction_satellite_100and900_05042021.RData")

b1 <- log10(-result.satellite$output$logscore$b1)
b2 <- log10(-result.satellite$output$logscore$b2)
b3 <- log10(-result.satellite$output$logscore$b3)
cc <- log10(-result.satellite$output$logscore$cc)
eu <- log10(-result.satellite$output$logscore$eu)

vdat3 <- data.frame(index = seq(length(b1)), m = result.allrandom$setting$m, approx_1 = b1, approx_2 = b2, approx_3 = b3, approx_4 = cc, approx_5 = eu)

par(mar = c(5.1, 4.1 + 0.9, 4.1, 2.1))

plot(c(5, 10, 20, 30, 40), b1, type = 'o', ylim = range(c(b1, b2, b3, cc, eu)) + c(0, 0.05), col = "#984EA3", lwd = 3, xlab = "m", ylab = "log10(logscore)", main = "")
points(c(5, 10, 20, 30, 40), b2, type = 'o', col = "#4DAF4A", lwd = 3)
points(c(5, 10, 20, 30, 40), b3, type = 'o', col = "#377EB8", lwd = 5)
points(c(5, 10, 20, 30, 40), cc, type = 'o', col = "#E41A1C", lwd = 3)
points(c(5, 10, 20, 30, 40), eu, type = 'o', col = "gray70", lwd = 3)
legend("topright", legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "E-MM + E-NN"), col = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), lwd = 3)

#####

rm(b1, b2, b3, cc, eu)

#####

vis <- vis_prediction(vdat1 = vdat1, vdat2 = vdat2, vdat3 = vdat3, legend = c("T-ord + T-NN", "T-ord + E-NN", "T-ord + C-NN", "C-MM + C-NN", "E-MM + E-NN"), color = c("#984EA3", "#4DAF4A", "#377EB8", "#E41A1C", "gray70"), shape = c(18, 15, 17, 16, 8))
ggplot2::ggsave("prediction_05042021.pdf", vis, width = 15.2, height = 5.7)
