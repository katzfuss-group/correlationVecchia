####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
####################################################################################

rm(list = ls())

#########################################################################

set.seed(10112021)

library(dplyr)
library(GpGp)
library(correlationVecchia)

#########################################################################

load("DATA/joint_CRCM_NCEP_10152021_estimates.RData")

round(fit.b1$betahat, 5) ; round(fit.b1$covparms, 5)
round(fit.b2$betahat, 5) ; round(fit.b2$covparms, 5)
round(fit.b3$betahat, 5) ; round(fit.b3$covparms, 5)
round(fit.b4$betahat, 5) ; round(fit.b4$covparms, 5)
round(fit.b5$betahat, 5) ; round(fit.b5$covparms, 5)
round(fit.b6$betahat, 5) ; round(fit.b6$covparms, 5)
round(fit.b7$betahat, 5) ; round(fit.b7$covparms, 5)
round(fit.cc$betahat, 5) ; round(fit.cc$covparms, 5)

#########################################################################

colnames(df.joint)

locs        <- train.joint[, c("x1", "x2", "t", "d")] %>% as.matrix()
z           <- train.joint$z
X.train     <- cbind(1, train.joint$d) %>% as.matrix()

locs.pred   <- test.joint[, c("x1", "x2", "t", "d")] %>% as.matrix()
z.pred      <- test.joint$z
X.pred      <- cbind(1, test.joint$d) %>% as.matrix()

ns_obs      <- train.joint %>% pull(d) %>% table() %>% as.numeric()
ns_pred     <- test.joint %>% pull(d) %>% table() %>% as.numeric()

rm(train.joint, test.joint)

#########################################################################

comparison <- function(m)
{
  pre.b2  <- predictions_bs_mulv(approx = 2, ns_obs, ns_pred, fit = fit.b2, locs_pred = locs.pred, m = m, joint = FALSE, nsims = 0, predvar = TRUE, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.b4  <- predictions_bs_mulv(approx = 4, ns_obs, ns_pred, fit = fit.b4, locs_pred = locs.pred, m = m, joint = FALSE, nsims = 0, predvar = TRUE, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.b5  <- predictions_bs_mulv(approx = 5, ns_obs, ns_pred, fit = fit.b5, locs_pred = locs.pred, m = m, joint = FALSE, nsims = 0, predvar = TRUE, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.b6  <- predictions_bs_mulv(approx = 6, ns_obs, ns_pred, fit = fit.b6, locs_pred = locs.pred, m = m, joint = FALSE, nsims = 0, predvar = TRUE, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.b7  <- predictions_bs_mulv(approx = 7, ns_obs, ns_pred, fit = fit.b7, locs_pred = locs.pred, m = m, joint = FALSE, nsims = 0, predvar = TRUE, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.cc  <- predictions_cv_mulv(            ns_obs, ns_pred, fit = fit.cc, locs_pred = locs.pred, m = m, joint = FALSE, nsims = 0, predvar = TRUE, X_pred = X.pred, scale = 'parms', tol = 1e-8)

  mspe.b2 <- mean( (pre.b2$means - z.pred)^2 )
  mspe.b4 <- mean( (pre.b4$means - z.pred)^2 )
  mspe.b5 <- mean( (pre.b5$means - z.pred)^2 )
  mspe.b6 <- mean( (pre.b6$means - z.pred)^2 )
  mspe.b7 <- mean( (pre.b7$means - z.pred)^2 )
  mspe.cc <- mean( (pre.cc$means - z.pred)^2 )

  mls.b2  <- mean(correlationVecchia::logscore(pre.b2$means, pre.b2$vars, z.pred))
  mls.b4  <- mean(correlationVecchia::logscore(pre.b4$means, pre.b4$vars, z.pred))
  mls.b5  <- mean(correlationVecchia::logscore(pre.b5$means, pre.b5$vars, z.pred))
  mls.b6  <- mean(correlationVecchia::logscore(pre.b6$means, pre.b6$vars, z.pred))
  mls.b7  <- mean(correlationVecchia::logscore(pre.b7$means, pre.b7$vars, z.pred))
  mls.cc  <- mean(correlationVecchia::logscore(pre.cc$means, pre.cc$vars, z.pred))

  return(list(mspe = c(mspe.b2, mspe.b4, mspe.b5, mspe.b6, mspe.b7, mspe.cc), mls = c(mls.b2, mls.b4, mls.b5, mls.b6, mls.b7, mls.cc)))
}

ms <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
output <- list()
for(i in 1:length(ms)) output[[i]] <- comparison(ms[i])

#########################################################################

output.m10 <- output[[1]]

ftn <- function(x) log10(-x)

par(mar=c(7.1, 5.1, 4.1, 2.1), mfrow = c(1, 2))
barplot(sqrt(c(output.m10$mspe[1], output.m10$mspe[2], output.m10$mspe[3], output.m10$mspe[4], output.m10$mspe[5], output.m10$mspe[6])), ylim = sqrt(range(c(output.m10$mspe[1], output.m10$mspe[2], output.m10$mspe[3], output.m10$mspe[4], output.m10$mspe[5], output.m10$mspe[6]))), names.arg = c("S-E-MM + J-E-NN", "S-E-MM + C-NN", "T-ord + T-NN", "T-ord + J-E-NN", "T-ord + C-NN", "C-MM + C-NN"), horiz = FALSE, las = 2, cex.names = 0.75, col = c(1, 2, 3, 4, 5, 6), main = "RMSPE")
barplot(ftn(c(output.m10$mls[1], output.m10$mls[2], output.m10$mls[3], output.m10$mls[4], output.m10$mls[5], output.m10$mls[6])), ylim = range(ftn(c(output.m10$mls[1], output.m10$mls[2], output.m10$mls[3], output.m10$mls[4], output.m10$mls[5], output.m10$mls[6]))), names.arg = c("S-E-MM + J-E-NN", "S-E-MM + C-NN", "T-ord + T-NN", "T-ord + J-E-NN", "T-ord + C-NN", "C-MM + C-NN"), horiz = FALSE, las = 2, cex.names = 0.75, col = c(1, 2, 3, 4, 5, 6), main = "Log10-scale logscore")

#########################################################################

mspe.b2 <- c()
mspe.b4 <- c()
mspe.b5 <- c()
mspe.b6 <- c()
mspe.b7 <- c()
mspe.cc <- c()
mls.b2 <- c()
mls.b4 <- c()
mls.b5 <- c()
mls.b6 <- c()
mls.b7 <- c()
mls.cc <- c()
for(i in 1:length(ms)) {

  mspe.b2[i] <- output[[i]]$mspe[1]
  mspe.b4[i] <- output[[i]]$mspe[2]
  mspe.b5[i] <- output[[i]]$mspe[3]
  mspe.b6[i] <- output[[i]]$mspe[4]
  mspe.b7[i] <- output[[i]]$mspe[5]
  mspe.cc[i] <- output[[i]]$mspe[6]

  mls.b2[i] <- output[[i]]$mls[1]
  mls.b4[i] <- output[[i]]$mls[2]
  mls.b5[i] <- output[[i]]$mls[3]
  mls.b6[i] <- output[[i]]$mls[4]
  mls.b7[i] <- output[[i]]$mls[5]
  mls.cc[i] <- output[[i]]$mls[6]
}

#########################################################################

par(mfrow = c(1, 2))

plot(ms, sqrt(mspe.b2), ylim = sqrt(range(c(mspe.b2, mspe.b5, mspe.cc))), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "RMSPE", main = "Root Mean squared prediction error")
points(ms, sqrt(mspe.b5), type = "o", col = c("#FF7F00"), pch = 17)
points(ms, sqrt(mspe.cc), type = "o", col = c("#E41A1C"), pch = 18)

ftn <- function(x) log10(-x)

plot(ms, ftn(mls.b2), ylim = ftn(range(c(mls.b2, mls.b5, mls.cc))), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "log10(logscore)", main = "Log10-scale Logarhithmic Score")
points(ms, ftn(mls.b5), type = "o", col = c("#FF7F00"), pch = 17)
points(ms, ftn(mls.cc), type = "o", col = c("#E41A1C"), pch = 18)

par(mfrow = c(1, 3))

plot(ms, sqrt(mspe.b2), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "RMSPE", main = "S-E-MM + J-E-NN")
plot(ms, sqrt(mspe.b5), type = "o", col = c("#FF7F00"), pch = 17, xlab = "increasing m", ylab = "RMSPE", main = "S-E-MM + C-NN")
plot(ms, sqrt(mspe.cc), type = "o", col = c("#E41A1C"), pch = 18, xlab = "increasing m", ylab = "RMSPE", main = "C-MM + C-NN")

plot(ms, ftn(mls.b2), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "log10(logscore)", main = "S-E-MM + J-E-NN")
plot(ms, ftn(mls.b5), type = "o", col = c("#FF7F00"), pch = 17, xlab = "increasing m", ylab = "log10(logscore)", main = "S-E-MM + C-NN")
plot(ms, ftn(mls.cc), type = "o", col = c("#E41A1C"), pch = 18, xlab = "increasing m", ylab = "log10(logscore)", main = "C-MM + C-NN")

par(mfrow = c(1, 1))

#########################################################################

save(locs, z, X.train, locs.pred, z.pred, X.pred, ns_obs, ns_pred,
     ms, output,
     mspe.b2, mspe.b4, mspe.b5, mspe.b6, mspe.b7, mspe.cc,
     mls.b2, mls.b4, mls.b5, mls.b6, mls.b7, mls.cc,
     file = "DATA/joint_CRCM_NCEP_10152021_prediction.RData")

