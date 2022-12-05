####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to compare performance of approximations when there are no inputs.
###
####################################################################################

rm(list = ls())
set.seed(09212021)

##

gc()

memory.size() # memory.limit()

##

library(correlationVecchia)

###################################################################

output                <- generate_hnm(mu0 = 0, sig20 = 1, k = 12, mlt = 2)

length(output$y)

covmat                <- covmat_hnm(output)
cormat                <- cormat_hnm(output)

# fields::imagePlot(cormat)

###################################################################

initial.pt            <- NULL
ms                    <- c(3, 5, 10, 15, 20, 25, 30, 35, 40, 45)

###################################################################

approx <- list()
for(i in 1:length(ms)) {

  approx[[i]]           <- cvecchia_specify_noneuc(cormat = cormat, initial.pt = initial.pt, m = ms[i])
}

n.approx              <- length(approx)
kls1                  <- rep(NA, n.approx)
for(i in 1:n.approx) {

  covmat.ord            <- covmat[approx[[i]]$ord, approx[[i]]$ord]
  U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord)$U

  revord                <- order(approx[[i]]$ord)
  covmat.hat            <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]

  kls1[i]               <- kldiv(covmat, covmat.hat, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps)
}

###################################################################

approx <- list()
for(i in 1:length(ms)) {

  approx[[i]]           <- lvecchia_specify_noneuc(cormat = cormat, m = ms[i])
}

n.approx              <- length(approx)
kls2                  <- rep(NA, n.approx)
for(i in 1:n.approx) {

  covmat.ord            <- covmat[approx[[i]]$ord, approx[[i]]$ord]
  U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord)$U

  revord                <- order(approx[[i]]$ord)
  covmat.hat            <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]

  kls2[i]               <- kldiv(covmat, covmat.hat, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps)
}

###################################################################

yrange <- log10(range(c(kls1, kls2)))

plot(ms, log10(kls1), ylim = yrange, type = "o", col = "#E41A1C", pch = 16, lwd = 2, cex = 1.25, xlab = "increasing m", ylab = "log10(KL)", main = "Non-spatial example: Hierarchical normal model")
lines(ms, log10(kls2), type = "o", col = "#4DAF4A", pch = 17, lwd = 2, cex = 1.25, xlab = "increasing m", ylab = "log10(KL)", main = "Non-spatial example: Hierarchical normal model")

axis(2)
axis(side = 1, at = ms, labels = ms)

abline(v = ms, lty = 2, col = "black")

###################################################################

nsim <- 200

kls3 <- list()
for(j in 1:nsim) {

  # approx <- list()
  # for(i in 1:length(ms)) {
  #
  #   approx[[i]]           <- rvecchia_specify_nonE(cormat = cormat, m = ms[i])
  # }

  approx                <- rvecchias_specify_noneuc(cormat = cormat, ms = ms)

  n.approx              <- length(approx)
  kls                   <- rep(NA, n.approx)
  for(i in 1:n.approx) {

    covmat.ord            <- covmat[approx[[i]]$ord, approx[[i]]$ord]
    U                     <- GPvecchia::createU(approx[[i]], c(1), 0, covmodel = covmat.ord)$U

    revord                <- order(approx[[i]]$ord)
    covmat.hat            <- as.matrix(solve(Matrix::tcrossprod(U)))[revord,revord]

    kls[i]                <- kldiv(covmat, covmat.hat, method.modify = NULL, pivot = FALSE, tol = .Machine$double.eps)
  }

  kls3[[j]] <- kls
}

###################################################################

# png(file = "nonspatial example.png", width = 700, height = 400)

yrange <- log10(range(c(kls1, kls2, unlist(kls3))))

plot(ms, log10(kls1), ylim = yrange, type = "o", col = "#E41A1C", pch = 16, lwd = 2, cex = 1.25, xlab = "increasing m", ylab = "log10(KL)", main = "Non-spatial example: Hierarchical normal model")
lines(ms, log10(kls2), type = "o", col = "#4DAF4A", pch = 17, lwd = 2, cex = 1.25)

for(j in 1:nsim) lines(ms, log10(kls3[[j]]), type = "o", col = scales::alpha("#377EB8", 0.05), pch = 15, lwd = 2, cex = 1)

axis(2)
axis(side = 1, at = ms, labels = ms)

abline(v = ms, lty = 2, col = "black")

legend("bottomleft", legend = c("RVecchia", "LVecchia", "CVecchia"), lty = 1, lwd = 2, col = c("#377EB8", "#4DAF4A", "#E41A1C"))

# dev.off()

###################################################################

save(output, cormat, ms, kls1, kls2, kls3, file = "simout_noneuc_10012021.RData")

