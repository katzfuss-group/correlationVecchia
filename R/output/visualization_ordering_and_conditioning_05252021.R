library(correlationVecchia)
library(dplyr)
library(ggplot2)

rm(list = ls())

# set.seed(20200524)
set.seed(05252021)

n         <- 20^2
sqrtn     <- as.integer(sqrt(n))
m         <- 6
visn      <- 40
pivot     <- 38
covparms  <- c(1, 0.1, 4)

# n         <- 20^2
# sqrtn     <- as.integer(sqrt(n))
# m         <- 6
# visn      <- 25
# pivot     <- 19
# covparms  <- c(1, 0.1, 4)

locs      <- as.matrix(expand.grid(seq(from = 0.1, to = 0.9, length.out = sqrtn), seq(from = 0.1, to = 0.9, length.out = sqrtn)))
noise     <- matrix(0.05 * runif(2 * sqrtn^2, min = 0, max = 1), nrow = sqrtn^2, ncol = 2)
noise     <- noise - mean(noise)
locs      <- locs + noise
n         <- nrow(locs)

plot(locs[, 1], locs[, 2], pch = 18, col = "gray80", xlab = "x", ylab = "y", main = "Euclidean MN on Original domain")

cormat    <- cov_expo_aniso(locs = locs, covparms = covparms)
ord       <- order_maxmin_euclidean(locs = locs, initial.pt = NULL)

locsord   <- locs[ord, , drop = FALSE]
cormat    <- cormat[ord, ord]

cond.sets <- conditioning_nn_Rcpp(m = m, d = 1 - cormat) + 1

text(x = locsord[seq(visn), 1], y = locsord[seq(visn), 2], labels = seq(visn))

idx <- cond.sets[pivot, ] ;

points(locsord[idx[-1], 1], locsord[idx[-1], 2], pch = 1, cex = 4, lwd = 3, col = "green")
points(locsord[idx[1], 1], locsord[idx[1], 2], pch = 0, cex = 4, lwd = 3, col = "red")

plot(locs[, 1], locs[, 2], pch = 18, col = "gray80", xlab = "x", ylab = "y", main = "Correlation-based MN on Original domain")

cormat    <- cov_expo_aniso(locs = locs, covparms = covparms)
ord       <- order_maxmin_correlation_inverseDist(locs = locs, d.inv = cormat, initial.pt = ord[1])

locsord   <- locs[ord, , drop = FALSE]
cormat    <- cormat[ord, ord]

cond.sets <- conditioning_nn_Rcpp(m = m, d = 1 - cormat) + 1

text(x = locsord[seq(visn), 1], y = locsord[seq(visn), 2], labels = seq(visn))

idx <- cond.sets[pivot, ] ;

points(locsord[idx[-1], 1], locsord[idx[-1], 2], pch = 1, cex = 4, lwd = 3, col = "blue")
points(locsord[idx[1], 1], locsord[idx[1], 2], pch = 0, cex = 4, lwd = 3, col = "red")

locsord.new <- locsord
locsord.new[, 1] <- locsord[, 1] * covparms[3]

sum((cov_expo_aniso(locsord, covparms) - cov_expo_iso(locsord.new, covparms[1:2]))^2) # It must be very close to 0!

plot(locsord.new[, 1], locsord.new[, 2], pch = 18, col = "gray80", xlab = "x", ylab = "y", main = "Correlation-based MN on Transformed domain")

text(x = locsord.new[seq(visn), 1], y = locsord.new[seq(visn), 2], labels = seq(visn))
points(locsord.new[idx[-1], 1], locsord.new[idx[-1], 2], pch = 1, cex = 4, lwd = 3, col = "blue")
points(locsord.new[idx[1], 1], locsord.new[idx[1], 2], pch = 0, cex = 4, lwd = 3, col = "red")

