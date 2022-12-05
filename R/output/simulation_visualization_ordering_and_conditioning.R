####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script generates Figure 1 (visualization_MN.pdf).
###
####################################################################################

library(correlationVecchia)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

rm(list = ls())

set.seed(05252021)

### setting ####################################################################################################################

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

### Locs #######################################################################################################################

locs      <- as.matrix(expand.grid(seq(from = 0.1, to = 0.9, length.out = sqrtn), seq(from = 0.1, to = 0.9, length.out = sqrtn)))
noise     <- matrix(0.05 * runif(2 * sqrtn^2, min = 0, max = 1), nrow = sqrtn^2, ncol = 2)
noise     <- noise - mean(noise)
locs      <- locs + noise
n         <- nrow(locs) # Redefine n!

locs.new  <- locs
locs.new[1, ] <- locs.new[1, ] * covparms[3]

### Function ###################################################################################################################

visout_MN <- function(locs, covparms, cormat, m, method)
{
  if( method %in% c("e", "euc", "euclidean", "E-MN", "e-MN") ) {

    ord <- order_maxmin_euclidean(locs = locs, initial.pt = NULL)

    locsord   <- locs[ord, , drop = FALSE]
    cormat    <- cormat[ord, ord]

    ord       <- order_maxmin_euclidean(locs = locs, initial.pt = NULL)

    locsord   <- locs[ord, , drop = FALSE]

    cond.sets <- GpGp::find_ordered_nn_brute(locsord, m)

  } else if( method %in% c("c", "cor", "correlation", "C-MN", "c-MN") ) {

    idx       <- order_maxmin_euclidean(locs = locs, initial.pt = NULL)[1]
    ord       <- order_maxmin_correlation(locs = locs, d.inv = cormat, initial.pt = idx)

    locsord   <- locs[ord, , drop = FALSE]
    cormat    <- cormat[ord, ord]

    cond.sets <- conditioning_m_Rcpp(m = m, d = 1 - cormat) + 1

  } else {

    stop("Check the argument method!")
  }

  locsord.new <- locsord
  locsord.new[, 1] <- locsord.new[, 1] * covparms[3]

  return(list(locsord = locsord, locsord.scaled = locsord.new, ord = ord, cond.sets = cond.sets))
}

### E-MN #######################################################################################################################

cormat      <- cov_expo_iso(locs = locs, covparms = covparms[1:2])
visout.euc  <- visout_MN(locs = locs, covparms, cormat = cormat, m = m, method = "E-MN")

### C-MN #######################################################################################################################

cormat      <- cov_expo_aniso(locs = locs, covparms = covparms)
visout.cor  <- visout_MN(locs = locs, covparms, cormat = cormat, m = m, method = "C-MN")

### Visualization ##############################################################################################################

plot(locs[, 1], locs[, 2], pch = 18, col = "gray80", xlab = "time", ylab = "space", main = "E-MN on Original domain")

text(x = visout.euc$locsord[seq(visn), 1], y = visout.euc$locsord[seq(visn), 2], labels = seq(visn))

idx <- visout.euc$cond.sets[pivot, ] ;

points(visout.euc$locsord[idx[-1], 1], visout.euc$locsord[idx[-1], 2], pch = 1, cex = 4, lwd = 3, col = "green")
points(visout.euc$locsord[idx[1], 1], visout.euc$locsord[idx[1], 2], pch = 0, cex = 4, lwd = 3, col = "red")

locsord.new <- visout.euc$locsord.scaled

sum((cov_expo_aniso(visout.euc$locsord, covparms) - cov_expo_iso(locsord.new, covparms[1:2]))^2) # It must be very close to 0!

plot(locsord.new[, 1], locsord.new[, 2], pch = 18, col = "gray80", xlab = "time", ylab = "space", main = "E-MN on Transformed domain")

text(x = locsord.new[seq(visn), 1], y = locsord.new[seq(visn), 2], labels = seq(visn))
points(locsord.new[idx[-1], 1], locsord.new[idx[-1], 2], pch = 1, cex = 4, lwd = 3, col = "green")
points(locsord.new[idx[1], 1], locsord.new[idx[1], 2], pch = 0, cex = 4, lwd = 3, col = "red")

### Visualization ##############################################################################################################

plot(locs[, 1], locs[, 2], pch = 18, col = "gray80", xlab = "time", ylab = "space", main = "C-MN on Original domain")

text(x = visout.cor$locsord[seq(visn), 1], y = visout.cor$locsord[seq(visn), 2], labels = seq(visn))

idx <- visout.cor$cond.sets[pivot, ] ;

points(visout.cor$locsord[idx[-1], 1], visout.cor$locsord[idx[-1], 2], pch = 1, cex = 4, lwd = 3, col = "blue")
points(visout.cor$locsord[idx[1], 1], visout.cor$locsord[idx[1], 2], pch = 0, cex = 4, lwd = 3, col = "red")

locsord.new <- visout.cor$locsord.scaled

sum((cov_expo_aniso(visout.cor$locsord, covparms) - cov_expo_iso(locsord.new, covparms[1:2]))^2) # It must be very close to 0!

plot(locsord.new[, 1], locsord.new[, 2], pch = 18, col = "gray80", xlab = "time", ylab = "space", main = "C-MN on Transformed domain")

text(x = locsord.new[seq(visn), 1], y = locsord.new[seq(visn), 2], labels = seq(visn))
points(locsord.new[idx[-1], 1], locsord.new[idx[-1], 2], pch = 1, cex = 4, lwd = 3, col = "blue")
points(locsord.new[idx[1], 1], locsord.new[idx[1], 2], pch = 0, cex = 4, lwd = 3, col = "red")

### ggplot2 ####################################################################################################################



idx             <- visout.euc$cond.sets[pivot, ]
grp             <- rep(1, n)
grp[seq(visn)]  <- 2
grp[idx]        <- 3
grp[idx[1]]     <- 4

visdat1  <- data.frame(x = visout.euc$locsord[, 1], y = visout.euc$locsord[, 2])
visdat2   <- visdat1[grp == 2 | grp == 3 | grp == 4, ] ; visdat2$ID <- seq(visn)
visdat3   <- visdat1[grp == 3, ]
visdat4   <- visdat1[grp == 4, ]

eps <- 0.07

##

p0 <- ggplot() + geom_point(data = visdat1, aes(x, y), colour = "gray90") +
  geom_point(data = visdat2, aes(x, y), colour = "gray90") + # geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  geom_point(data = visdat3, aes(x, y), colour = "#377EB8", shape = 0, size = 8, stroke = 2) +
  geom_point(data = visdat4, aes(x, y), colour = "#E41A1C", shape = 1, size = 10, stroke = 2) + geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  theme_bw() + xlim(0+eps, 1-eps) + ylim(0+eps, 1-eps) + xlab("") + ylab("") + theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14))
p0

ggplot2::ggsave("visualization_E-MN.pdf", p0, width = 5.7, height = 5.7)

##

p1 <- ggplot() + geom_point(data = visdat1, aes(x, y), colour = "gray90") +
  geom_point(data = visdat2, aes(x, y), colour = "gray90") + # geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  geom_point(data = visdat3, aes(x, y), colour = "#377EB8", shape = 0, size = 8, stroke = 2) +
  geom_point(data = visdat4, aes(x, y), colour = "#E41A1C", shape = 1, size = 10, stroke = 2) + geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  theme_bw() + xlim(0+eps, 1-eps) + ylim(0+eps, 1-eps) + xlab("time") + ylab("space") + theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14))
p1

idx             <- visout.euc$cond.sets[pivot, ]
grp             <- rep(1, n)
grp[seq(visn)]  <- 2
grp[idx]        <- 3
grp[idx[1]]     <- 4

visdat1  <- data.frame(x = visout.euc$locsord.scaled[, 1], y = visout.euc$locsord.scaled[, 2])
visdat2   <- visdat1[grp == 2 | grp == 3 | grp == 4, ] ; visdat2$ID <- seq(visn)
visdat3   <- visdat1[grp == 3, ]
visdat4   <- visdat1[grp == 4, ]

eps <- 0.075

p2 <- ggplot() + geom_point(data = visdat1, aes(x, y), colour = "gray90") +
  geom_point(data = visdat2, aes(x, y), colour = "gray90") + # geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  geom_point(data = visdat3, aes(x, y), colour = "#377EB8", shape = 0, size = 8, stroke = 2) +
  geom_point(data = visdat4, aes(x, y), colour = "#E41A1C", shape = 1, size = 10, stroke = 2) + geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  theme_bw() + xlim(covparms[3]*(0+eps), covparms[3]*(1-eps)) + ylim(0+eps, 1-eps) + xlab("scaled time") + ylab("scaled space") + theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14))
p2

idx             <- visout.cor$cond.sets[pivot, ]
grp             <- rep(1, n)
grp[seq(visn)]  <- 2
grp[idx]        <- 3
grp[idx[1]]     <- 4

visdat1  <- data.frame(x = visout.cor$locsord[, 1], y = visout.cor$locsord[, 2])
visdat2   <- visdat1[grp == 2 | grp == 3 | grp == 4, ] ; visdat2$ID <- seq(visn)
visdat3   <- visdat1[grp == 3, ]
visdat4   <- visdat1[grp == 4, ]

eps <- 0.07

p3 <- ggplot() + geom_point(data = visdat1, aes(x, y), colour = "gray90") +
  geom_point(data = visdat2, aes(x, y), colour = "gray90") + # geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  geom_point(data = visdat3, aes(x, y), colour = "#4DAF4A", shape = 0, size = 8, stroke = 2) +
  geom_point(data = visdat4, aes(x, y), colour = "#E41A1C", shape = 1, size = 10, stroke = 2) + geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  theme_bw() + xlim(0+eps, 1-eps) + ylim(0+eps, 1-eps) + xlab("time") + ylab("space") + theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14))
p3

idx             <- visout.cor$cond.sets[pivot, ]
grp             <- rep(1, n)
grp[seq(visn)]  <- 2
grp[idx]        <- 3
grp[idx[1]]     <- 4

visdat1  <- data.frame(x = visout.cor$locsord.scaled[, 1], y = visout.cor$locsord.scaled[, 2])
visdat2   <- visdat1[grp == 2 | grp == 3 | grp == 4, ] ; visdat2$ID <- seq(visn)
visdat3   <- visdat1[grp == 3, ]
visdat4   <- visdat1[grp == 4, ]

eps <- 0.07

p4 <- ggplot() + geom_point(data = visdat1, aes(x, y), colour = "gray90") +
  geom_point(data = visdat2, aes(x, y), colour = "gray90") + # geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  geom_point(data = visdat3, aes(x, y), colour = "#4DAF4A", shape = 0, size = 8, stroke = 2) +
  geom_point(data = visdat4, aes(x, y), colour = "#E41A1C", shape = 1, size = 10, stroke = 2) + geom_text(data = visdat2, aes(x, y), label = seq(visn), size = 5) +
  theme_bw() + xlim(covparms[3]*(0+eps), covparms[3]*(1-eps)) + ylim(0+eps, 1-eps) + xlab("scaled time") + ylab("scaled space") + theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14))
p4

p.all <- grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2, widths = c(1, covparms[3]))

ggplot2::ggsave("visualization_MN.pdf", p.all, width = 15.2, height = 5.7)


