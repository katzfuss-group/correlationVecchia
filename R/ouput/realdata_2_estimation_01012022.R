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
library(correlationVecchia)

#########################################################################

load("DATA/joint_CRCM_NCEP_10112021.RData")

df.joint %>% glimpse()

#########################################################################

# remove some locations (spacewise)

# splt.swise <- splitData(df.joint = df.joint, method = "spacewise")

#########################################################################

# remove some time points (timewise)

# splt.twise <- splitData(df.joint = df.joint, method = "timewise")

#########################################################################

# remove some spacetime blocks (blockwise)

splt.bwise <- splitData(df.joint = df.joint, method = "blockwise", size.blck = 5^2)

#########################################################################

train.joint     <- splt.bwise$df.train
test.joint      <- splt.bwise$df.test

#########################################################################

inputs      <- train.joint[, c("x1", "x2", "t", "d")]
inputs.lst  <- split(inputs[, , drop = FALSE], inputs$d)

z           <- train.joint[, c("z")]
mean.z      <- mean(z)
sd.z        <- sd(z)

# z           <- (z - mean.z) / sd.z
X           <- as.matrix(data.frame(rep(1, length(z)), inputs$d))

m           <- 50

##

inputs      <- as.matrix(inputs)
inputs.lst  <- lapply(inputs.lst, as.matrix)

## (~ 6 hrs)

Sys.time()

fit.b1 <- fit_scaled_bs_mulv(approx = 1, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), find.vcf = FALSE)

Sys.time()

fit.b2 <- fit_scaled_bs_mulv(approx = 2, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), find.vcf = FALSE)

Sys.time()

fit.b3 <- fit_scaled_bs_mulv(approx = 3, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), find.vcf = FALSE)

Sys.time()

fit.b4 <- fit_scaled_bs_mulv(approx = 4, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), find.vcf = FALSE)

Sys.time()

fit.b5 <- fit_scaled_bs_mulv(approx = 5, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), find.vcf = FALSE)

Sys.time()

fit.b6 <- fit_scaled_bs_mulv(approx = 6, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), find.vcf = FALSE)

Sys.time()

fit.b7 <- fit_scaled_bs_mulv(approx = 7, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), find.vcf = FALSE)

Sys.time()

fit.cc <- fit_scaled_cv_mulv(            y = z, inputs = inputs,     X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), find.vcf = FALSE)

Sys.time()

##

round(fit.b1$covparms, 4)
round(fit.b2$covparms, 4)
round(fit.b3$covparms, 4)
round(fit.b4$covparms, 4)
round(fit.b5$covparms, 4)
round(fit.b6$covparms, 4)
round(fit.b7$covparms, 4)
round(fit.cc$covparms, 4)

#########################################################################

save(df.joint, train.joint, test.joint, mean.z, sd.z, m, fit.b1, fit.b2, fit.b3, fit.b4, fit.b5, fit.b6, fit.b7, fit.cc, file = "DATA/joint_CRCM_NCEP_10152021_estimates.RData")

