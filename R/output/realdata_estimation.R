####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script is to perform Fisher scoring for real data analysis.
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

# m           <- 20 # 50

##

inputs      <- as.matrix(inputs)
inputs.lst  <- lapply(inputs.lst, as.matrix)

#########################################################################

fit_all <- function(m)
{
  # Sys.time()
  #
  # ptm.b1 <- proc.time()
  # fit.b1 <- fit_scaled_bs_mulv(approx = 1, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), group = TRUE, find.vcf = FALSE, print.level = 2)
  # ptm.b1 <- proc.time() - ptm.b1

  # Sys.time()

  ptm.b2 <- proc.time()
  fit.b2 <- fit_scaled_bs_mulv(approx = 2, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), group = TRUE, find.vcf = FALSE, print.level = 2)
  ptm.b2 <- proc.time() - ptm.b2

  # Sys.time()

  ptm.b3 <- proc.time()
  fit.b3 <- fit_scaled_bs_mulv(approx = 3, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), group = TRUE, find.vcf = FALSE, print.level = 2)
  ptm.b3 <- proc.time() - ptm.b3

  # Sys.time()
  #
  # ptm.b4 <- proc.time()
  # fit.b4 <- fit_scaled_bs_mulv(approx = 4, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), group = TRUE, find.vcf = FALSE, print.level = 2)
  # ptm.b4 <- proc.time() - ptm.b4

  # Sys.time()

  ptm.b5 <- proc.time()
  fit.b5 <- fit_scaled_bs_mulv(approx = 5, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), group = FALSE, find.vcf = FALSE, print.level = 2)
  ptm.b5 <- proc.time() - ptm.b5

  # Sys.time()
  #
  # ptm.b6 <- proc.time()
  # fit.b6 <- fit_scaled_bs_mulv(approx = 6, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), group = TRUE, find.vcf = FALSE, print.level = 2)
  # ptm.b6 <- proc.time() - ptm.b6

  # Sys.time()

  ptm.b7 <- proc.time()
  fit.b7 <- fit_scaled_bs_mulv(approx = 7, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), group = TRUE, find.vcf = FALSE, print.level = 2)
  ptm.b7 <- proc.time() - ptm.b7

  # Sys.time()

  ptm.b8 <- proc.time()
  fit.b8 <- fit_scaled_bs_mulv(approx = 8, y = z, inputs = inputs.lst, X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), group = TRUE, find.vcf = FALSE, print.level = 2)
  ptm.b8 <- proc.time() - ptm.b8

  # Sys.time()

  ptm.cc <- proc.time()
  fit.cc <- fit_scaled_cv_mulv(            y = z, inputs = inputs,     X = X, ms = c(m), trend = "X", var.ini = sd.z^2, ranges.ini = c(0.6, 0.5, 0.1, 2.0), nu = 0.75, nug = 1e-4, n.est = nrow(inputs), group = TRUE, find.vcf = FALSE, print.level = 2)
  ptm.cc <- proc.time() - ptm.cc

  # Sys.time()

  return( list(fit.b2, fit.b3, fit.b5, fit.b7, fit.b8, fit.cc, ptm.b2, ptm.b3, ptm.b5, ptm.b7, ptm.b8, ptm.cc) )
}
#########################################################################

Sys.time() # "2022-11-26 01:02:20 CST"
fit.m5 <- fit_all(m = 5)

Sys.time() # "2022-11-26 01:24:40 CST"
fit.m10 <- fit_all(m = 10)

Sys.time() # "2022-11-26 02:19:51 CST"
fit.m15 <- fit_all(m = 15)

Sys.time()
fit.m20 <- fit_all(m = 20)

Sys.time() # "2022-11-26 07:35:43 CST"
fit.m25 <- fit_all(m = 25)

Sys.time() # "2022-11-26 11:53:28 CST"
fit.m30 <- fit_all(m = 30)

Sys.time() # "2022-11-26 18:29:48 CST"
fit.m35 <- fit_all(m = 35)

Sys.time() # "2022-11-27 07:31:21 CST"
fit.m40 <- fit_all(m = 40)

Sys.time() # "2022-11-27 18:59:01 CST"
fit.m45 <- fit_all(m = 45)

Sys.time() # "2022-11-28 08:02:27 CST"
fit.m50 <- fit_all(m = 50)

Sys.time() # "2022-11-29 02:50:49 CST"

#########################################################################

save(df.joint, train.joint, test.joint, mean.z, sd.z, fit.m5, fit.m10, fit.m15, fit.m20, fit.m25, fit.m30, fit.m35, fit.m40, fit.m45, fit.m50, file = "DATA/joint_CRCM_NCEP_10152021_estimates.RData")
