####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
####################################################################################

rm(list = ls())

library(correlationVecchia)
library(microbenchmark)

set.seed(08122020)

##### functions for only ordering + conditioning #####

ordcond_packages <- function(locs, m)
{
  ord         <- GPvecchia::order_maxmin_exact(locs = locs)
  
  cond.sets   <- GpGp::find_ordered_nn(locs = locs[ord, ], m = m)
  
  return(list(ord = ord, cond = cond.sets))
}

ordcond_bruteforce <- function(locs, m, initial.pt)
{
  ord         <- order_maxmin_euclidean(locs = locs, initial.pt = initial.pt)
  
  cond.sets   <- conditioning_nn(m = m, d = fields::rdist(locs[ord, ]))
  
  return(list(ord = ord, cond = cond.sets))
}

ordcond_florian <- function(locs, rho, initial.pt)
{
  output      <- rho_based_sortSparse(locs = locs, rho = rho, initial.pt = initial.pt, distype = 'euclidean', covmodel = 'cov_expo_iso', covparms = c(1, 0.1))
  
  output$sparsity   <- as(output$sparsity, "TsparseMatrix")
  maxsize           <- max(table(output$sparsity@j))
  output$sparsity   <- cbind(rev(output$sparsity@j), rev(output$sparsity@i))
  
  cond.sets   <- conditioning_Rcpp(output$sparsity[, 1], output$sparsity[, 2], output$P-1, maxsize, locs, 'euclidean', 'cov_expo_iso', c(1, 0.1))
  cond.sets   <- cond.sets + 1
  
  return(list(ord = output$P, cond = cond.sets))
}

##### small n #####

n       <- 100
locs    <- matrix(runif(n * 2), n, 2)

output_packages     <- ordcond_packages(locs = locs, m = 5)
output_bruteforce   <- ordcond_bruteforce(locs, m = 5, initial.pt = output_packages$ord[1])
output_florian      <- ordcond_florian(locs = locs, rho = 2.65, initial.pt = output_packages$ord[1])

all.equal(output_packages$ord, output_bruteforce$ord)
all.equal(output_packages$ord, output_florian$ord)

sum(!is.na(output_packages$cond))
sum(!is.na(output_bruteforce$cond))
sum(!is.na(output_florian$cond))

output_packages$cond[seq(10), ]
output_bruteforce$cond[seq(10), ]
output_florian$cond[seq(10), ]

microbenchmark::microbenchmark(ordcond_packages(locs = locs, m = 5),
                               ordcond_bruteforce(locs, m = 5, initial.pt = output_packages$ord[1]),
                               ordcond_florian(locs = locs, rho = 2.65, initial.pt = output_packages$ord[1]),
                               times = 10)

##### large n #####

n       <- 10000
locs    <- matrix(runif(n * 2), n, 2)

output_packages     <- ordcond_packages(locs = locs, m = 100)
output_florian      <- ordcond_florian(locs = locs, rho = 11.25, initial.pt = output_packages$ord[1])

all.equal(output_packages$ord, output_florian$ord)

sum(!is.na(output_packages$cond))
sum(!is.na(output_florian$cond))

microbenchmark::microbenchmark(ordcond_packages(locs = locs, m = 100),
                               ordcond_florian(locs = locs, rho = 11.25, initial.pt = output_packages$ord[1]),
                               times = 5)

####################################################################################
###   old code
####################################################################################

# n <- 20^2
# d <- 2
# locs <- matrix(runif(n * d), nrow = n, ncol = d)
# 
# rho <- 1.5
# m <- 5
# 
# initial.pt <- 1
# covparms <- c(1, 0.1)
# 
# # to compare it ordcond with sortSparse_Rcpp(locs, rho, initial.pt, covmodel, covparms)
# corr_ordcond_timeCheck <- function(locs, m, initial.pt = NULL, covmodel, covparms) 
# {
#   locs  <- as.matrix(locs)
#   p     <- ncol(locs)
#   n     <- nrow(locs)
#   
#   rho   <- correlationVecchia:::.correlation(locs = locs, covmodel = covmodel, covparms = covparms, abs.corr = FALSE)
#   ord   <- order_maxmin_correlation_inverseDist(locs = locs, d.inv = rho, initial.pt = initial.pt)
#   
#   locsord     <- locs[ord, , drop = FALSE]
#   
#   if(is.matrix(covmodel)) covmodel <- covmodel[ord, ord]
#   
#   rho         <- rho[ord, ord]
#   cond.sets   <- conditioning_nn(m = m, d = 1 - rho)
#   
#   return(list(ord = ord, cond = cond.sets))
# }
# 
# #######################################################################################################################################################
# 
# # save(locs, n, d, rho, m, initial.pt, covparms, file = "ComparisonLocations.RData")
# # load(file = "ComparisonLocations.RData")
# # 
# # out0 <- GPvecchia::order_maxmin_exact(locs = locs)
# # out1 <- rho_based_sortSparse(locs = locs, rho, initial.pt = out0[1], covmodel = "cov_expo_iso", covparms = covparms)$P
# # out2 <- corrvecchia_specify_knownCovparms(locs = locs, m = m, initial.pt = out0[1], covmodel = cov_expo_iso, covparms = covparms)$ord
# # 
# # out0 == out1
# # out0 == out2
# # 
# # distance <- 100
# # for(i in 1:596) distance <- min(distance, sqrt((sum(locs[out0[i], ]-locs[706, ])^2)))
# # distance
# # 
# # distance <- 100
# # for(i in 1:596) distance <- min(distance, sqrt((sum(locs[out1[i], ]-locs[706, ])^2)))
# # distance
# # 
# # distance <- 100
# # for(i in 1:596) distance <- min(distance, sqrt((sum(locs[out2[i], ]-locs[706, ])^2)))
# # distance
# # 
# # all.equal(rho_based_sortSparse(locs = locs, rho, initial.pt = initial.pt, covmodel = "cov_expo_iso", covparms = covparms)$P,
# #           corrvecchia_specify_knownCovparms(locs = locs, m = m, initial.pt = initial.pt, covmodel = cov_expo_iso, covparms = covparms)$ord)
# 
# microbenchmark(rho_based_sortSparse(locs = locs, rho, initial.pt = initial.pt, covmodel = "cov_expo_iso", covparms = covparms),
#                corrvecchia_specify_knownCovparms(locs = locs, m = m, initial.pt = initial.pt, covmodel = cov_expo_iso, covparms = covparms),
#                times = 1)
# 
# #######################################################################################################################################################
# 
# rhovec <- c(1, 2, 3, 4)
# output.rho_1 <- rho_based_sortSparse(locs = locs, rho = rhovec[1], initial.pt = initial.pt, distype = "euclidean", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
# output.rho_2 <- rho_based_sortSparse(locs = locs, rho = rhovec[2], initial.pt = initial.pt, distype = "euclidean", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
# output.rho_3 <- rho_based_sortSparse(locs = locs, rho = rhovec[3], initial.pt = initial.pt, distype = "euclidean", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
# output.rho_4 <- rho_based_sortSparse(locs = locs, rho = rhovec[4], initial.pt = initial.pt, distype = "euclidean", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
# 
# # rhovec <- c(1.00, 1.20, 1.40, 1.60)
# # output.rho_1 <- rho_based_sortSparse(locs = locs, rho = rhovec[1], initial.pt = initial.pt, distype = "correlation", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
# # output.rho_2 <- rho_based_sortSparse(locs = locs, rho = rhovec[2], initial.pt = initial.pt, distype = "correlation", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
# # output.rho_3 <- rho_based_sortSparse(locs = locs, rho = rhovec[3], initial.pt = initial.pt, distype = "correlation", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
# # output.rho_4 <- rho_based_sortSparse(locs = locs, rho = rhovec[4], initial.pt = initial.pt, distype = "correlation", covmodel = "cov_expo_iso", covparms = c(1, 0.1))
# 
# # output.rho_1 <- rho_based_sortSparse(locs = locs, rho = 1, initial.pt = initial.pt, covmodel = "cov_expo_iso", covparms = c(1, 1))
# # output.rho_2 <- rho_based_sortSparse(locs = locs, rho = 2, initial.pt = initial.pt, covmodel = "cov_expo_iso", covparms = c(1, 1))
# # output.rho_3 <- rho_based_sortSparse(locs = locs, rho = 3, initial.pt = initial.pt, covmodel = "cov_expo_iso", covparms = c(1, 1))
# # output.rho_4 <- rho_based_sortSparse(locs = locs, rho = 4, initial.pt = initial.pt, covmodel = "cov_expo_iso", covparms = c(1, 1))
# 
# # colSums(as.matrix(output.rho_1$sparsity))
# # colSums(as.matrix(output.rho_2$sparsity))
# # colSums(as.matrix(output.rho_3$sparsity))
# # colSums(as.matrix(output.rho_4$sparsity))
# 
# summary(colSums(as.matrix(output.rho_1$sparsity)))
# summary(colSums(as.matrix(output.rho_2$sparsity)))
# summary(colSums(as.matrix(output.rho_3$sparsity)))
# summary(colSums(as.matrix(output.rho_4$sparsity)))
# 
# par(mfrow = c(2, 2))
# plot(colSums(as.matrix(output.rho_1$sparsity)), xlab = "index", ylab = "size of conditioning set", main = paste0("rho = ", rhovec[1]))
# plot(colSums(as.matrix(output.rho_2$sparsity)), xlab = "index", ylab = "size of conditioning set", main = paste0("rho = ", rhovec[2]))
# plot(colSums(as.matrix(output.rho_3$sparsity)), xlab = "index", ylab = "size of conditioning set", main = paste0("rho = ", rhovec[3]))
# plot(colSums(as.matrix(output.rho_4$sparsity)), xlab = "index", ylab = "size of conditioning set", main = paste0("rho = ", rhovec[4]))
# par(mfrow = c(1, 1))
# 
# 
# output.m_5 <- corrvecchia_specify_knownCovparms_2(locs = locs, m = 5, covmodel = cov_expo_iso, covparms = c(1, 0.1))
# 
# rowSums(!is.na(output.m_5$U.prep$revNNarray)) - 1
# 
# #######################################################################################################################################################
# 
# nvec <- c(100, 200, 300, 400, 500, 600, 700, 800, 900)
# d <- 2
# initial.pt <- 1
# covparms <- c(1, 0.1)
# 
# times <- 10
# 
# output.mbased <- list()
# output.rhobased <- list()
# timePerformance <- list()
# 
# for(i in 1:length(nvec)) {
#   
#   n     <- nvec[i]
#   
#   rho   <- 1.5
#   m     <- 5
#   
#   locs                    <- matrix(runif(n * d), nrow = n, ncol = d)
#   output.mbased[[i]]      <- corr_ordcond_timeCheck(locs = locs, m = m, initial.pt = initial.pt, covmodel = cov_expo_iso, covparms = covparms)
#   output.rhobased[[i]]    <- sortSparse_Rcpp(x = locs, rho = rho, initInd = initial.pt - 1, fstr = "cov_expo_iso", covparms = covparms)
#   
#   timePerformance[[i]]    <- microbenchmark(corr_ordcond_timeCheck(locs = locs, m = m, initial.pt = initial.pt, covmodel = cov_expo_iso, covparms = covparms), 
#                                             sortSparse_Rcpp(x = locs, rho = rho, initInd = initial.pt, fstr = "cov_expo_iso", covparms = covparms),
#                                             times = times)
# }
# 
# time.mbased <- c()
# time.rhobased <- c()
# for(i in 1:length(nvec)) {
#   timevec <- sort(timePerformance[[i]]$time)
#   time.rhobased[i] <- mean(timevec[seq(times)])
#   time.mbased[i] <- mean(timevec[seq(times) + times])
# }
# 
# time.rhobased <- time.rhobased / 1e6 # milliseconds
# time.mbased <- time.mbased / 1e6 # milliseconds
# 
# par(mfrow = c(1, 2))
# plot(nvec, time.mbased, type = 'o', col = 'blue', ylab = 'time (millisecond)', xlab = 'sample size', main = 'Time Performance (Ordering + Conditioning)')
# points(nvec, time.rhobased, type = 'o', col = 'red')
# legend("topleft", legend = c("Brute-force Implementation", "Florian's algorithm"), col = c("red", "blue"), lty=1, cex = 0.8)
# 
# plot(nvec, time.rhobased, col = 'red', type = 'o', ylab = 'time (millisecond)', xlab = 'sample size', main = 'Time Performance (Ordering + Conditioning)')
# par(mfrow = c(1, 1))
# 
# #######################################################################################################################################################
# 
# n <- 100^2
# d <- 2
# locs <- matrix(runif(n * d), nrow = n, ncol = d)
# 
# rho <- 1.1
# 
# initial.pt <- 1
# covparms <- c(1, 0.1)
# 
# microbenchmark(fastcorrvecchia_specify_knownCovparms(locs = locs, rho = rho, initial.pt = initial.pt, covmodel = "cov_expo_iso", covparms = covparms), times = 1)
# 
# microbenchmark(rho_based_sortSparse(locs = locs, rho = rho, initial.pt = initial.pt, covmodel = "cov_cpp", covparms = covparms), times = 10)
# 
# output.simple <- rho_based_sortSparse(locs = locs, rho = 3, initial.pt = initial.pt, covmodel = "cov_cpp", covparms = covparms)






