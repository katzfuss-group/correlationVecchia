####################################################################
####
####  Author: Myeongjong (MJ) Kang (kmj.stat@gmail.com)
####
####  Description: This example shows that the GpGp::find_ordered_nn() function has an error. But we need to find reasons.
####
####################################################################

set.seed(29913)

source("2_corrvecchia/corrvecchia.R")

n             <- 15^2
m             <- 10
locs          <- matrix(runif(n * 2, 0, 1), n, 2)

nn1 <- GpGp::find_ordered_nn(locs, m)
nn2 <- GpGp::find_ordered_nn_brute(locs, m)

D <- as.matrix(dist(locs, method = "euclidean"))
nn3 <- conditioning_nn(m, D)

head(nn1 == nn2, 50)
head(nn2 == nn3, 50)

nn1[26, ]
nn2[26, ]
nn3[26, ]

D[26, 18]
D[26, 21]

library(microbenchmark)
microbenchmark(GpGp::find_ordered_nn(locs, m), GpGp::find_ordered_nn_brute(locs, m), conditioning_nn(m, D), times = 10)
