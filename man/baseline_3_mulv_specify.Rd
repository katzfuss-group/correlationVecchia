% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baseline.R
\name{baseline_3_mulv_specify}
\alias{baseline_3_mulv_specify}
\title{The third baseline approximation for multivariate GPs (S-E-MM + S-E-NN)}
\usage{
baseline_3_mulv_specify(locs, m)
}
\arguments{
\item{locs}{A matrix of locations}

\item{m}{A size of conditioning sets}
}
\value{
vecchia.approx
}
\description{
The third baseline approximation for multivariate GPs (S-E-MM + S-E-NN)
}
\examples{
locs1  <- matrix(runif(20), 10, 2)
locs2  <- matrix(runif(30), 15, 2)
locs   <- list(locs1 = locs1, locs2 = locs2)
m      <- 5

# baseline 1 = separate E-MM + half-and-half E-NN (S-E-MM + HH-E-NN)
# baseline 2 = separate E-MM + joint E-NN (S-E-MM + J-E-NN)
# baseline 3 = separate E-MM + separate E-NN (S-E-MM + S-E-NN)
# baseline 4 = separate E-MM + C-NN (S-E-MM + C-NN)

# Note 1 = For all baseline approximations, the first process is not relevant to the second one
# Note 2 = For the third baseline approximations, two processes are totally irrelevant

### ordering
ord1 <- baseline_1_mulv_specify(locs = locs, m = m)$ord
ord2 <- baseline_2_mulv_specify(locs = locs, m = m)$ord
ord3 <- baseline_3_mulv_specify(locs = locs, m = m)$ord
ord4 <- baseline_4_mulv_specify(locs = locs, m = m,
                                covmodel = cov_flexMatern_biv,
                                sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2),
                                nu.mat = matrix(0.5, 2, 2),
                                alpha.mat = matrix(1, 2, 2))$ord

all.equal(ord1, ord2) ; all.equal(ord1, ord3) ; all.equal(ord1, ord4)

### conditioning
ind <- seq(from = m, to = 1, by = -1)
cond1 <- baseline_1_mulv_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
cond2 <- baseline_2_mulv_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
cond3 <- baseline_3_mulv_specify(locs = locs, m = m)$U.prep$revNNarray[, ind]
cond4 <- baseline_4_mulv_specify(locs = locs, m = m,
                                 covmodel = cov_flexMatern_biv,
                                 sigma.mat = matrix(c(1, 0.5, 0.5, 1), 2, 2),
                                 nu.mat = matrix(0.5, 2, 2),
                                 alpha.mat = matrix(1, 2, 2))$U.prep$revNNarray[, ind]

cond1 == cond2 # similar
cond1 == cond3 # different
cond3 == cond4 # similar
}
