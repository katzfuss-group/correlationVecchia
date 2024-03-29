% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covmodels.R
\name{cov_flexMatern_mulv}
\alias{cov_flexMatern_mulv}
\title{Flexible Matern covariance function for Multivariate Processes}
\usage{
cov_flexMatern_mulv(
  locs,
  sigma.mat = NULL,
  p = NULL,
  rho = NULL,
  nu.mat = NULL,
  nu.vec = NULL,
  alpha.mat = NULL,
  alpha = NULL
)
}
\arguments{
\item{locs}{A matrix of locations or a list of matrices of locations}

\item{sigma.mat}{A numerical \code{p} by \code{p} matrix of collocated covariance coefficients at \code{NULL} by default}

\item{p}{A number of processes at \code{NULL} by default. Please use this argument with the next argument rho}

\item{rho}{A coefficient for defining a matrix of collocated covariance coefficients as an autoregressive covariance matrix. Please use this argument with the previous argument p}

\item{nu.mat}{A numerical \code{p} by \code{p} matrix of smoothness parameters at \code{NULL} by default}

\item{nu.vec}{A numeric vector at \code{NULL} by default. If nu.vec is used instead of nu.mat, the \code{(i, j)}-th entry a matrix of smoothness parametes is defined by \code{(nu_i + nu_j)/2}}

\item{alpha.mat}{A numerical \code{p} by \code{p} matrix of scale parameters}

\item{alpha}{Numeric at NULL by default. If alpha is used instead of alpha.mat, a matrix of scale parameters is defined by a constant matrix of alpha}
}
\value{
A flexible Matern covariance matrix
}
\description{
Flexible Matern covariance function for Multivariate Processes
}
\examples{
# Toy example for the univariate isotropic covariance function
n <- 5 ; d <- 2 ; p <- 1
locs <- matrix(runif(n * d), n, d)

covmat <- cov_flexMatern_mulv(locs = locs,
                              sigma.mat = diag(p),
                              nu.mat = matrix(0.5, p, p),
                              alpha.mat = matrix(1, p, p))

covmat
exp(-fields::rdist(locs))

isSymmetric(covmat)
eigen(covmat)$values


# Example for the bivariate Matern covariance function
n1 <- 4 ; n2 <- 3 ; d <- 2 ; p <- 2
locs1 <- matrix(runif(n1 * d), n1, d) ; locs2 <- matrix(runif(n2 * d), n2, d)
sigma.mat <- diag(p) ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)

covmat <- cov_flexMatern_mulv(locs = list(locs1 = locs1, locs2 = locs2),
                              sigma.mat = sigma.mat,
                              nu.mat = nu.mat,
                              alpha.mat = alpha.mat)

isSymmetric(covmat)
eigen(covmat)$values

covmat2 <- cov_flexMatern_biv(locs1 = locs1, locs2 = locs2,
                              sigma.mat = sigma.mat,
                              nu.mat = nu.mat,
                              alpha.mat = alpha.mat)

identical(covmat, covmat2)

# Example for comparison between two differnt versions
#      of covariance matrix of multivariate gaussian process

n <- 3 ; d <- 2 ; p <- 2
locs <- matrix(runif(n * d), n, d)
sigma.mat <- diag(p) ; nu.mat <- matrix(0.5, p, p) ; alpha.mat <- matrix(1, p, p)

covmat <- cov_flexMatern_mulv(locs = locs,
                              sigma.mat = sigma.mat,
                              nu.mat = nu.mat,
                              alpha.mat = alpha.mat)
covmat_each <- cov_GK_mulv(locs = locs,
                           sigma.mat = sigma.mat,
                           nu.mat = nu.mat,
                           alpha.mat = alpha.mat)
covmat_each_bruteforce <- cov_GK_mulv_bruteforce(locs = locs,
                                                 sigma.mat = sigma.mat,
                                                 nu.mat = nu.mat,
                                                 alpha.mat = alpha.mat)

covmat_matched <- matrix(NA, n * p, n * p)
for(i in 1:p){
  for(j in 1:p){
    covmat_matched[seq(from = 1 + n * (i - 1), to = n * i, by = 1),
                   seq(from = 1 + n * (j - 1), to = n * j, by = 1)] <-
                   covmat_each[seq(from = i, to = i + p * (n - 1), by = p),
                               seq(from = j, to = j + p * (n - 1), by = p)]
  }
}

covmat
covmat_each
covmat_each_bruteforce

identical(covmat_each, covmat_each_bruteforce)
identical(covmat, covmat_matched)
}
\references{
Bachoc, François, and Reinhard Furrer. "On the smallest eigenvalues of covariance matrices of multivariate spatial processes." Stat 5.1 (2016): 102-107.
}
