% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covmodels.R
\name{cov_expo_aniso}
\alias{cov_expo_aniso}
\title{Anisotropic exponential covariance function}
\usage{
cov_expo_aniso(locs, covparms, a = NULL)
}
\arguments{
\item{locs}{A numerical matrix with \code{n} rows and \code{d} columns. Each row of locs gives a point in R^d}

\item{covparms}{A numerical vector with covariance parameters. It can be of two different forms: One is the 3-dim'l vector (variance, range, degree of anisotropy) and the other is simply 2-dim'l vector (variance, range). In the latter case, a vector of anisotropy must be specified}

\item{a}{A vector of anisotropy. At \code{NULL} by default}
}
\value{
A matrix with \code{n} rows and \code{n} columns, with the \code{(i, j)} entry containing the anisotropic exponenital covariance between observations \code{locs[i, ]} and \code{locs[j, ]}
}
\description{
From a location matrix \code{locs} and a vector with covariance parameters \code{covparms}, this function returns an anisotropic exponential covariance matrix.
}
\section{Parametrization}{
 The covariance parameters are a variance \eqn{\sigma^2}, range \eqn{r}, and degree (vector) of anisotropy \eqn{\alpha}. The general form of the covariance is \deqn{ C(x, y) = \sigma^2 exp( || A ( x - y ) || / r)} where \eqn{x} and \eqn{y} are locations in R^d. If \eqn{\alpha} is a number, then A is a diagonal matrix \eqn{diag( \sqrt \alpha , 1 , ... , 1 )}. On the other hand, if \eqn{\alpha} is a vector, then A is \eqn{diag( \sqrt \alpha_1, \sqrt \alpha_2, ... , \sqrt \alpha_d )}.
}

\examples{
# grid locations
cov_expo_aniso(locs = as.matrix(expand.grid(c(0.25, 0.75), c(0.25, 0.75))),
               covparms = c(1, 0.1, 10))
cov_expo_aniso(locs = as.matrix(expand.grid(c(0.25, 0.75), c(0.25, 0.75))),
               covparms = c(1, 0.1), a = c(10, 1)) # equal to the above

# randomly selected locations
cov_expo_aniso(locs = matrix(runif(8), 4, 2),
               covparms = c(1, 0.1, 10))
cov_expo_aniso(locs = matrix(runif(8), 4, 2),
               covparms = c(1, 0.1), a = c(10, 1)) # equal to the above
}
