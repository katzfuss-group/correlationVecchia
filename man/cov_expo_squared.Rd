% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covmodels.R
\name{cov_expo_squared}
\alias{cov_expo_squared}
\title{Squared exponential covariance function}
\usage{
cov_expo_squared(locs, covparms)
}
\arguments{
\item{locs}{A numerical matrix with \code{n} rows and \code{d} columns. Each row of locs gives a point in R^d}

\item{covparms}{A numerical vector with covariance parameters in the form (variance, range)}
}
\value{
A matrix with \code{n} rows and \code{n} columns, with the \code{(i, j)} entry containing the squared exponential covariance between observations \code{locs[i, ]} and \code{locs[j, ]}
}
\description{
From a location matrix \code{locs} and a vector with covariance parameters \code{covparms}, this function returns a squared exponential covariance matrix.
}
\examples{
locs <- matrix(runif(10), 5, 2)
all.equal(cov_expo_squared(locs = locs, covparms = c(1, 1)),
          dnorm(fields::rdist(locs)) * sqrt(2 * pi))
}
