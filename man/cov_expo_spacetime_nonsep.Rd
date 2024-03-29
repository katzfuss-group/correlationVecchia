% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covmodels.R
\name{cov_expo_spacetime_nonsep}
\alias{cov_expo_spacetime_nonsep}
\title{Non-separable space-time exponential covariance function}
\usage{
cov_expo_spacetime_nonsep(locs, covparms)
}
\arguments{
\item{locs}{A matrix of locations = (1st coordinate of spatial location, 2nd coordinate of spatial location, temporal location)}

\item{covparms}{A numeric vector of covariance parameters = (sigma^2, kappa, a, c)}
}
\value{
A matrix with \code{n} rows and \code{n} columns, with the \code{(i, j)} entry containing the non-separable spatio-temporal covariance between observations \code{locs[i, ]} and \code{locs[j, ]} described in the paper below
}
\description{
From a location matrix \code{locs} and a vector with covariance parameters \code{covparms}, this function returns a non-separable space-time exponential covariance matrix.
}
\examples{
set.seed(5152021)

locs <- matrix(runif(30), 10, 3)
locs

covmat <- cov_expo_spacetime_nonsep(locs = locs, covparms = c(1, 0.75, 50, 25))
fields::image.plot(covmat)

covmat <- cov_expo_spacetime_nonsep(locs = locs, covparms = c(1, 0.75, 500, 2.5))
fields::image.plot(covmat)
}
\references{
Datta, Abhirup, et al. "Nonseparable dynamic nearest neighbor Gaussian process models for large spatio-temporal data with an application to particulate matter analysis." The annals of applied statistics 10.3 (2016): 1286.
}
