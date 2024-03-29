% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{conditioning_rho_Rcpp}
\alias{conditioning_rho_Rcpp}
\title{conditioning based on Florian's algorithm (C++ version)}
\usage{
conditioning_rho_Rcpp(indvec, condvec, P, maxsize, x, distype, fstr, covparms)
}
\arguments{
\item{indvec}{indvec}

\item{condvec}{condvec}

\item{P}{P}

\item{maxsize}{maxsize}

\item{x}{x}

\item{distype}{euclidean or correlation}

\item{fstr}{covariance function name}

\item{covparms}{A numerical vector with covariance parameters}
}
\value{
Matrix. Check it out!
}
\description{
conditioning based on Florian's algorithm (C++ version)
}
