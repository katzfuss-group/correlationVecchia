% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sortSparse_Rcpp}
\alias{sortSparse_Rcpp}
\title{Florian's Algorithm without unmeasured inputs (C++ version)}
\usage{
sortSparse_Rcpp(x, rho, initInd, distype, fstr, covparms)
}
\arguments{
\item{x}{A numerical matrix of inputs (locations)}

\item{rho}{hyperparameter rho (instead of m)}

\item{initInd}{initial input number}

\item{distype}{euclidean or correlation}

\item{fstr}{covariance function name}

\item{covparms}{A numerical vector with covariance parameters}
}
\value{
List. Check it out!
}
\description{
Florian's Algorithm without unmeasured inputs (C++ version)
}
