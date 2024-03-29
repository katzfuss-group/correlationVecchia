% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generation.R
\name{generate_hnm}
\alias{generate_hnm}
\title{Generating a realization of hierarchical normal model (non-Euclidean example)}
\usage{
generate_hnm(mu0, sig20, k, sig2s = NULL, n = NULL, mlt = 2)
}
\arguments{
\item{mu0}{Mean of the root normal model}

\item{sig20}{Variance of the root normal model}

\item{k}{Depth of the hierarchical normal model}

\item{sig2s}{Numeric vector of length of k. The i-th entry is the variance of the normal models in level i. As \code{NULL} (rep(1, k)) by default.}

\item{n}{Number of leaves. \code{NULL} (mlt^k) as default.}

\item{mlt}{Number of children for each node}
}
\value{
\code{generate_hnm} returns
\itemize{
\item{\code{mu}: } A numeric vector of means of the parent normal models
\item{\code{y}: } A numeric vector of leaves
\item{\code{mu0}: } Mean of the root normal model
\item{\code{sig20}: } Variance of the root normal model
\item{\code{k}: } Depth of the hierarchical normal model
\item{\code{sig2s}: } Numeric vector of length of k. The i-th entry is the variance of the normal models in level i
\item{\code{n}: } Number of leaves
\item{\code{mlt}: } Number of children for each node
}
}
\description{
Generating a realization of hierarchical normal model (non-Euclidean example)
}
\examples{
output <- generate_hnm(mu0 = 0, sig20 = 1, k = 3, mlt = 2)

length(output$y)

covmat <- covmat_hnm(output)
cormat <- cormat_hnm(output)

fields::imagePlot(cormat)
}
