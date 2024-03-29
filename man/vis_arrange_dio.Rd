% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{vis_arrange_dio}
\alias{vis_arrange_dio}
\title{Visualization}
\usage{
vis_arrange_dio(
  vdat1,
  vdat2,
  legend,
  color,
  shape,
  ftn = log10,
  xlim = c("m", "m"),
  ylim = c("log10(KL)", "log10(KL)"),
  xlab = list(xl1 = NULL, xl2 = NULL),
  ylab = list(yl1 = NULL, yl2 = NULL),
  alpha = 0.7,
  size.point = 2,
  size.line = 1,
  size.legend = 14,
  size.lab = 14,
  size.text = 12,
  size.margin = c(5.5, 20, 5.5, 5.5)
)
}
\arguments{
\item{vdat1}{First dataset consisting of the columns: index, variable (for example, m), approx_1, approx_2, ...}

\item{vdat2}{Second dataset consisting of the columns: index, variable (for example, a), approx_1, approx_2, ...}

\item{legend}{legend}

\item{color}{color}

\item{shape}{shape}

\item{ftn}{scale function}

\item{xlim}{xlim}

\item{ylim}{ylim}

\item{xlab}{xlab}

\item{ylab}{ylab}

\item{alpha}{alpha}

\item{size.point}{size.point}

\item{size.line}{size.line}

\item{size.legend}{size.legend}

\item{size.lab}{size.lab}

\item{size.text}{size.text}

\item{size.margin}{size.margin}
}
\value{
output of grid.arrange(). Use plot(vis_arrange(...)) to plot it
}
\description{
Visualization
}
\examples{
1 + 1
}
