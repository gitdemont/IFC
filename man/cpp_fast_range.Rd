% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cpp_fast_range}
\alias{cpp_fast_range}
\title{Use Rcpp for Range}
\arguments{
\item{x_}{a Nullable NumericVector.}
}
\value{
a NumericVector.
}
\description{
Determines range of numeric vector
}
\details{
the behaviour is the same as R base::range(x_, na.rm = TRUE, finite = TRUE) without creating warnings
}
\keyword{internal}
