% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cpp_pnt_in_gate}
\alias{cpp_pnt_in_gate}
\title{Point in Gate}
\arguments{
\item{pnts}{NumericMatrix, a 2-columns matrix with (x and y) coordinates of the points of interest.}

\item{gate}{NumericMatrix, a 2-columns matrix defining polygon vertices or ellipse boundaries.}

\item{algorithm}{int, used for computation. Default is 1.\cr
1: Trigonometry.\cr
2: Special case = axes-aligned rectangle.\cr
3: Special case = axes-aligned ellipse.}

\item{epsilon}{double, epsilon threshold value. Default is 0.000000000001}
}
\description{
This function checks if points lie in a polygon or ellipse.
}
\keyword{internal}
