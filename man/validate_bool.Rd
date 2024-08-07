% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strings_utils.R
\name{validate_bool}
\alias{validate_bool}
\title{Boolean Expression Validation}
\usage{
validate_bool(x = "", all_names = "")
}
\arguments{
\item{x}{a string vector representing the boolean expression to be validated. Default is ""}

\item{all_names}{a character vector of scalars which are allowed to be part of the the boolean expression.}
}
\value{
x is returned if no exception is raised during validation process.
}
\description{
Helper to check if a boolean vector is valid.
}
\keyword{internal}
