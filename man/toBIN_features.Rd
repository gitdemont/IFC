% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toBIN_features.R
\name{toBIN_features}
\alias{toBIN_features}
\title{IFC_features Raw Conversion}
\usage{
toBIN_features(
  features,
  w_con = raw(),
  endianness = .Platform$endian,
  verbose = FALSE,
  display_progress = TRUE,
  title_progress = "",
  ...
)
}
\arguments{
\item{features}{an `IFC_features` object.}

\item{w_con}{a connection opened for writing. Default is raw().}

\item{endianness}{The endian-ness ("big" or "little") of the target system for the file. Default is .Platform$endian.\cr
Endianness describes the bytes order of data stored within the files. This parameter may not be modified.}

\item{verbose}{whether to display message about current action. Default is FALSE.}

\item{display_progress}{whether to display a progress bar. Default is TRUE.}

\item{title_progress}{character string, giving the title of the progress bar. Default is "".}

\item{...}{other arguments to be passed.}
}
\value{
a raw vector of features binaries.
}
\description{
Helper to convert features (`IFC_features` object) to raw vector.
}
\keyword{internal}
