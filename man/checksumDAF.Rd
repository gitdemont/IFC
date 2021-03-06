% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checksumDAF.R
\name{checksumDAF}
\alias{checksumDAF}
\title{DAF File Checksum}
\usage{
checksumDAF(fileName, endianness = .Platform$endian, ...)
}
\arguments{
\item{fileName}{path to file.}

\item{endianness}{The endian-ness ("big" or "little") of the target system for the file. Default is .Platform$endian.\cr
Endianness describes the bytes order of data stored within the files. This parameter may not be modified.}

\item{...}{other arguments to be passed.}
}
\description{
This function returns CIF checksum computed from images values found in DAF
Checksum is the sum of img IFDs (Image Field Directory) offsets of objects 0, 1, 2, 3 and 4.
}
\keyword{internal}
