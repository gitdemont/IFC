% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildImage.R
\name{buildImage}
\alias{buildImage}
\title{IFC Image Coercion}
\usage{
buildImage(physicalChannel = 1, ..., BF = FALSE, MODE = 1)
}
\arguments{
\item{physicalChannel}{the channel. Default is 1. Allowed are [1-12].}

\item{...}{other arguments to be passed. See members in returned value.}

\item{BF}{should 'physicalChannel' channel be considered as brightfield. Default is FALSE.}

\item{MODE}{collection mode (as retrieved by getInfo) determining the range. Default is 1.}
}
\value{
a list whose members are "name", "color", "physicalChannel", "xmin", "xmax",
"xmid", "ymid", "scalemin", "scalemax", "tokens", "baseimage", "function", "saturation".
}
\description{
Helper to build a list image values to allow export
}
\keyword{internal}
