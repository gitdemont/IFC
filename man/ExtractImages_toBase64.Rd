% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtractImages_toBase64.R
\name{ExtractImages_toBase64}
\alias{ExtractImages_toBase64}
\title{Shorcut for Batch Images Extraction to Base64}
\usage{
ExtractImages_toBase64(
  ...,
  objects,
  offsets,
  display_progress = TRUE,
  mode = c("rgb", "gray")[1]
)
}
\arguments{
\item{...}{arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'bypass'(=TRUE).\cr
If 'param' is provided 'export'(="base64") and 'mode' will be overwritten.\cr
If 'offsets' are not provided extra arguments can also be passed with ... \code{\link{getOffsets}}.\cr
/!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.}

\item{objects}{integer vector, IDEAS objects ids numbers to use.
This argument is not mandatory, if missing, the default, all objects will be used.}

\item{offsets}{object of class `IFC_offset`. 
This argument is not mandatory but it may allow to save time for repeated image export on same file.}

\item{display_progress}{whether to display a progress bar. Default is TRUE.}

\item{mode}{(\code{\link{objectParam}} argument) color mode export. Either "rgb", "gray". Default is "rgb".}
}
\value{
A list of base64 encoded images corresponding to objects extracted.
}
\description{
Function to shortcut extraction, normalization and eventually colorization of images to matrix ! excludes mask.
}
\details{
arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractImages_toBase64}} input arguments.
}
