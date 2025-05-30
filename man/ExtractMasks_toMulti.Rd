% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtractToFuns.R
\name{ExtractMasks_toMulti}
\alias{ExtractMasks_toMulti}
\title{Shortcut for Batch Masks Extraction to Multichannel Tiff}
\usage{
ExtractMasks_toMulti(
  ...,
  objects,
  offsets,
  display_progress = TRUE,
  write_to,
  overwrite = FALSE
)
}
\arguments{
\item{...}{arguments to be passed to \code{\link{objectExtract}} with the exception of \code{'ifd'} and \code{'bypass'}(=\strong{TRUE}).\cr
\strong{/!\\} If not any of \code{'fileName'}, \code{'info'} and \code{'param'} can be found in \code{'...'} then \code{attr(offsets, "fileName_image")} will be used as \code{'fileName'} input parameter to pass to \code{\link{objectParam}}.}

\item{objects}{integer vector, IDEAS objects ids numbers to use.
This argument is not mandatory, if missing, the default, all objects will be used.}

\item{offsets}{object of class `IFC_offset`.
This argument is not mandatory but it may allow to save time for repeated image export on same file.\cr
If \code{'offsets'} are not provided, extra arguments can also be passed with \code{'...'} to \code{\link{getOffsets}}.}

\item{display_progress}{whether to display a progress bar. Default is TRUE.}

\item{write_to}{used to compute respectively exported file name.\cr
Exported \code{"multi"} extension will be deduced from this pattern. Allowed export are \code{".tif"}, \code{".tiff"}.
Placeholders, if found, will be substituted:\cr
-\code{\%d}: with full path directory,\cr
-\code{\%p}: with first parent directory,\cr
-\code{\%e}: with extension (without leading .),\cr
-\code{\%s}: with shortname (i.e. basename without extension),\cr
-\code{\%o}: with object_id,\cr
A good trick is to use: \code{"\%d/\%s/\%s_\%o.tiff"}.}

\item{overwrite}{whether to overwrite file or not. Default is \code{FALSE}.}
}
\value{
It invisibly returns a list of exported file path of corresponding to objects extracted.
}
\description{
Function to shortcut extraction of masks to multichannel tiff ! excludes image.
}
\details{
If \code{'param'} is provided in \code{'...'}:\cr
-\code{'param$export'<-"multi"}, \code{'param$mode'<-"raw"} and \code{'param$overwrite'<-'overwrite'} will be overwritten.\cr
-if \code{'write_to'} is not missing, \code{'param$write_to'<-'write_to'} will be overwritten. Otherwise, \code{'param$write_to'} will be used \strong{only} if \code{'param$export'} was \code{"multi"}.\cr\cr
\code{'write_to'} has to be provided if \code{'param'} can't be found in \code{'...'} or if \code{'param$export'} was not \code{"multi"}.
}
\note{
Arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractMasks_toMulti}} input arguments.
}
