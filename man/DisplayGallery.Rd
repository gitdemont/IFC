% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DisplayGallery.R
\name{DisplayGallery}
\alias{DisplayGallery}
\title{Gallery Display}
\usage{
DisplayGallery(
  ...,
  objects,
  offsets,
  image_type = "img",
  layout,
  name = "DisplayGallery",
  caption = FALSE,
  pageLength = 10L,
  pdf_pageSize = "A2",
  pdf_pageOrientation = "landscape",
  pdf_image_dpi = 96,
  extract_max = 10,
  sampling = FALSE,
  display_progress = TRUE,
  mode = c("rgb", "gray")[1]
)
}
\arguments{
\item{...}{arguments to be passed to \code{\link{objectExtract}} with the exception of \code{'ifd'} and \code{'bypass'}(=\strong{TRUE}).\cr
If \code{'param'} is provided \code{'export'}(=\strong{"base64"}) and the above parameters will be overwritten.\cr
If \code{'offsets'} are not provided extra arguments can also be passed with \code{...} \code{\link{getOffsets}}.\cr
\strong{/!\\} If not any of \code{'fileName'}, \code{'info'} and \code{'param'} can be found in \code{'...'} then \code{attr(offsets, "fileName_image")} will be used as \code{'fileName'} input parameter to pass to \code{\link{objectParam}}.}

\item{objects}{integer vector, IDEAS objects ids numbers to use.
This argument is not mandatory, if missing, the default, all objects will be used.}

\item{offsets}{object of class `IFC_offset`. 
This argument is not mandatory but it may allow to save time for repeated image export on same file.}

\item{image_type}{image_type of desired offsets. Either \code{"img"} or \code{"msk"}. Default is \code{"img"}.}

\item{layout}{a character vector of [acquired channels + 'composite' images] members to export. Default is missing to export everything.\cr
Note that members can be missing to be removed from final display.\cr
Note that members not found will be automatically removed and a warning will be thrown.}

\item{name}{id of the datatable container. Default is \code{"DisplayGallery"}.}

\item{caption}{whether to display caption name or not. Default is \code{FALSE}.}

\item{pageLength}{integer, number of objects to display per page. Default is \code{10}.}

\item{pdf_pageSize}{string, page dimension when exporting to pdf. Default is \code{"A2"}.}

\item{pdf_pageOrientation}{string, page orientation when exporting to pdf. Default is \code{"landscape"}. Allowed are \code{"landscape"} or \code{"portrait"}.}

\item{pdf_image_dpi}{integer, desired image resolution. Default is \code{96}, for full resolution.}

\item{extract_max}{maximum number of objects to extract. Default is \code{10}. Use \code{+Inf} to extract all.}

\item{sampling}{whether to sample objects or not. Default is \code{FALSE}.}

\item{display_progress}{whether to display a progress bar. Default is \code{TRUE}.}

\item{mode}{(\code{\link{objectParam}} argument) color mode export. Either \code{"rgb"} or \code{"gray"}. Default is \code{"rgb"}.}
}
\value{
it invisibly returns a list whose members are:\cr
-data, data for DT::datatable(),\cr
-args, associated arguments to pass to DT::datatable().
}
\description{
Displays gallery of `IFC_img` / `IFC_msk` objects
}
\details{
arguments of \code{\link{objectExtract}} will be deduced from \code{\link{DisplayGallery}} input arguments.\cr
Please note that PDF export link will be available if \code{'write_to'} does not result in a \code{"bmp"}.\cr
Please note that viewing PDF with gallery exported as \code{"tiff"} may depend on browser capabilities.\cr
Please note that a warning may be sent if gallery to display contains large amount of data. This is due to use of datatable() from \pkg{DT}.\cr
\verb{
In instance$preRenderHook(instance) :
It seems your data is too big for client-side DataTables. You may consider server-side processing: http://rstudio.github.io/DT/server.html
}
For these reasons, it may be better to use \code{"png"} extension to display images.
}
\examples{
if(requireNamespace("IFCdata", quietly = TRUE)) {
  ## use a cif file
  file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
  cif <- ExtractFromXIF(fileName = file_cif)
  info <- getInfo(fileName = file_cif, from = "analysis")
  ## randomly show at most 10 "img" objects from file
  DisplayGallery(info = info, image_type = "img", extract_max = 10,
                 sampling = TRUE, write_to = "example.bmp")
} else {
  message(sprintf('Please run `install.packages("IFCdata", repos = "\%s", type = "source")` \%s',
                  'https://gitdemont.github.io/IFCdata/',
                  'to install extra files required to run this example.'))
}
}
