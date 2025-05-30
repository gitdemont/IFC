% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExportToGallery.R
\name{ExportToGallery}
\alias{ExportToGallery}
\title{Gallery Export}
\usage{
ExportToGallery(
  ...,
  objects,
  offsets,
  image_type = "img",
  layout,
  export = c("file", "matrix", "base64")[2],
  write_to,
  base64_id = FALSE,
  base64_att = "",
  overwrite = FALSE,
  main = "",
  add_channels = TRUE,
  add_ids = 1,
  add_lines = 2,
  bg_color = "grey20",
  dpi = 300,
  scale = list(),
  extract_max = 10,
  sampling = FALSE,
  display_progress = TRUE
)
}
\arguments{
\item{...}{arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'bypass'(=TRUE).\cr
If 'param' is provided 'mode'(="rgb") and the above parameters will be overwritten.\cr
If 'offsets' are not provided extra arguments can also be passed with ... \code{\link{getOffsets}}.\cr
/!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.}

\item{objects}{integer vector, IDEAS objects ids numbers to use.
This argument is not mandatory, if missing, the default, all objects will be used.}

\item{offsets}{object of class `IFC_offset`. 
This argument is not mandatory but it may allow to save time for repeated image export on same file.}

\item{image_type}{image_type of desired offsets. Either "img" or "msk". Default is "img".}

\item{layout}{a character vector of [acquired channels + 'composite' images] members to export. Default is missing to export everything.
Note that members can be missing to be removed from final gallery export.
Note that members not found will be automatically removed and a warning will be thrown.}

\item{export}{export format. Either "file", "matrix", "base64". Default is "matrix".}

\item{write_to}{used when 'export' is "file" or "base64" to compute respectively filename or base64 id attribute.
Exported type will be deduced from this pattern. Allowed export are '.bmp', '.jpg', '.jpeg', '.png', '.tif', '.tiff'.\cr
'.pdf' can also be used when 'export' is not "base64".
Note that '.bmp' are faster but not compressed producing bigger data.\cr
Placeholders, if found, will be substituted:\cr
-\%d: with full path directory\cr
-\%p: with first parent directory\cr
-\%e: with extension (without leading .)\cr
-\%s: with shortname (i.e. basename without extension)\cr
-\%o: with objects (at most 10, will be collapse with "_", if more than one).\cr
-\%c: with channel_id (will be collapse with "_", if more than one, composite in any will be bracketed).
A good trick is to use:\cr
-"\%d/\%s_gallery_Obj[\%o]_Ch[\%c].tiff", when 'export' is "file"\cr
-"\%s_gallery.bmp", when 'export' is "base64".\cr
Note that if missing and 'export' is not "file", 'write_to' will be set to "\%s_gallery.bmp".}

\item{base64_id}{whether to add id attribute to base64 exported object. Default is TRUE.\cr
Only applied when 'export' is "base64".}

\item{base64_att}{attributes to add to base64 exported object. Default is "".\cr
Only applied when export is "base64". For example, use "class=draggable".\cr
Note that id (if base64_id is TRUE) and width and height are already used.}

\item{overwrite}{whether to overwrite file or not. Default is FALSE.}

\item{main}{main title that will be displayed on top center of the image.
If too large it will be clipped.}

\item{add_channels}{whether to add channels names. Default is TRUE.}

\item{add_ids}{integer, indice of column to mark objects ids number. Default is 1.
If add_ids < 1, no ids are added.}

\item{add_lines}{integer, size of separating lines between objects. Default is 1.
If add_lines < 1, no separating lines are added.}

\item{bg_color}{background color for main, channels and separating lines. Default is "grey20".}

\item{dpi}{integer, the resolution of the image in DPI (dots per inch). Default is 300.\cr
Please note that whatever this parameter is final resolution will be 96 dpi.\cr
However image will be scaled according this parameter and magnification factor will be equal to this parameter divided by 96.}

\item{scale}{a named list whose members are 'size', 'style', 'color', 'xoff', 'yoff'. Default is list() to draw no scale. Otherwise,\cr
-'size' positive integer. Scale's bar size in micro-meter. Default is '7'.\cr
This parameter can't be lesser than 6px nor higher than image width + scale text.\cr
-'style' a character string. Scale's bar style, either "dash" or "line". Default is "dash".\cr
-'color' a character string. color of the scale. Default is "white".\cr
-'xoff' positive integer. x offset in image to draw scale, starting from bottom left corner.\cr
-'yoff' positive integer. y offset in image to draw scale, starting from bottom left corner.\cr
\code{'scale'=TRUE}, can eventually be used to draw scale with default parameters.}

\item{extract_max}{maximum number of objects to extract. Default is 10. Use +Inf to extract all.}

\item{sampling}{whether to sample objects or not. Default is FALSE.}

\item{display_progress}{whether to display a progress bar. Default is TRUE.}
}
\value{
Depending on 'export':\cr
-"matrix", a rgb array,\cr
-"base64", a data-uri string,\cr
-"file", invisibly returns path of exported file.\cr
with ids attribute corresponding to the objects exported.
}
\description{
Exports gallery of `IFC_img` / `IFC_msk` objects
}
\details{
arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExportToGallery}} input arguments.
TRICK: for exporting only ONE 'objects', set 'add_channels' = FALSE, 'add_ids' >= 1, 'force_width' = FALSE, 'dpi' = 96; this allows generating image with its original size incrusted with its id number.
}
