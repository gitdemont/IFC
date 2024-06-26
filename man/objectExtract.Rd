% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objectExtract.R
\name{objectExtract}
\alias{objectExtract}
\title{Object Extraction}
\source{
For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}\cr
\verb{
BSD implementations of Bio-Formats readers and writers
%%
Copyright (C) 2005 - 2017 Open Microscopy Environment:
  - Board of Regents of the University of Wisconsin-Madison
  - Glencoe Software, Inc.
  - University of Dundee
%%
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
}
}
\usage{
objectExtract(ifd, param, verbose = FALSE, bypass = FALSE, ...)
}
\arguments{
\item{ifd}{list of sub elements of IFD data information extracted by \code{\link{getIFD}}. This argument can't be missing.}

\item{param}{object of class `IFC_param`, containing extraction parameters defined by \code{\link{objectParam}}.\cr
This argument is not mandatory but it may allow to save time for repeated image export on same file.
If this argument is missing, \code{\link{objectExtract}} will use named arguments in \code{'...'} as input to \code{\link{objectParam}} to control object extraction.\cr
However, if \code{'param'} is provided, \code{'...'} will be \strong{ignored}.}

\item{verbose}{whether to display information (use for debugging purpose). Default is \code{FALSE}.}

\item{bypass}{whether to bypass checks on \code{'ifd'} and \code{'param'}. Default is \code{FALSE}.}

\item{...}{
  Arguments passed on to \code{\link[=objectParam]{objectParam}}
  \describe{
    \item{\code{info}}{object of class `IFC_info`, rich information extracted by \code{\link{getInfo}}. 
This argument is not mandatory but it may allow to save time for repeated image export on same file.
If missing, the default, \code{'info'} will be extracted thanks to \code{'...'}.}
    \item{\code{mode}}{color mode export. Either \code{"rgb"}, \code{"gray"} or \code{"raw"}. Default is \code{"raw"}.
Note that \code{"raw"} is only possible when \code{'export'} is \code{"matrix"} or \code{"multi"}.}
    \item{\code{export}}{format mode export. Either \code{"file"}, \code{"matrix"}, \code{"base64"}, or \code{"multi"}. Default is \code{"matrix"}.}
    \item{\code{write_to}}{used when export is not \code{"matrix"} to compute exported file name or base64 id attribute.\cr
Exported file extension and base64 MIME type will be deduced from this pattern. Allowed export are \code{".bmp"}, \code{".jpg"}, \code{".jpeg"}, \code{".png"}, \code{".tif"}, \code{".tiff"}.
Note that \code{".bmp"} are faster but not compressed producing bigger data.\cr
Placeholders, if found, will be substituted:\cr
-\code{\%d}: with full path directory,\cr
-\code{\%p}: with first parent directory,\cr
-\code{\%e}: with extension (without leading .),\cr
-\code{\%s}: with shortname (i.e. basename without extension),\cr
-\code{\%o}: with object_id,\cr
-\code{\%c}: with channel_id (not possible when \code{'export'} is \code{"multi"}).\cr
A good trick is to use:\cr
-\code{"\%d/\%s/\%s_\%o_\%c.tiff"}, when \code{'export'} is \code{"file"},\cr
-\code{"\%d/\%s/\%s_\%o.tiff"}, when \code{'export'} is \code{"multi"},\cr
-\code{"\%o_\%c.bmp"}, when \code{'export'} is \code{"base64"}.\cr
Note that if missing and \code{'export'} is not \code{"file"}, \code{'write_to'} will be set to \code{"\%o_\%c.bmp"}.}
    \item{\code{base64_id}}{whether to add id attribute to base64 exported object. Default is \code{FALSE}.\cr
Only applied when export is \code{"base64"}.}
    \item{\code{base64_att}}{attributes to add to base64 exported object. Default is \code{""}.\cr
Only applied when export is \code{"base64"}. For example, use \code{"class='draggable'"}.\cr
Note that \code{id} (if \code{'base64_id'} is \code{TRUE}) and \code{width} and \code{height} are already used.}
    \item{\code{overwrite}}{only apply when \code{'export'} is \code{"file"} whether to overwrite file or not. Default is \code{FALSE}.}
    \item{\code{composite}}{character vector of image composite. Default is \code{""}, for no image composite.\cr
Should be like \code{"1.05/2.4/4.55"} for a composition of 5 perc. of channel 1, 40 perc. of channel 2 and 50 perc. of channel 55.\cr
Note that channels should have been acquired and final image composition should be 100 perc., otherwise an error is thrown.\cr
Note that each composite will be appended after \code{'selection'}.\cr
Note that composite will be forced to \code{""} when \code{'export'} is \code{"multi"}.}
    \item{\code{selection}}{physical channels to extract.\cr
Note that this parameter will be ordered.\cr
Default is \code{"all"} to extract all acquired channels.\cr
Use \code{"none"} to only extract composite.}
    \item{\code{random_seed}}{a list of elements to pass to \link[base]{set.seed} or a single value, interpreted as an integer, or NULL to be used when \code{'add_noise'} is set to \code{TRUE}. Default is \code{NULL}.
Note that \code{NA_integer_} or \code{list(seed = NA_integer_)} can be used to not call \link[base]{set.seed} at all.}
    \item{\code{size}}{a length 2 integer vector of final dimensions of the image, height 1st and width 2nd. Default is \code{c(0,0)} for no change.}
    \item{\code{force_width}}{whether to use information in \code{'info'} to fill size. Default is \code{TRUE}.
When set to \code{TRUE}, width of \code{'size'} argument will be overwritten.}
    \item{\code{removal}}{removal method: Either \code{"none"}, \code{"raw"}, \code{"clipped"}, \code{"masked"}, \code{"MC"}.\cr
-\code{"none"}, to keep image as is, no mask extraction will be performed resulting in faster extraction,\cr
-\code{"raw"}, to keep image as is, it provides a convenient way to retrieve \code{"raw"} value for the mask,\cr
-\code{"clipped"}, to remove clipped object from image,\cr
-\code{"masked"}, to only keep masked object from image,\cr
-\code{"MC"}, to keep MC masked object from image.
This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to \code{'selection'} and \code{'composite'} parameters.}
    \item{\code{add_noise}}{if \code{TRUE} adds normal noise to background using \pkg{Rcpp}. Default is \code{TRUE}.\cr
Note that it is better to set it to \code{FALSE} when \code{'removal'} is \code{"masked"} or \code{"MC"}. Doing so will allow to place masked object in a zero filled background,
otherwise background will still be filled with noise.
This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to \code{'selection'} and \code{'composite'} parameters.}
    \item{\code{full_range}}{only apply when \code{'mode'} is not \code{"raw"}, if \code{'full_range'} is \code{TRUE}, then object range will be considered as 0 to 4095, it is like \code{"raw"} \code{'mode'} but resulting in [0,4095] normalization to [0,1]. Default is \code{FALSE}.\cr
This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to \code{'selection'} and \code{'composite'} parameters.}
    \item{\code{force_range}}{only apply when \code{'mode'} is not \code{"raw"}, if \code{'force_range'} is \code{TRUE}, then range will be adjusted to object range in \code{[-4095,+inf]} resulting in normalization to [0,1]. Default is \code{FALSE}.\cr
This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to \code{'selection'} and \code{'composite'} parameters.\cr
Note that this parameter takes the precedence over \code{'full_range'}.}
    \item{\code{spatial_correction}}{only apply on RIF file, whether to apply spatial correction. Default is \code{FALSE}.}
  }}
}
\value{
A list (for every extracted objects) of list (for every exported channels) depending on '\code{param$export}' parameter:\cr
-\code{"matrix"}, a matrix when '\code{param$mode}' is set to \code{"raw"} or \code{"gray"} OR an array when '\code{param$mode}' is \code{"rgb"},\cr
-\code{"base64"}, a data-uri string,\cr
-\code{"file"}, an invisible file path corresponding to the location of exported file(s).\cr
-\code{"multi"}, an invisible file path corresponding to the location of exported file(s).
}
\description{
Extracts / Decompress objects stored in RIF or CIF Files.
}
\details{
When a mask is detected, \code{'add_noise'}, \code{'full_range'} and \code{'force_range'} are set to \code{FALSE} and range used will be forced to \code{[0,3]}.\cr\cr
Range of image is controlled by \code{'Images'} information from supplied \code{'info'} or as extracted by \code{\link{getInfo}} and will be returned as \code{'channels'} by \code{\link{objectParam}}.
In case \code{'mode'} is not \code{"raw"}, '\code{channels$xmin}', '\code{channels$xmax}', '\code{channels$gamma}' will be used for object extraction by \code{\link{objectExtract}} unless any of \code{'force_range'} or \code{'full_range'} is \code{TRUE}.\cr\cr
Experimental (as of v0.2.0.501): once returned by \code{\link{objectParam}}, those '\code{channels$xmin}' and '\code{channels$xmax}' can be manually adjusted to \code{]0,1[} so as to be used as \code{'probs'} argument to \link[stats]{quantile} to allow quantile normalization during object extraction (\code{\link{objectExtract}}) afterwards.
}
\note{
\code{'...'} will be used to compute \code{'param'} \strong{only} if \code{'param'} is \strong{not} provided.\cr
\strong{/!\\} When \code{'param'} is not provided, if not any of \code{'fileName'}, \code{'info'} can be found in \code{'...'}, then \code{attr(ifd, "fileName_image")} will be used as \code{'fileName'} input parameter to pass to \code{\link{objectParam}}.
}
\examples{
if(requireNamespace("IFCdata", quietly = TRUE)) {
  ## use a cif file
  file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
  cif_offs <- getOffsets(fileName = file_cif, fast = TRUE)
  ## extract infomation
  info <- getInfo(fileName = file_cif, from = "analysis")
  ## retrieve number of objects stored
  nobj <- as.integer(info$objcount)
  ## randomly subset the offsets of at most 5 "img" objects
  sel = sample(0:(nobj-1), min(5, nobj))
  sub_offs <- subsetOffsets(cif_offs, objects = sel, image_type = "img")
  ## read IFDs from these "img" objects
  IFDs <- getIFD(fileName = file_cif, offsets = sub_offs)
  ## extract raw data of these"img" objects to matrix
  raw = objectExtract(ifd = IFDs, info = info, mode = "raw", 
                      export = "matrix")
  ## extract base64 "rgb" colorized version of these "img" objects to base64
  b64 = objectExtract(ifd = IFDs, info = info, mode = "rgb", 
                      export = "base64", base64_id = TRUE,
                      write_to = "example_\%o_\%c.bmp")
  ## use DisplayGallery to show the first "img" objects and play with ... extra parameters
  ## force_range, add_noise, selection, composite, see objectParam
  DisplayGallery(info = info, offsets = cif_offs, objects = sel,
                 base64_id = TRUE, write_to = "example_\%o_\%c.bmp",
                 force_range = c(FALSE,TRUE,FALSE,TRUE), add_noise = FALSE,
                 selection = c(1,2,4,6), composite = "1.7/4.3")
} else {
  message(sprintf('Please run `install.packages("IFCdata", repos = "\%s", type = "source")` \%s',
                  'https://gitdemont.github.io/IFCdata/',
                  'to install extra files required to run this example.'))
}
}
