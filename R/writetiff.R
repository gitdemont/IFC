################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2024 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2020                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,             #
#                     Jean-Pierre Marolleau, Loïc Garçon,                      #
#                     INSERM, UPD, CHU Amiens                                  #
#                                                                              #
# DISCLAIMER:                                                                  #
# -You are using this package on your own risk!                                #
# -We do not guarantee privacy nor confidentiality.                            #
# -This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or #
# contributors be liable for any direct, indirect, incidental, special,        #
# exemplary, or consequential damages (including, but not limited to,          #
# procurement of substitute goods or services; loss of use, data, or profits;  #
# or business interruption) however caused and on any theory of liability,     #
# whether in contract, strict liability, or tort (including negligence or      #
# otherwise) arising in any way out of the use of this software, even if       #
# advised of the possibility of such damage.                                   #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with IFC. If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################

################################################################################
#             functions described hereunder are experimental                   #
#              inputs and outputs may change in the future                     #
################################################################################

#' @title Image Setup for TIFF Export
#' @name setupimage
#' @description
#' Perpares image for TIFF export.
#' @param image a matrix or array or representing the image:\cr
#' -2D matrix [h,w] will be written as single grayscale [h,w] IFD (Image Field Directory),\cr
#' -3D array [h,w,c] will be written as single IFD of multichannel [h,w,c],\cr
#' -4D array [h,w,c,f] will be written as multiple (multiframe) IFDs of multichannel [h,w,c],\cr
#' with h=height, w=width, c=channel, f=frame.\cr
#' Eventually, a raw vector with `dims`, `what` and `comp` attributes can be directly passed asis (allowing the use \code{compression} other than deflate).
#' @param what bits mode used to store image. Allowed are \code{"uint8"}, \code{"int8"}, \code{"uint16"}, \code{"int16"}, \code{"uint32"}, \code{"int32"}, \code{"float"} and \code{"double"}.
#' @param compression whether image should be lossless compressed with deflate algorithm. Default is FALSE.
#' @param endianness The endian-ness ("big" or "little") of the return object. Default is .Platform$endian.\cr
#' @return It returns a raw vector with `dims`, `what` and `comp` attributes.
#' @keywords internal
setupimage <- function(image, what, compression, endianness = .Platform$endian) {
  if(all(c("dims", "what", "comp") %in% names(attributes(image)))) return(image) # allow to pass other compression format
  what = match.arg(arg = what, choices = c("uint8","int8","uint16","int16","uint32","int32","float","double"), several.ok = FALSE)
  compression = na.omit(compression)
  stopifnot(length(compression) == 1,
            length(endianness) == 1 && endianness %in% c("little","big"))
  
  # determine type and bits
  z = switch(what, "int8" = 1, "uint8" = 2, "int16" = 3, "uint16" = 4, "int32" = 5, "uint32" = 6, "float" = 7, "double" = 8, 0)
  ans = cpp_cast_image(image, z, !identical(endianness, .Platform$endian))
  a = attributes(ans)
  
  if(compression == 8) {
    ans = memCompress(ans, type = "gzip")
  } else {
    if(compression != 1) warning("compression[",compression,"] not handled, will be set to 1")
    compression = 1
  }
  attributes(ans) <- a
  attr(ans, "comp") = compression
  return(ans)
}

#' @title Tiff Writter
#' @name writemulti
#' @description
#' Writes TIFF.
#' @param image a matrix or array or representing the image:\cr
#' -2D matrix [h,w] will be written as single grayscale [h,w] IFD (Image Field Directory),\cr
#' -3D array [h,w,c] will be written as single IFD of multichannel [h,w,c] grayscale or RGB if \code{rgb = TRUE},\cr
#' -4D array [h,w,c,f] will be written as multiple (multiframe) IFD of multichannel [h,w,c] grayscale or RGB if \code{rgb = TRUE},\cr
#' with h=height, w=width, c=channel, f=frame.
#' @param write_to file name or a raw vector. Default is raw(0).
#' @param tags list of extra tags to be included to IFD. Expecting a list whose sub-elements are list containing:\cr
#' -'tag' uint16_t, tag number, it should be [1-65535],\cr
#' -'typ' uint16_t typ number it should be [1-12],\cr
#' -'map' vector of values to write, it should not be empty and should be even for 'typ' 5 and 10.\cr
#' Default provides an example as how to include \code{"N/A"} image description (tag=270) and date (tag=306)
#' @param as.rgb whether to write image as RGB. Default is FALSE. It will be possible only for array where 3rd dimension length is 3.
#' @param compression whether image should be lossless compressed with deflate algorithm. Default is FALSE.
#' @param what bits mode used to store image. Default is \code{"uint16"}. Allowed are \code{"uint8"}, \code{"int8"}, \code{"uint16"}, \code{"int16"}, \code{"uint32"}, \code{"int32"}, \code{"float"} and \code{"double"}.
#' @param endianness The endian-ness ("big" or "little") of the return object. Default is .Platform$endian.\cr
#' Endianness describes the bytes order of data stored within the files.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @param verbosity quantity of information displayed when verbose is TRUE; 1: normal, 2: rich. Default is 1.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   ## extract image values from cif file
#'   RAW = ExtractImages_toMatrix(fileName = file_cif, mode = "raw", objects = 3,
#'                                force_width = FALSE, size = c(0,0))
#'   ## create multichannel tiff file
#'   file_tif = writemulti(image = array(unlist(RAW[[1]], use.names = FALSE, recursive = TRUE),
#'                                       dim = c(dim(RAW[[1]][[1]]), length(RAW[[1]]))),
#'                         write_to = tempfile(fileext = ".tif"))
#'   ## inspect exported file
#'   IFD = getIFD(file_tif)
#'   sapply(IFD[[1]]$tags, FUN = function(x) x)
#'   ## create multiframe tiff file
#'   file_tif = writemulti(image = array(unlist(RAW[[1]], use.names = FALSE, recursive = TRUE),
#'                                       dim = c(dim(RAW[[1]][[1]]), 1, length(RAW[[1]]))),
#'                         write_to = tempfile(fileext = ".tif"))
#'   ## inspect exported file
#'   IFDs = c(getIFD(file_tif),getIFD(file_tif, "all"))
#'   lapply(IFDs, FUN = function(i) sapply(i$tags, FUN = function(x) x))
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return It invisibly returns full path of exported file or a raw vector when \code{write_to = raw(0)} otherwise.
#' @export
writemulti <- function(image, write_to = raw(0),
                       tags = list(list(tag = 270, typ = 2, map = "N/A"), # ImageDescription
                                   list(tag = 306, typ = 2, map = "N/A")),# Date
                       as.rgb = FALSE, compression = TRUE,  what = "int16", endianness = .Platform$endian,
                       verbose = FALSE, verbosity = 1) {
  # TODO maybe add something to include mask ?
  stopifnot(typeof(write_to) %in% c("character", "raw"),
            length(as.rgb) == 1 && as.rgb %in% c(TRUE,FALSE),
            length(verbose) == 1 && verbose %in% c(TRUE,FALSE),
            length(verbosity) == 1 && verbosity %in% c(1,2),
            length(endianness) == 1 && endianness %in% c("little","big"))
  if(typeof(write_to) != "raw") {
    stopifnot(length(write_to) == 1)
    write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
  }
  dd = dim(image)
  if(all(c("dims", "what", "comp") %in% names(attributes(image)))) {
    as.is = TRUE
    dd = rev(attr(image, "dims")[-1])
  } else {
    stopifnot(length(compression) == 1 && compression %in% c(TRUE,FALSE))
    as.is = FALSE
  }
  if((length(dd) < 2) || (length(dd) > 4)) stop("expecting [h,w], [h,w,c] or [h,w,c,f] image")
  nb_chan = ifelse(length(dd) <= 2, 1, dd[3])
  
  # TIFF Magic Number & endian +  IFD starting position
  pos = 8L
  P1 = cpp_uint32_to_raw(pos)
  if(identical(endianness, "little")) {
    M1 = as.raw(c(0x49,0x49))
    M2 = as.raw(c(0x2a,0x00))
  } else {
    M1 = as.raw(c(0x4D,0x4D))
    M2 = as.raw(c(0x00,0x2a))
  }
  if(!identical(endianness, .Platform$endian)) P1 = rev(P1)
  ans = c(M1,M2,P1)
  
  # verbose info
  verb = ifelse(verbose & (verbosity==2), TRUE, FALSE)
  
  # frame / channels / rgb
  if(as.rgb && nb_chan != 3) {
    as.rgb = FALSE
    warning("can't create RGB tiff with ", nb_chan, " channels != 3")
  }
  d = dd
  while(length(d) < 4) d = c(d, 1)
  if(verbose) cat("Writing ", what, ifelse(as.rgb, " RGB", " Grayscale"), " image of ",
                  paste0(paste0(c("height=", "width=", "channel(s)=", "frame(s)="), d), collapse=", "), "\n", sep = "")
  if(as.is || length(dd) <= 3 || (length(dd) == 3 && as.rgb)) {
    ans = c(ans, cpp_writeIFD(setupimage(image, what, ifelse(compression, 8 ,1), endianness), tags, pos, endianness, as.rgb, TRUE, verb))
  } else {
    for(iframe in seq_len(dd[4])) {
      ret = cpp_writeIFD(setupimage(image[,,,iframe], what, ifelse(compression, 8, 1), endianness), tags, pos, endianness, as.rgb, iframe == dd[4], verb)
      pos = attr(ret, "offset")
      ans = c(ans, ret)
    }
  }
  if(typeof(write_to) == "raw") return(ans)
  writeBin(object = ans, con = write_to)
  return(invisible(write_to))
}

#' @title Tiff Writter
#' @name writetiff
#' @description
#' Writes TIFF from [0-1] Normalized Image.
#' @inheritParams writemulti
#' @param ... other arguments to be passed to \code{\link{writemulti}}.
#' @param what bits mode used to store image. Default is \code{"uint8"}. Allowed are \code{"uint8"}, \code{"int8"}, \code{"uint16"}, \code{"int16"}, \code{"uint32"}, \code{"int32"}, \code{"float"} and \code{"double"}.
#' @details it is aimed to simplify the use of \code{\link{writemulti}} by requiring minimal input from user.\cr
#' /!\ Note that image should range from 0 to 1.
#' /!\ Note that \code{as.rgb} that will be automatically filled with \code{TRUE} if \code{image} is 3D array with 3rd dimension being of length 3. So, if an RGB output is not desired, you should provide a 4D array.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   ## extract image values from cif file as [0-1] RGB 3D array
#'   img = ExtractImages_toMatrix(fileName = file_cif, mode = "rgb", objects = 3, selection = 3,
#'                                force_width = FALSE, size = c(0,0))
#'   ## create RGB tiff file
#'   file_tif = writetiff(image = img[[1]][[1]], write_to = tempfile(fileext = ".tif"))
#'   ## inspect exported file
#'   IFD = getIFD(file_tif)
#'   sapply(IFD[[1]]$tags, FUN = function(x) x)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @inherit writemulti return
#' @export
writetiff <- function(image, ..., what = "uint8") {
  dots = list(...)
  tmp = grepl("as.rgb", names(dots), fixed = TRUE)
  if(any(tmp)) dots = dots[!tmp]
  what = match.arg(arg = what, choices = c("uint8","int8","uint16","int16","uint32","int32","float","double"), several.ok = FALSE)
  vmax = switch(what, "uint8" = 255, "int" = 127, "uint16" = 65535, "int16" = 32767, "uint32" = 4294967295, "int32" = 2147483647, 1)
  d = dim(image)
  is.rgb = (length(d) == 3) && (d[3] == 3)
  if(substr(what,1,1) == "u") {
    img = image * vmax
    img[img < 0] <- 0
    do.call(args = c(dots, list(image = img, as.rgb = is.rgb, what = what)), what = writemulti) 
  } else {
    do.call(args = c(dots, list(image = image * vmax, as.rgb = is.rgb, what = what)), what = writemulti)
  }
}
