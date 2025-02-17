################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2020 Yohann Demont                                             #
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

#' @title Image Field Directory Full Tag Retrieval
#' @description
#' Retrieves full tag value from IFDs (Image Field Directory) extracted by \code{\link{getIFD}}.
#' @param IFD an object of class `IFC_ifd_list` extracted by \code{\link{getIFD}}.
#' @param which scalar, integer (index) or the name of 'IFD' sub-element to extract 'tag' from. Default is 1 to extract 'tag' from the first member of 'IFD'.
#' @param tag scalar, integer (index) or the name of the IFD[[which]] of the desired 'tag'.
#' @param raw whether to return tag as a raw vector. Default is FALSE.
#' @source TIFF 6.0 specifications archived from web \url{https://web.archive.org/web/20211209104854/https://www.adobe.io/open/standards/TIFF.html}
#' @details It may be usefull to extract all information contained in a specific 'tag' since \code{\link{getIFD}} is designed to be run with argument trunc_bytes so as to only extract essential bytes to run faster and save memory.
#' Nonetheless, thanks to \code{\link{getFullTag}} users will still be able to get full extraction of specific tag.
#' @return the full value of the corresponding IFD tag.
#' @export
getFullTag <- function(IFD, which = 1, tag = "256", raw = FALSE) {
  if(!("IFC_ifd_list"%in%class(IFD))) stop("'IFD' object is not of class `IFC_ifd_list`")
  assert(raw, len=1, alw=c(TRUE, FALSE))
  endian = cpp_checkTIFF(enc2native(attr(x = IFD, which = "fileName_image")))
  if(typeof(which) == "character") {
    if(length(which) != 1) stop("'which' should be of length 1")
    if(!(which %in% names(ifd$tags))) {
      warning("can't find 'which'[",which,"] in 'IFD'")
      return(NULL)
    }
  } else {
    which = na.omit(as.integer(which))
    if(length(which) != 1) stop("'which' should be of length 1")
    if((which < 1) || (which > length(IFD))) {
      warning("can't find 'which'[",which,"] in 'IFD'")
      return(NULL)
    }
  }
  ifd = IFD[[which]]
  if(typeof(tag) == "character") {
    if(length(tag) != 1) stop("'tag' should be of length 1")
    if(!(tag %in% names(ifd$tags))) {
      warning("can't find 'tag'[",tag,"] in 'IFD[[",which,"]]'")
      return(NULL)
    }
  } else {
    tag = na.omit(as.integer(tag))
    if(length(tag) != 1) stop("'tag' should be of length 1")
    if((tag < 1) || (tag > length(ifd$tags))) {
      warning("can't find 'tag'[",tag,"] in IFD[[",which,"]] tags")
      return(NULL)
    }
  }
  if(ifd$tags[[tag]]$byt > 4) {
    toread = file(description = attr(IFD, "fileName_image"), open = "rb")
    on.exit(close(toread))
    seek(toread, where = ifd$tags[[tag]]$val, origin = "start")
    if(raw) return(readBin(toread, n = ifd$tags[[tag]]$byt, what = "raw", signed = FALSE))
    switch(ifd$tags[[tag]]$typ,
           { # 1 BYTE, 1 Byte
             return(readBin(toread, n = ifd$tags[[tag]]$byt, what = "raw", signed = FALSE))
           },
           { # 2 ASCII, 1 Byte
             if(length(ifd$tags[[tag]]$byt) == 0) return(character(0))
             tmpcon <- rawConnection(readBin(toread, n = ifd$tags[[tag]]$byt, what = "raw", signed = FALSE))
             on.exit(close(tmpcon), add = TRUE)
             return(paste(suppressWarnings(readLines(tmpcon, skipNul = FALSE, encoding = "UTF-8")), sep = "", collapse = "\n"))
           },
           { # 3 SHORT 2 bytes
             return(readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 2, signed = FALSE, endian = endian))
           },
           { # 4 LONG, 4 bytes
             return(readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, endian = endian))
           },
           { # 5 RATIONAL = 2 LONG, 1st numerator, 2nd denominator
             if(length(ifd$tags[[tag]]$byt) == 0) return(numeric(0))
             foo = readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, endian = endian) 
             if(ifd$tags[[tag]]$len %% 2) {
               warning(paste0("IFD 'tag'[",tag,"] is of 'typ'[5] is expected to be even but number of bytes is odd.\neven was not divided by odd in returned value."))
               return(foo)
             } 
             odd = seq(from = 1, to = ifd$tags[[tag]]$len, by = 2)
             return(foo[odd] / foo[-odd])
           },
           { # 6 SBYTE, 1 Byte
             return(readBin(toread, n = ifd$tags[[tag]]$byt, what = "raw", signed = TRUE))
           },
           { # 7 UNDEFINED, 1 Byte
             return(readBin(toread, n = ifd$tags[[tag]]$byt, what = "raw"))
           },
           { # 8 SSHORT, 2 bytes
             return(readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, endian = endian))
           },
           { # 9 SLONG, 4 bytes
             return(readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, endian = endian))
           },
           { # 10 SRATIONAL, 2 SLONG, 1st numerator, 2nd denominator
             if(length(ifd$tags[[tag]]$byt) == 0) return(numeric(0))
             foo = readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, endian = endian) 
             if(ifd$tags[[tag]]$len %% 2) {
               warning(paste0("IFD 'tag'[",tag,"] is of 'typ'[10] is expected to be even but number of bytes is odd.\neven was not divided by odd in returned value."))
               return(foo)
             } 
             odd = seq(from = 1, to = ifd$tags[[tag]]$len, by = 2)
             return(foo[odd] / foo[-odd])
           },
           { # 11 FLOAT, 4 bytes
             return(readBin(toread, n = ifd$tags[[tag]]$len, what = "double", size = 4, endian = endian))
           },
           { # 12 DOUBLE, 8 bytes
             return(readBin(toread, n = ifd$tags[[tag]]$len, what = "double", size = 8, endian = endian))
           })
  } else {
    if(raw) {
      if(ifd$tags[[tag]]$byt == 0) return(raw())
      v = c(09,10,11,12)[1:ifd$tags[[tag]]$byt]
      if(endian != .Platform$endian) v = rev(v)
      return(ifd$tags[[tag]]$raw[v])
    } else {
      return(ifd$tags[[tag]]$map) 
    }
  }
}
