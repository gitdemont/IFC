#' @title Image Field Directory Full Tag Getter
#' @description
#' Retrieves full tag value from IFDs (Image Field Directory) extracted by \code{\link{getIFD}}.
#' @param IFD an object of class `IFC_ifd_list` extracted by \code{\link{getIFD}}.
#' @param which scalar, integer (index) or the name of 'IFD' sub-element to extract 'tag' from. Default is 1 to extract 'tag' from the first member of 'IFD'.
#' @param tag scalar, integer (index) or the name of the IFD[[which]] of the desired 'tag'.
#' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
#' @details It may be usefull to extract all information contained in a specific 'tag' since \code{\link{getIFD}} is designed to be run with argument trunc_bytes so as to only extract essential bytes to run faster and save memory.
#' Nonetheless, thanks to \code{\link{getFullTag}} users will still be able to get full extraction of specific tag.
#' @export
getFullTag <- function(IFD, which = 1, tag = "256") {
  if(!("IFC_ifd_list"%in%class(IFD))) stop("'IFD' object is not of class `IFC_ifd_list`")
  endian = cpp_checkTIFF(attr(x = IFD, which = "fileName_image"))
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
  toread = file(description = attr(IFD, "fileName_image"), open = "rb")
  on.exit(close(toread))
  seek(toread, where = ifd$tags[[tag]]$val, origin = "start")
  switch(ifd$tags[[tag]]$typ,
         { return(readBin(toread, n = ifd$tags[[tag]]$byt, what = "raw")) # 1 BYTE
         },
         { return(readChar(toread, nchars = ifd$tags[[tag]]$byt, useBytes = TRUE)) # 2 ASCII
         },
         { return(readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 2, signed = FALSE, endian = endian)) # 3 SHORT 2 bytes
         },
         { return(readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, signed = FALSE, endian = endian)) # 4 LONG, 4 bytes
         },
         { # 5 RATIONAL = 2 LONG, 1st numerator, 2nd denominator
           if(length(ifd$tags[[tag]]$byt) == 0) return(numeric(0))
           foo = readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, signed = FALSE, endian = endian) 
           if(ifd$tags[[tag]]$len %% 2) {
             warning(paste0("IFD 'tag'[",tag,"] is of 'typ'[5] is expected to be even but number of bytes is odd.\neven was not divided by odd in returned value."))
             return(foo)
           } 
           odd = seq(from = 1, to = ifd$tags[[tag]]$len, by = 2)
           return(foo[odd] / foo[-odd])
         },
         { return(readBin(toread, n = ifd$tags[[tag]]$byt, what = "raw")) # 6 SBYTE
         },
         { return(readBin(toread, n = ifd$tags[[tag]]$byt, what = "raw")) # 7 UNDEFINED, 1 Byte
         },
         { return(readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, signed = TRUE, endian = endian)) # 8 SSHORT, 2 bytes
         },
         { return(readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, signed = TRUE, endian = endian)) # 9 SLONG, 4 bytes
         },
         { # 10 SRATIONAL, 2 SLONG, 1st numerator, 2nd denominator
           if(length(ifd$tags[[tag]]$byt) == 0) return(numeric(0))
           foo = readBin(toread, n = ifd$tags[[tag]]$len, what = "integer", size = 4, signed = TRUE, endian = endian) 
           if(ifd$tags[[tag]]$len %% 2) {
             warning(paste0("IFD 'tag'[",tag,"] is of 'typ'[10] is expected to be even but number of bytes is odd.\neven was not divided by odd in returned value."))
             return(foo)
           } 
           odd = seq(from = 1, to = ifd$tags[[tag]]$len, by = 2)
           return(foo[odd] / foo[-odd])
         },
         { return(readBin(toread, n = ifd$tags[[tag]]$len, what = "double", size = 4, endian = endian)) # 11 FLOAT, 4 bytes
         },
         { return(readBin(toread, n = ifd$tags[[tag]]$len, what = "double", size = 8, endian = endian)) # 12 DOUBLE, 8 bytes
         })
}
