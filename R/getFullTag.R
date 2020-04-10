#' @title Image Field Directory Full Tag Getter
#' @description
#' Retrieves full tag value from IFDs (Image Field Directory) extracted by \code{\link{getIFD}}.
#' @param IFD an object of class `IFC_ifd_list` extracted by \code{\link{getIFD}}.
#' @param which scalar, integer (index) or the name of 'IFD' sub-element to extract 'tag' from. Default is 1 to extract 'tag' from the first member of 'IFD'.
#' @param tag scalar, integer (index) or the name of the IFD[[which]] of the desired 'tag'.
#' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
#' @details It may be usefull to extract all information contained in a specific 'tag' since \code{\link{getIFD}} is designed to be run with argument trunc_bytes so as to only extract essential bytes to run faster and save memory.
#' Nonetheless, thanks to \code{\link{getFullTag}} users will still be able to get full extraction of specific tag.
#' @return If 'tag' is of type 2, it returns a character string.\cr
#' If 'tag' is of type 1, 6 or 7 it returns a raw vector.\cr
#' For other types, it will return mapped value corresponding to TIFF 6.0 specifications.
#' @export
getFullTag <- function(IFD, which = 1, tag = "256") {
  if(!("IFC_ifd_list"%in%class(IFD))) stop("'IFD' object is not of class `IFC_ifd_list`")
  if(typeof(which) == "character") {
    if(length(which) != 1) stop("'which' should be of length 1")
    if(!(which %in% names(ifd$tags))) {
      warning("can't find 'which' in 'IFD'")
      return(NULL)
    }
  } else {
    which = na.omit(as.integer(which))
    if(length(which) != 1) stop("'which' should be of length 1")
    if((which < 1) || (which > length(IFD))) {
      warning("can't find 'which' in 'IFD'")
      return(NULL)
    }
  }
  ifd = IFD[[which]]
  if(typeof(tag) == "character") {
    if(length(tag) != 1) stop("'tag' should be of length 1")
    if(!(tag %in% names(ifd$tags))) {
      warning("can't find 'tag' in 'IFD[[which]]'")
      return(NULL)
    }
  } else {
    tag = na.omit(as.integer(tag))
    if(length(tag) != 1) stop("'tag' should be of length 1")
    if((tag < 1) || (tag > length(ifd$tags))) {
      warning("can't find 'tag' in IFD[[which]] tags")
      return(NULL)
    }
  }
  tag = ifd$tags[[tag]]
  if(tag$typ == 2) if(tag$len == nchar(tag$map)) return(tag$map)
  if(tag$off) {
    #TODO add checksum
    toread = file(description = attr(IFD, "fileName_image"), open = "rb")
    on.exit(close(toread))
    seek(toread, where = tag$val, origin = "start")
    if(tag$typ == 2) { 
      return(readChar(toread, nchars = tag$byt, useBytes = TRUE))
    } else {
      return(readBin(toread, n = tag$byt, what = "raw"))
    }
  } else {
    return(tag$val)
  }
}
