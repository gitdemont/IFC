#' @title Image Field Directory Full Tag Getter
#' @description
#' Retrieves full tag value from IFDs (Image Field Directory) extracted by \code{\link{getIFD}}.
#' @param fileName path to file..
#' @param ifd only one sub element of IFD data information extracted by \code{\link{getIFD}}.
#' @param tag numeric value or character value of the desired tag present ifd.
#' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
#' @details It is usefull to extract all information contained in a specific ifd tag since \code{\link{getIFD}} is designed to be run with argument trunc_bytes so as to only extract essential bytes to run faster and save memory.
#' Nonetheless, thanks to \code{\link{getFullTag}} users will still be able to get full extraction of specific tag.
#' @return If tag is of type 2, it returns a character string.\cr
#' If tag is of type 1, 6 or 7 it returns a raw vector.\cr
#' For other types, it will return mapped value corresponding to TIFF 6.0 specifications.
#' @export
getFullTag <- function(fileName, ifd, tag) {
  if(any(class(tag) == "character")) if(!(tag%in%names(ifd$tags))) stop("can't find 'tag' in 'ifd' tags names")
  if(any(class(tag) == "numeric")) if((tag<1) | (tag>length(ifd$tags))) stop("can't find 'tag' in 'ifd' tags")
  if(!("IFC_ifd"%in%class(ifd))) stop("'ifd' object is not of class IFC_ifd")
  tag = ifd$tags[[tag]]
  if(tag$typ == 2) if(tag$len == nchar(tag$map)) return(tag$map)
  if(tag$off) {
    #TODO add checksum
    toread = file(description = fileName, open = "rb")
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
