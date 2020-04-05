#' @title IFC Files Checksum
#' @description
#' This function returns RIF/CIF checksum.
#' Checksum is the sum of img IFDs (Image Field Directory) offsets of objects 0, 1, 2, 3 and 4.
#' @param fileName path to file.
#' @param ... arguments to pass to \code{\link{checksumDAF}} or \code{\link{checksumXIF}}.
#' @details if fileName is a DAF file, then CIF checksum is computed from images values found in DAF.
#' @export
checksumIFC <- function(fileName, ...) {
  dots=list(...)
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf", "rif", "cif"))
  if(file_extension == "daf") return(checksumDAF(fileName = fileName, ...))
  return(checksumXIF(fileName = fileName, ...))
}