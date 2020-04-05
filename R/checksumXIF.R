#' @title RIF/CIF File Checksum 
#' @description 
#' This function returns RIF/CIF checksum.
#' Checksum is the sum of img IFDs (Image Field Directory) offsets of objects 0, 1, 2, 3 and 4.
#' @param fileName path to file.
#' @param ... other arguments to be passed.
#' @keywords internal
checksumXIF = function(fileName, ...) {
  # TODO ask AMNIS how checksum is computed
  return(cpp_checksum(fileName))
}
