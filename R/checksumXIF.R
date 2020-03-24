#' @title RIF/CIF File Checksum 
#' @description 
#' This function returns XIF checksum..
#' @param fileName path to file.
#' @keywords internal
checksumXIF = function(fileName) {
  # TODO ask AMNIS how checksum is computed
  return(cpp_checksum(fileName))
}
