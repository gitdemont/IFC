#' @title IFC Files Generic Reader
#' @description
#' Extracts IFC data from IFC files no matter if they are DAF, RIF or CIF.
#' @param fileName path to file.
#' @param ... arguments to pass to \code{\link{ExtractFromDAF}} or \code{\link{ExtractFromXIF}}.
#' @details see \code{\link{ExtractFromDAF}} or \code{\link{ExtractFromXIF}}.
#' @examples
# #' \dontrun{
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a rif file, but you can also read daf or cif
#'   file_rif <- system.file("extdata", "example.rif", package = "IFCdata")
#'   rif <- readIFC(fileName = file_rif)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
# #' }
#' @export
readIFC <- function(fileName, ...) {
  dots=list(...)
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf", "rif", "cif"))
  if(file_extension == "daf") return(ExtractFromDAF(fileName = fileName, ...))
  return(ExtractFromXIF(fileName = fileName, ...))
}
