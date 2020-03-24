#' @title IFC Files Generic Writer
#' @description
#' Writes IFC data to DAF or subsets RIF/CIF Files.
#' @param fileName path to file.
#' @param ... arguments to pass to \code{\link{ExportToDAF}} or \code{\link{ExportToXIF}}.
#' @details see \code{\link{ExportToDAF}} or \code{\link{ExportToXIF}}.
#' @examples
# #' \dontrun{
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   tmp <- tempdir(check = TRUE)
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   ## create a tagged population named test with 1st object
#'   pop <- buildPopulation(name = "test", type = "T", obj = 0)
#'   writeIFC(file_daf, export_to = paste0(tmp, "\\test_write.daf"),
#'            overwrite = TRUE, pops = list(pop))
#'   ## use a rif file, but you can also use a cif
#'   file_rif <- system.file("extdata", "example.rif", package = "IFCdata")
#'   writeIFC(fileName = file_rif, export_to = paste0(tmp, "\\test_write.rif"), 
#'            overwrite = TRUE, objects = 0)
#' } else {
#'   message(sprintf('Please type `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
# #' }
#' @export
writeIFC <- function(fileName, ...) {
  dots=list(...)
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf", "rif", "cif"))
  if(file_extension == "daf") return(ExportToDAF(fileName = fileName, ...))
  return(ExportToXIF(fileName = fileName, ...))
}
