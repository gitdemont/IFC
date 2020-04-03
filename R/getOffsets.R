#' @title RIF/CIF File Image Field Directories Offsets Extraction
#' @description
#' Extracts offsets of the IFDs (Image Field Directories) within a XIF file.
#' @param fileName path to file.
#' @param fast whether to fast extract objects or not. Default is TRUE.\cr
#' Meaning that offsets will be extracting expecting that raw object are stored in ascending order.\cr
#' A message will be thrown since fast extraction method does not ensure correct mapping between objects and offsets.\cr
#' If set to FALSE, all object_ids will be scanned from 'fileName' to ensure extraction of desired offsets.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
#' @details Offsets are byte positions of IFDs found within RIF or CIF file. For more details see TIFF specifications.
#' @examples
# #' \dontrun{
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   system.time(offsets_fast <- getOffsets(fileName = file_cif, fast = TRUE))
#'   system.time(offsets_slow <- getOffsets(fileName = file_cif, fast = FALSE))
#'   identical(offsets_fast, offsets_slow)   
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
# #' }
#' @return an integer vector of class `IFC_offset` of IFD offsets found in XIF file.
#' If no offsets is found an error is thrown.
#' @export
getOffsets <- function(fileName, fast = TRUE, display_progress = TRUE, verbose = FALSE) {
  fast = as.logical(fast); assert(fast, len = 1, alw = c(TRUE, FALSE))
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  fileName = normalizePath(fileName, winslash = "/", mustWork = FALSE)
  obj_count = as.integer(getInfo(fileName, warn = FALSE, force_default = TRUE, display_progress = FALSE)$objcount)
  if(fast) {
    offsets = cpp_getoffsets_noid(fileName, obj_count = obj_count, display_progress = display_progress, verbose = verbose)
    offsets_1 = offsets[1]
    offsets = offsets[-1]
    if(length(offsets) != obj_count * 2) stop("Number of offsets found is different from expected object count.")
    message("Offsets were extracted from XIF file with fast method.\nCorrect mapping between offsets and objects ids is not guaranteed.")
  } else {
    offsets = as.data.frame(do.call(what = "cbind", args = cpp_getoffsets_wid(fileName, obj_count = obj_count, display_progress = display_progress, verbose = verbose)), stringsAsFactors = FALSE)
    offsets_i = offsets[offsets$TYPE == 2, ]
    offsets_m = offsets[offsets$TYPE == 3, ]
    
    ORD_i = order(offsets_i$OBJECT_ID)
    offsets_img = offsets_i$OFFSET[ORD_i]
    ORD_m = order(offsets_m$OBJECT_ID)
    offsets_msk = offsets_m$OFFSET[ORD_m]
    
    offsets_1 = offsets[offsets$TYPE == 1, ]
    offsets_1 = offsets_1$OFFSET
    
    if(length(offsets_img) != length(offsets_msk)) stop("Offsets contain different numbers of 'img' and 'msk'.")
    if(length(offsets_img) != obj_count) stop("Number of offsets found is different from expected object count.")
    offsets = as.integer(apply(cbind(offsets_img, offsets_msk), 1, FUN=function(i) i))
  }
  if(length(offsets) != 0) {
    N = nchar(sprintf("%1.f",abs(obj_count-1)))
    names(offsets) = c(paste0(c("img_", "msk_"), rep(sprintf(paste0("%0",N,".f"), 0:(obj_count-1)), each = 2)))
    attr(offsets, "all") = offsets
    attr(offsets, "fileName_image") = fileName
    attr(offsets, "checksum") = checksumXIF(fileName)
    attr(offsets, "class") = c("IFC_offset")
    attr(offsets, "first") = offsets_1 
    return(offsets)
  }
  stop(paste0("No IFD offsets found in\n", fileName))
}
