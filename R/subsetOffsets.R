#' @title IFC_offset Subsetting
#' @description
#' Subsets IFC_offset
#' @param offsets object of class `IFC_offset` to subset.
#' @param objects integer vector, indices of objects to extract.
#' @param image_type image_type of desired offsets. Default is c("img", "msk"). Allowed are "img" and/or "msk".
#' @examples
# #' \dontrun{
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   ## extract offsets
#'   offsets <- getOffsets(fileName = file_cif)
#'   ## subset offsets of the 4 first "img" objects
#'   sub_offs <- subsetOffsets(offsets = offsets, objects = 0:3, image_type = "img")
#'   ## show subsetted offsets' structure
#'   str(sub_offs)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
# #' }
#' @return A class IFC_offset integer vector or empty list if objects are outside of offsets.
#' @export
subsetOffsets <- function(offsets, objects, image_type = c("img", "msk")) {
  if(missing(offsets)) stop("'offsets' can't be missing")
  if(missing(objects)) stop("'objects' can't be missing")
  assert(offsets, cla = "IFC_offset")
  assert(image_type, alw = c("img", "msk"))
  objects = as.integer(objects); objects = objects[(objects >= 0) & (objects < length(attr(offsets, "all")) / 2) & (is.finite(objects))]
  if(length(objects) == 0) {
    warning("subsetOffsets: No objects to subset, check the objects you provided.", immediate. = TRUE, call. = FALSE)
    foo = list()
  } else {
    if(all(c("img", "msk") %in% image_type)) {
      in_offsets = unlist(lapply(objects, FUN=function(x) c(2*x+1,2*(x+1))))
    } else {
      if("img" %in% image_type) {
        in_offsets = 2 * objects + 1
      } else {
        in_offsets = 2 * (objects + 1)
      }
    }
    foo = attr(offsets, "all")[in_offsets]
  }
  attr(foo, "all") = attr(offsets, "all")
  attr(foo, "first") = attr(offsets, "first")
  attr(foo, "fileName_image") = attr(offsets, "fileName_image")
  attr(foo, "checksum") = attr(offsets, "checksum")
  class(foo) = c("IFC_offset")
  return(foo)
}
