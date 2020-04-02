#' @title IFC_pops Indices Getter
#' @description
#' Retrieves indices of objects belonging to a population.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param pop a population name from 'obj'. Default is "".
#' If left as is or not found an error is thrown displaying all available population in 'obj'.
#' @examples
# #' \dontrun{
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#'   obj <- popsGetIndices(obj = daf, pop = names(daf$pops)[length(daf$pops)])
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
# #' }
#' @return An integer vector is returned
#' @export
popsGetIndices <- function(obj, pop = "") {
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  if(length(pop) != 1) stop("'pop' should be of length 1")
  N = names(obj$pops)
  if(!all(pop%in%N)) stop(paste0("pop:[",pop,"] was not found in DAF object, valid names are:\n", paste0(paste("-", N), collapse = "\n")))
  return(as.integer(which(obj$pops[[pop]][["obj"]])-1))
}
