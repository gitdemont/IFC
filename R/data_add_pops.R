#' @title Pops Adder to IFC_data Object
#' @description
#' Adds populations to an already existing IFC_data object.
#' @param obj an IFC_data object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param pops a list of population(s) to add to obj. Each element of this list will be coerced by \code{\link{buildPopulation}}.
#' @param pnt_in_poly_algorithm algorithm used to determine if object belongs to a polygon region or not. Default is 1.\cr
#' Note that for the moment only 1(Trigonometry) is available.
#' @param pnt_in_poly_epsilon epsilon to determine if object belongs to a polygon region or not. It only applies when algorithm is 1. Default is 1e-12.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @details A warning will be thrown if a provided population is already existing in obj.\cr
#' In such a case this population will not be added to obj.\cr
#' If any input population is not well defined and can't be created then an error will occur.
#' @param ... Other arguments to be passed.
#' @source For pnt_in_poly_algorithm, Trigonometry, is an adaptation of Jeremy VanDerWal's code \url{http://github.com/jjvanderwal/SDMTools}
#' @examples
# #' \dontrun{
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#'   ## copy 1st population from existing daf
#'   pop <- daf$pops[[1]]
#'   if(length(pop) != 0) {
#'     pop_copy <- pop
#'     ## modify name, obj and type of copied population
#'     pop_copy$name <- paste0(pop_copy$name,"_copy")
#'     pop_copy$obj <- (which(pop_copy$obj)-1)[1]
#'     pop_copy$type <- "T"
#'     ## create new object with this new population
#'     dafnew <- data_add_pops(obj = daf, pops = list(pop_copy))
#'   }
#' } else {
#'   message(sprintf('Please type `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
# #' }
#' @return an IFC_data object with pops added.
#' @export
data_add_pops <- function(obj, pops, pnt_in_poly_algorithm = 1, pnt_in_poly_epsilon = 1e-12, display_progress = TRUE, ...) {
  dots = list(...)
  assert(obj, cla = "IFC_data")
  
  # try to coerce inputs to compatible daf format
  pops = lapply(pops, FUN=function(x) do.call(what=buildPopulation, args=x))
  names(pops) = sapply(pops, FUN=function(x) x$name)

  # removes duplicated inputs
  tmp = duplicated(names(pops))
  if(any(tmp)) {
    warning(paste0("duplicated pops automatically removed: ", names(pops)[tmp]), immediate. = TRUE, call. = FALSE)
    pops = pops[!tmp]
  }
  
  exported_pops = sapply(pops, FUN=function(pop) {
    if(pop$name%in%names(obj$pops)) {
      warning(paste0(pop$name, "\nnot exported: trying to export an already defined population"), immediate. = TRUE, call. = FALSE)
      return(FALSE)
    }
    return(TRUE)
  })
  exported_pops = pops[exported_pops]
  names(exported_pops) = sapply(exported_pops, FUN=function(x) x$name)
  
  exported_pops = popsCompute(pops = c(obj$pops, exported_pops), 
                              regions = obj$regions, 
                              features = obj$features, 
                              pnt_in_poly_algorithm = pnt_in_poly_algorithm, 
                              pnt_in_poly_epsilon = pnt_in_poly_epsilon, 
                              display_progress = display_progress, 
                              title_progress = basename(obj$fileName),
                              bypass = FALSE, ...)

  obj$pops = exported_pops
  obj_count = as.integer(obj$description$ID$objcount)
  if(nrow(obj$stats)!=0) {
    obj$stats = data.frame(stringsAsFactors = FALSE, check.rows = FALSE, check.names = FALSE, t(sapply(names(exported_pops), FUN=function(p) {
      count = sum(exported_pops[[p]]$obj)
      base = exported_pops[[p]]$base
      type = exported_pops[[p]]$type
      if(base=="") base = "All"
      parent = sum(exported_pops[[base]]$obj)
      c("type" = type, "parent" = base, "count" = count, "perc_parent" = count/parent*100, "perc_tot" = count/obj_count*100)
    })))
    obj$stats[,3] = as.numeric(obj$stats[,3])
    obj$stats[,4] = as.numeric(obj$stats[,4])
    obj$stats[,5] = as.numeric(obj$stats[,5])
  }
  return(obj)
}
