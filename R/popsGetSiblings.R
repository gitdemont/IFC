#' @title IFC_pops Sibling Population Identification
#' @description
#' Gives names of graphical pops's siblings in a `IFC_data` object.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param pops graphical populations names to get siblings of.
#' @export
popsGetSiblings <- function(obj, pops) {
  if(missing(obj)) stop("'obj' can't be missing")
  if(missing(pops)) stop("'pops' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  if(is.null(pops)) stop("'pops' argument can't be NULL")
  N = names(obj$pops)
  if(!all(pops%in%N)) stop(paste0("pops was not found in DAF object, valid names are:\n", paste0(N, collapse=", ")))
  lapply(obj$pops[pops], FUN=function(p) {
    if(p$type!="G") return(p$name)
    map = sapply(obj$pops, FUN=function(m) {
      if(is.null(p$fy)) {
        return(c(m$base==p$base, ifelse(is.null(m$fx), FALSE, m$fx==p$fx), is.null(m$fy)))
      } else {
        return(c(m$base==p$base, ifelse(is.null(m$fx), FALSE, m$fx==p$fx), ifelse(is.null(m$fy), FALSE, m$fy==p$fy)))
      }
    })
    return(N[apply(map, 2, all)])
  })
}
