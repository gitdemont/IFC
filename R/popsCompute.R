#' @title IFC_pops Computation
#' @description
#' Function used to compute IFC_pops\cr
#' It requires pops, regions and features.
#' @param pops list of populations.
#' @param regions list of regions.
#' @param features data.frame of features.
#' @param pnt_in_poly_algorithm algorithm used to determine if object belongs to a polygon region or not. Default is 1.\cr
#' Note that for the moment only 1(Trigonometry) is available.
#' @param pnt_in_poly_epsilon epsilon to determine if object belongs to a polygon region or not. It only applies when algorithm is 1. Default is 1e-12.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param pb_title character string, giving the window title on the dialog box. Default is "".
#' @param bypass logical to avoid passing pops to buildPopulation() checking. Default is FALSE.
#' @source For pnt_in_poly_algorithm, Trigonometry, is an adaptation of Jeremy VanDerWal's code \url{http://github.com/jjvanderwal/SDMTools}
#' @export
popsCompute <- function(pops, regions, features, pnt_in_poly_algorithm = 1, pnt_in_poly_epsilon = 1e-12, display_progress = TRUE, pb_title = "", bypass = FALSE) {
  bypass = as.logical(bypass); assert(bypass, len=1, alw=c(TRUE,FALSE))
  if(!bypass) {
    pops = lapply(pops, FUN=function(p) do.call(what = "buildPopulation", args = p))
    class(pops) = "IFC_pops"
  } 
  pops = popsGetAffiliation(pops)
  pops = popsOrderNodes(pops)
  return(popsWithin(pops = pops, regions = regions, features = features, pnt_in_poly_algorithm = pnt_in_poly_algorithm, pnt_in_poly_epsilon = pnt_in_poly_epsilon, display_progress = display_progress, pb_title = pb_title))
}
