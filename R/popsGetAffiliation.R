#' @title IFC_pops Affiliation Finder
#' @description
#' Helper that extracts populations dependencies/affiliations.
#' @param pops list of populations
#' @keywords internal
popsGetAffiliation <- function(pops, operators = c("And","Or","Not","(",")")) {
  assert(pops, cla = "IFC_pops")
  K = class(pops)
  pops = sapply(pops, USE.NAMES = TRUE, simplify = FALSE, FUN = function(p) {
    if("C" %in% p$type) {
      p$split = splitn(definition = p$definition, all_names = names(pops), operators = operators)
      p$names = setdiff(p$split, operators)
    } else {
      p$names = ""
    }
    return(p)
  })
  class(pops) = c(K, "Affiliated")
  return(pops)
}
