#' @title IFC_pops Levels Dependency Determination
#' @description
#' Helper that extract population levels.
#' @param pops list of populations.
#' @keywords internal
popsGetLevels <- function(pops) {
  assert(pops, cla = c("IFC_pops","Affiliated","Ordered"))
  l = length(pops)
  lev = rep(1, l)
  names(lev) = names(pops)
  if(l > 1) { 
    for(i in 2:l) {
      pop = pops[[i]]
      lev[i] = lev[i] + max(lev[setdiff(c(pop$base, pop$names), "")])
    }
  }
  return(lev)
}
