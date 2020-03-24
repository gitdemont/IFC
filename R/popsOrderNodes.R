#' @title IFC_pops Population Level Ordering
#' @description
#' Helper that sort populations so that populations that depend from other are placed after.
#' @param pops list of populations.
#' @keywords internal
popsOrderNodes <- function(pops) {
  assert(pops, cla = c("IFC_pops","Affiliated"))
  K = class(pops)
  i=1; l=length(pops)
  while(i<l) {
    pop=pops[[i]]
    index=setdiff(c(pop$base,pop$names),"")
    index=unlist(lapply(index, function(x) which(x==names(pops))))
    index=index[index>i]
    if(length(index)!=0) {
      pops = c(pops[index],pops[setdiff(1:l,index)])
      i=1
    } else {
      i=i+1
    }
  }
  class(pops) = c(K, "Ordered")
  return(pops)
}
