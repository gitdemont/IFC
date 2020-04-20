#' @title IFC_regions XML conversion
#' @description 
#' Helper to convert regions (`IFC_regions` object) to XML nodes.
#' @param regions an `IFC_regions` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a xml_node.
#' @keywords internal
toXML2_regions = function(regions, verbose = FALSE) {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating regions node")
  assert(regions, cla = "IFC_regions")
  if(length(regions)==0) return(xml_new_node(name = "Regions", text = ""))
  xml_new_node(name = "Regions", .children = lapply(regions, FUN=function(i_reg) {
    if(i_reg$color=="Cyan4") i_reg$color <- "Teal"
    if(i_reg$lightcolor=="Cyan4") i_reg$lightcolor <- "Teal"
    if(i_reg$color=="Green4") i_reg$color <- "Green"
    if(i_reg$lightcolor=="Green4") i_reg$lightcolor <- "Green"
    if(i_reg$color=="Chartreuse") i_reg$color <- "Lime"
    if(i_reg$lightcolor=="Chartreuse") i_reg$lightcolor <- "Lime"
    xml_new_node(name = "Region",
               attrs = i_reg[!grepl("^x$|^y$", names(i_reg))],
               .children = lapply(1:length(i_reg[["x"]]), FUN = function(i_coord) {
                 xml_new_node(name = "axy", attrs = list(x = cpp_num_to_string(i_reg[["x"]][i_coord]), y = cpp_num_to_string(i_reg[["y"]][i_coord])))
               }))
  }))
}
