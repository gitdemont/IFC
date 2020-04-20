#' @title IFC_masks XML Conversion
#' @description 
#' Helper to convert masks (`IFC_masks` object) to XML nodes.
#' @param masks an `IFC_masks` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a xml_node.
#' @keywords internal
toXML2_masks = function(masks, verbose = FALSE) {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating masks node")
  assert(masks, cla = "IFC_masks")
  if(length(masks)==0) return(xml_new_node(name = "masks", text = ""))
  xml_new_node(name = "masks", .children = lapply(1:nrow(masks), FUN=function(i) {
    xml_new_node(name = "mask", attrs = masks[i, ])
  }))
}
