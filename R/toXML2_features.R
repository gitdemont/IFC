#' @title IFC_features XML Conversion
#' @description 
#' Helper to convert features (`IFC_features` object) to XML nodes.
#' @param features an `IFC_features` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @return a xml_node.
#' @keywords internal
toXML2_features = function(features, verbose = FALSE, display_progress = TRUE, title_progress = "") {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating features node")
  assert(features, cla = "IFC_features")
  if(length(features)==0) return(xml_new_node(name = "FeatureValues", text = ""))
  xml_new_node(name = "FeatureValues", .children = lapply(1:length(features), FUN=function(i) {
    xml_new_node(name = "UDFValues", attrs = list("fid" = i-1, "fv" = paste0(features[[i]], collapse = "|")))
  }))
}
