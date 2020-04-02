#' @title IFC_features_def XML Conversion
#' @description 
#' Helper to convert features definition (`IFC_features_def` object) to XML nodes.
#' @param features an `IFC_features_def` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a XML::xmlNode.
#' @keywords internal
toXML2_features_def =function(features_def, verbose = verbose) {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating features definition node")
  assert(features_def, cla = "IFC_features_def")
  if(length(features_def)==0) return(xml_new_node(name = "DefinedFeatures", text = ""))
  xml_new_node(name = "DefinedFeatures", .children = lapply(features_def, FUN=function(i_def) {
    xml_new_node(name = "UDF", attrs = i_def)
  }))
}
