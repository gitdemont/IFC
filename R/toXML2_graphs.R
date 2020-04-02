#' @title IFC_graphs XML Conversion
#' @description 
#' Helper to convert graphs (`IFC_graphs` object) to XML nodes.
#' @param graphs an `IFC_graphs` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a xml_node.
#' @keywords internal
toXML2_graphs = function(graphs, verbose = FALSE) {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating graphs node")
  assert(graphs, cla = "IFC_graphs")
  if(length(graphs)==0) return(xml_new_node(name = "Displays", text = ""))
  graphs = lapply(graphs, FUN=function(i) {
    if(typeof(i) %in% c("integer","double")) { 
      return(cpp_num_to_string(i))
    } else {
      return(i)
    }
  })
  xml_new_node(name = "Displays", attrs = list(count = cpp_num_to_string(length(graphs)), layout="1", entriesPerRow="12", lightDarkMode="1"),
             .children = lapply(graphs, FUN=function(i_graph) {
               xml_new_node(name = "Graph", attrs = i_graph[!grepl("Legend|BasePop|GraphRegion|ShownPop", names(i_graph)) & sapply(i_graph, FUN=function(ele) length(ele)!=0)],
                          .children = c(lapply(i_graph[["Legend"]], FUN=function(i) xml_new_node(name = "Legend", attrs = i)),
                                        lapply(i_graph[["BasePop"]], FUN=function(i) xml_new_node(name = "BasePop", attrs = i)),
                                        lapply(i_graph[["GraphRegion"]], FUN=function(i) xml_new_node(name = "GraphRegion", attrs = i)),
                                        lapply(i_graph[["ShownPop"]], FUN=function(i) xml_new_node(name = "ShownPop", attrs = i))))
  }))
}
