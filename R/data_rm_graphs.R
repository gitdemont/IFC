################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2022 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2020                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,             #
#                     Jean-Pierre Marolleau, Loïc Garçon,                      #
#                     INSERM, UPD, CHU Amiens                                  #
#                                                                              #
# DISCLAIMER:                                                                  #
# -You are using this package on your own risk!                                #
# -We do not guarantee privacy nor confidentiality.                            #
# -This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or #
# contributors be liable for any direct, indirect, incidental, special,        #
# exemplary, or consequential damages (including, but not limited to,          #
# procurement of substitute goods or services; loss of use, data, or profits;  #
# or business interruption) however caused and on any theory of liability,     #
# whether in contract, strict liability, or tort (including negligence or      #
# otherwise) arising in any way out of the use of this software, even if       #
# advised of the possibility of such damage.                                   #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with IFC. If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################

#' @title Remove Graph from an IFC_data Object
#' @description
#' Removes graphs from an already existing `IFC_data` object.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param graphs an integer vector of graph(s) to remove within 'obj'.
#' @param list_only whether to return a list of elements that will be impacted by the removal. Default is TRUE.
#' If FALSE then modified object will be returned.
#' @param adjust_graph whether to try to adjust graph(s) when possible. Default is TRUE.\cr
#' -TRUE, graph(s) will be kept if possible using only regions, pops it depends that can be found in 'obj',\cr
#' -FALSE, graph(s) will be kept only if all features, regions, pops it refers to are found in 'obj',\cr
#' -NA, graph(s) will be removed no matter if features, regions, pops it refers to are found in 'obj'.
#' @param ... Other arguments to be passed.
#' @return an `IFC_data` object or a list of elements impacted by removal depending on 'list_only' parameter.
#' @export
data_rm_graphs <- function(obj, graphs, list_only=TRUE, adjust_graph=TRUE, ...) {
  dots = list(...)
  assert(obj, cla = "IFC_data")
  assert(list_only, len = 1, alw = c(TRUE,FALSE))
  to_remove_graphs = as.integer(graphs)
  l_default = list(masks = character(),
                   features = character(),
                   regions = character(),
                   pops = character(),
                   graphs = integer())
  
  if(length(obj$graphs) == 0) {
    warning("'obj' contains no graph", immediate. = TRUE, call. = FALSE)
    if(list_only) {
      return(l_default)
    } else {
      return(obj)
    }
  }
  
  # removes duplicated inputs
  tmp = duplicated(to_remove_graphs)
  if(any(tmp)) {
    warning(paste0("duplicated 'graphs' automatically removed:\n", paste0(paste0("\t- ", to_remove_graphs[tmp]), collapse = "\n")), immediate. = TRUE, call. = FALSE)
    to_remove_graphs = to_remove_graphs[!tmp]
  }
  
  # removes graphs not in obj
  assert(to_remove_graphs, alw = seq_along(obj$graphs), fun = "warning")
  to_remove_graphs = to_remove_graphs[to_remove_graphs %in% seq_along(obj$graphs)]
  
  # creates list
  if(list_only) {
    l_default$graphs = to_remove_graphs
    return(l_default)
  }
  
  # otherwise return 
  if(length(to_remove_graphs) != 0) {
    G = obj$graphs
    for(i in to_remove_graphs) {
      G[[i]] = adjustGraph(obj = obj, graph =  obj$graphs[[i]], adjust_graph = adjust_graph)
    }
    obj$graphs = G[sapply(G, length) != 0]
    class(obj$graphs) = "IFC_graphs"
  }
  return(obj)
}
