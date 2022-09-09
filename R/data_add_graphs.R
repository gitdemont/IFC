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

#' @title Add Graphs to IFC_data Object
#' @description
#' Adds graphs to an already existing `IFC_data` object.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param graphs a list of graph(s) to add to obj. Each element of this list will be coerced by \code{\link{buildGraph}}.
#' @param adjust_graph whether to try to adjust graph(s) when possible. Default is TRUE.\cr
#' -TRUE, graph(s) will be kept if possible using only regions, pops it depends that can be found in 'obj',\cr
#' -FALSE, graph(s) will be kept only if all features, regions, pops it refers to are found in 'obj',\cr
#' -NA, is not allowed and will throw an error.
#' @param ... Other arguments to be passed.
#' @return an IFC_data object with graphs added.
#' @export
data_add_graphs <- function(obj, graphs, adjust_graph=TRUE, ...) {
  assert(obj, cla = "IFC_data")
  adjust_graph = as.logical(adjust_graph); adjust_graph = assert(adjust_graph, len = 1, alw = c(TRUE,FALSE))
  L = length(obj$graphs)
  if(L == 0) obj$graphs = list()
  for(i in seq_along(graphs)) {
    obj$graphs = c(obj$graphs, graphs[i])
    obj$graphs = adjustGraph(obj = obj, selection = L + i, adjust_graph = adjust_graph)$graphs
    if(length(obj$graphs) != (L + i)) stop("can't add 'graphs[",i,"]' to 'obj'") 
  }
  class(obj$graphs) <- "IFC_graphs"
  return(obj)
}
