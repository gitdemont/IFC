################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2020 Yohann Demont                                             #
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

#' @title Remove Population from an IFC_data Object
#' @description
#' Removes populations from an already existing `IFC_data` object.
#' @inheritParams data_rm_graphs
#' @param pops a character vector of population names to remove within 'obj'. Note that "All" and "" are not allowed and will be excluded from 'pops' if present.
#' @return an `IFC_data` object or a list of elements impacted by removal depending on 'list_only' parameter.
#' @export
data_rm_pops <- function(obj, pops, list_only = TRUE, adjust_graph = TRUE, ...) {
  dots = list(...)
  assert(obj, cla = "IFC_data")
  assert(list_only, len = 1, alw = c(TRUE,FALSE))
  assert(adjust_graph, len = 1, alw = c(TRUE,FALSE))
  to_remove_pops = as.character(pops)
  assert(to_remove_pops, typ = "character")
  if(length(obj$pops) == 0) {
    warning("'obj' contains no population", immediate. = TRUE, call. = FALSE)
    if(list_only) {
      return(list(masks = character(),
                  features = character(),
                  regions = character(),
                  pops = character(),
                  graphs = integer()))
    } else {
      return(obj)
    }
  }
  
  # forbids removal of "All" and ""
  if(!list_only) {
    if(any(c("All", "") %in% to_remove_pops)) warning("\"All\" and \"\" are mandatory and can not be removed", immediate. = TRUE, call. = FALSE)
    to_remove_pops = setdiff(to_remove_pops, c("All", ""))
  }
  
  # removes duplicated inputs
  tmp = duplicated(to_remove_pops)
  if(any(tmp)) {
    warning(paste0("duplicated 'pops' automatically removed:\n", paste0(paste0("\t- ", to_remove_pops[tmp]), collapse = "\n")), immediate. = TRUE, call. = FALSE)
    to_remove_pops = to_remove_pops[!tmp]
  }
  
  # removes pops not in obj
  tmp = to_remove_pops %in% names(obj$pops)
  if(any(!tmp)) {
    warning(paste0("some 'pops' are not in 'obj$pops' and can't be removed:\n", paste0(paste0("\t- ", to_remove_pops[!tmp]), collapse = "\n")), immediate. = TRUE, call. = FALSE)
    to_remove_pops = to_remove_pops[tmp]
  }
  if(length(to_remove_pops) == 0) {
    warning("no population to remove in 'obj'", immediate. = TRUE, call. = FALSE)
    if(list_only) {
      return(list(masks = character(),
                  features = character(),
                  regions = character(),
                  pops = character(),
                  graphs = integer()))
    } else {
      return(obj)
    }
  }
  
  # search pops that depend on input pops
  if(length(obj$pops) > 0) for(i in 1:length(obj$pops)) {
    if(any(to_remove_pops %in% c(obj$pops[[i]]$base, obj$pops[[i]]$names))) to_remove_pops = c(to_remove_pops, obj$pops[[i]]$name)
  }
  to_remove_pops = unique(to_remove_pops)
  
  # search graphs that depend on input pops
  to_remove_graphs = integer()
  oo = obj
  oo$pops = obj$pops[setdiff(names(oo$pops), to_remove_pops)]
  for(i in seq_along(obj$graphs)) {
    if(length(adjustGraph(oo, obj$graphs[[i]], adjust_graph = adjust_graph)) == 0) to_remove_graphs = c(to_remove_graphs, i)
  }
  to_remove_graphs = unique(to_remove_graphs)
  
  # create list
  if(list_only) return(list(masks = character(),
                            features = character(),
                            regions = character(),
                            pops = to_remove_pops,
                            graphs = to_remove_graphs))
  
  pops_back = obj$pops
  obj$pops = list()
  obj = data_add_pops(obj, pops = pops_back[!(names(pops_back) %in% to_remove_pops)], ...)
  obj = data_rm_graphs(obj = obj, graphs = to_remove_graphs, list_only = list_only, adjust_graph = NA)
  return(obj)
}
