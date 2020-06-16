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

#' @title Graph Retrieval from Graphical IFC_pops
#' @description
#' Retrieves the graph a graphical population originate from
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param pops names of graphical populations present in 'obj'. Note that they should be siblings.
#' @param vis2D when original graph is not an histogram, whether to display it as "scatter" or "density". Default is "density".
#' @param all_siblings whether to add all 'pop' siblings in the graph. Default is FALSE.
#' @return a list of parameters needed to build an IFC graph.
#' @keywords internal
popsRetrieveGraph = function(obj, pops, vis2D = "density", all_siblings = FALSE) {
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(missing(pops)) stop("'pops' can't be missing")
  pops = unique(pops); pops = as.character(pops); assert(pops, alw = names(obj$pops))
  if(!all(sapply(obj$pops[pops], FUN=function(p) p$type=="G"))) stop("'pops' should be of type graphical")
  # vis2D = as.character(vis2D); assert(vis2D, len=1, alw=c("scatter","density"))
  all_siblings = as.logical(all_siblings); assert(all_siblings, len = 1, alw =c(TRUE, FALSE))
  
  siblings = popsGetSiblings(obj, pops)
  are_siblings = all(sapply(siblings, FUN=function(s) s == siblings[[1]]))
  if(!are_siblings) stop("'pops' should be siblings")
  
  # initializes variables
  if(all_siblings) {
    pops = siblings[[1]]
  } else {
    pops = pops
  }
  P = obj$pops[pops]
  SUB = obj$pops[[obj$pops[[pops[1]]]$base]]$obj
  R = lapply(P, FUN=function(p) obj$regions[[p$region]])
  foo = list()
  
  # start rebuilding original graph
  foo$f1 = P[[1]]$fx
  foo$xlogrange = R[[1]]$xlogrange
  xran = range(obj$features[SUB, foo$f1], unlist(lapply(R, FUN=function(r) c(r$x, r$cx))), na.rm = TRUE)
  if(foo$xlogrange == "P") {
    xran = xran + diff(xran) * c(-0.1,0.1)
  } else {
    xran = smoothLinLog(xran, hyper = as.numeric(foo$xlogrange))
    xran = xran + diff(xran) * c(-0.1,0.1)
    xran = inv_smoothLinLog(xran, hyper = as.numeric(foo$xlogrange))
  }
  foo$xmin = xran[1]
  foo$xmax = xran[2]
  foo$ShownPop = list()
  foo$title = P[[1]]$base
  if(length(P[[1]]$fy) == 0) {
    foo$type = "histogram"
  } else {
    foo$f2 = P[[1]]$fy
    foo$ylogrange = R[[1]]$ylogrange
    yran = range(obj$features[SUB, foo$f2], unlist(lapply(R, FUN=function(r) c(r$y,r$cy))), na.rm = TRUE)
    if(foo$ylogrange == "P") {
      yran = yran + diff(yran) * c(-0.1,0.1)
    } else {
      yran = smoothLinLog(yran, hyper = as.numeric(foo$ylogrange))
      yran = yran + diff(yran) * c(-0.1,0.1)
      yran = inv_smoothLinLog(yran, hyper = as.numeric(foo$ylogrange))
    }
    foo$ymin = yran[1]
    foo$ymax = yran[2]
    foo$type = vis2D
  }
  foo$BasePop = list(list(name = P[[1]]$base))
  foo$GraphRegion = lapply(R, FUN=function(r) list("name"=r$label))
  return(foo)
}
