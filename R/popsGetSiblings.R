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

#' @title IFC_pops Sibling Population Identification
#' @description
#' Gives names of graphical pops's siblings in a `IFC_data` object.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param pops graphical populations names to get siblings of.
#' @return names of population siblings.
#' @export
popsGetSiblings <- function(obj, pops) {
  if(missing(obj)) stop("'obj' can't be missing")
  if(missing(pops)) stop("'pops' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  if(is.null(pops)) stop("'pops' argument can't be NULL")
  N = names(obj$pops)
  if(!all(pops%in%N)) stop(paste0("pops was not found in 'obj', valid names are:\n", paste0(N, collapse=", ")))
  lapply(obj$pops[pops], FUN=function(p) {
    if(p$type!="G") return(p$name)
    map = sapply(obj$pops, FUN=function(m) {
      if(is.null(p$fy)) {
        return(c(m$base==p$base, ifelse(is.null(m$fx), FALSE, m$fx==p$fx), is.null(m$fy)))
      } else {
        return(c(m$base==p$base, ifelse(is.null(m$fx), FALSE, m$fx==p$fx), ifelse(is.null(m$fy), FALSE, m$fy==p$fy)))
      }
    })
    return(N[apply(map, 2, all)])
  })
}
