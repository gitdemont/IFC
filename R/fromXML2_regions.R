################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2025 Yohann Demont                                             #
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

#' @title XML Regions Parser
#' @description 
#' Helper to parse regions from XML document.
#' @param xml_doc an `xml_document` object.
#' @param caller a string, name of the calling function. Default is \code{character(0)}.
#' @param title_progress title displayed in progress bar. Default is \code{""}.
#' @param display_progress whether to display a progress bar. Default is \code{FALSE}.
#' @return a list of regions.
#' @keywords internal
fromXML2_regions <- function(xml_doc, caller = character(0), title_progress = "", display_progress = FALSE) {
  assert(xml_doc, cla="xml_document")
  regions=lapply(xml_attrs(xml_find_all(xml_doc, "//Region")), FUN=function(x) as.list(x))
  if(length(regions) != 0) {
    names(regions)=lapply(regions, FUN=function(x) x$label)
    regions_tmp=c("cx","cy")
    regions=lapply(regions, FUN=function(x) {replace(x, regions_tmp, lapply(x[regions_tmp], as.numeric))})
    regions_tmp=lapply(regions, FUN=function(i_region) {
      pat=paste0("//Region[@label='",i_region$label,"']//axy")
      axy=do.call(cbind, args = xml_attrs(xml_find_all(xml_doc, pat)))
      x=as.numeric(axy["x",]); y = as.numeric(axy["y",])
      nas = is.na(x) | is.na(y)
      if(any(nas)) warning(paste0(caller,"NAs in region['",collapse = ": "),i_region$label,"'] have been removed", call. = FALSE, immediate. = TRUE)
      x=x[!nas]; y=y[!nas]
      if(length(x) < ifelse(identical(i_region$type, "poly"), 1, 2)) stop("invalid vertices length for region['",i_region$label,"']")
      list(x=x,y=y)
    })
    regions=mapply(FUN = append, regions, regions_tmp, SIMPLIFY = FALSE)
    rm(regions_tmp)
    ##### changes unknown color names in regions and retrieves sync attribute if any
    for(i in seq_along(regions)) {
      sync = regions[[i]]$sync
      regions[[i]] = regions[[i]][setdiff(names(regions[[i]]), "sync")]
      attr(regions[[i]], "sync") = sync
      regions[[i]]$color = map_color(regions[[i]]$color)
      regions[[i]]$lightcolor = map_color(regions[[i]]$lightcolor)
    }
    return(sync_validate(regions, quietly = FALSE, caller = caller))
  }
  return(list())
}
