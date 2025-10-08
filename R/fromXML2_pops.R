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

#' @title XML Populations Parser
#' @description 
#' Helper to parse populations from XML document.
#' @param xml_doc an `xml_document` object.
#' @param caller a string, name of the calling function. Default is \code{character(0)}.
#' @param title_progress title displayed in progress bar. Default is \code{""}.
#' @param display_progress whether to display a progress bar. Default is \code{FALSE}.
#' @return a list of populations.
#' @keywords internal
fromXML2_pops <- function(xml_doc, caller = character(0), title_progress = "", display_progress = FALSE) {
  assert(xml_doc, cla="xml_document")
  pops=lapply(xml_attrs(xml_find_all(xml_doc, "//Pop")), FUN=function(x) as.list(x))
  if(length(pops)>0) {
    names(pops)=lapply(pops, FUN=function(x) x$name)
    if(display_progress) {
      pb_pops = newPB(min = 0, max = length(pops), initial = 0, style = 3)
      tryCatch({
        pops_=lapply(1:length(pops), FUN=function(i_pop) {
          setPB(pb_pops, value = i_pop, title = title_progress, label = "extracting tagged population objects")
          pat=paste0("//Pop[@name='",pops[[i_pop]]$name,"']//ob")
          list(obj=as.integer(unlist(xml_attrs(xml_find_all(xml_doc, pat)))))
        })
      }, error = function(e) {
        stop(e$message)
      }, finally = endPB(pb_pops))
    } else {
      pops_=lapply(1:length(pops), FUN=function(i_pop) {
        pat=paste0("//Pop[@name='",pops[[i_pop]]$name,"']//ob")
        list(obj=as.integer(unlist(xml_attrs(xml_find_all(xml_doc, pat)))))
      })
    }
    return(mapply(FUN = append, pops, pops_, SIMPLIFY = FALSE))
  }
  return(list())
}
