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

#' @title XML Graphs Parser
#' @description 
#' Helper to parse graphs from XML document.
#' @param xml_doc an `xml_document` object.
#' @param caller a string, name of the calling function. Default is \code{character(0)}.
#' @param title_progress title displayed in progress bar. Default is \code{""}.
#' @param display_progress whether to display a progress bar. Default is \code{FALSE}.
#' @return a list of graphs.
#' @keywords internal
fromXML2_graphs <- function(xml_doc, caller = character(0), title_progress = "", display_progress = FALSE) {
  assert(xml_doc, cla="xml_document")
  plots=lapply(xml_attrs(xml_find_all(xml_doc, "//Graph")), FUN=function(x) as.list(x))
  if(length(plots)!=0) {
    plots_tmp=lapply(plots, FUN=function(plot) {
      pat=paste0("//Graph[@xlocation='",plot$xlocation,"'][@ylocation='",plot$ylocation,"']")
      sapply(c("Legend","BasePop","GraphRegion","ShownPop"), simplify=FALSE, FUN=function(i_subnode){
        lapply(xml_attrs(xml_find_all(xml_doc, paste(pat,i_subnode,sep="//"))), FUN=function(x) as.list(x))
      })
    })
    plots=mapply(plots, plots_tmp, FUN = append, SIMPLIFY = FALSE)
    plots_tmp=c("xlocation","ylocation","scaletype","xmin","xmax","ymin","ymax","axislabelsfontsize","axistickmarklabelsfontsize",
                "graphtitlefontsize","regionlabelsfontsize","bincount","histogramsmoothingfactor","xsize","ysize","splitterdistance","maxpoints")
    plots=lapply(plots, FUN=function(x) {plots_tmp = plots_tmp[plots_tmp %in% names(x)];replace(x, plots_tmp, lapply(x[plots_tmp], as.numeric))})
    plot_order=sapply(plots, FUN=function(i_plot) as.numeric(i_plot[c("xlocation", "ylocation")]))
    return(plots[order(unlist(plot_order[1,]),unlist(plot_order[2,]))])
  }
  return(list())
}
