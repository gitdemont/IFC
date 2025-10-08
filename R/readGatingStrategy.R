################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2021 Yohann Demont                                             #
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

################################################################################
#             functions described hereunder are experimental                   #
#              inputs and outputs may change in the future                     #
################################################################################

#' @title Gating Strategy File Reader
#' @description
#' Extracts Gating Strategy from files.
#' @param fileName path to file. It should be a .ast, .cif, .daf, .ist, .rif or .xml file.
#' @return A named list of class `IFC_gating`, whose members are:\cr
#' -spillover, a list of spillover matrices found,\cr
#' -graphs, a list of graphical elements found,\cr
#' -pops, a list describing populations found,\cr
#' -regions, a list describing how regions are defined.
#' @param ... other arguments to be passed.
#' @keywords internal
readGatingStrategy <- function(fileName, ...) {
  dots=list(...)
  if(missing(fileName)) stop("'fileName' can't be missing")
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf","cif","rif","xml","ast","ist"))
  if(!file.exists(fileName)) stop(paste("can't find",fileName,sep=" "))
  title_progress = basename(fileName)
  display_progress = FALSE
  if(length(dots[["display_progress"]]) != 0) display_progress = dots[["display_progress"]]
  assert(display_progress, len=1, alw=c(TRUE,FALSE))
  
  fileName = enc2native(normalizePath(fileName, winslash = "/", mustWork = FALSE))
  if(file_extension == "xml") return(readGatingML(fileName, ...))
  if(file_extension %in% c("daf", "ast")) {
    assay = switch(file_extension,
                   daf = "/Assay",
                   ast = "/AssayTemplate")
    toskip=cpp_scanFirst(fileName, charToRaw(paste0('<',assay,'>')), start = 0, end = 0)
    if(toskip==0) stop(paste0(fileName, "\ndoes not seem to be well formatted: <",assay,"> not found"))
    toskip = toskip + nchar(paste0("<",assay,">")) - 1
    tmp=read_xml(readBin(con = fileName, what = "raw", n = toskip), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  } else {
    if(file_extension %in% c("cif", "rif")) {
      IFD = getIFD(fileName = fileName, offsets = "first", trunc_bytes = 8, force_trunc = FALSE, verbose = FALSE, verbosity = 0, bypass = TRUE)
      tmp=read_xml(as_list(xml_find_first(read_xml(getFullTag(IFD = IFD, which = 1, tag = "33027", raw = TRUE),
                                          options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN")), "//Imaging//DafFile"))[[1]],
                   options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
    } else { # for ist
      tmp=read_xml(as_list(xml_find_first(read_xml(fileName, 
                                                   options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN")), "//Imaging//DafFile"))[[1]],
                   options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
    }
    assay = "/AssayTemplate"
  }

  ##### extracts graphs information
  plots=fromXML2_graphs(tmp)
  class(plots) <- "IFC_graphs"
  
  ##### extracts regions information
  regions=fromXML2_regions(tmp, caller = "readGatingStrategy")
  class(regions) <- "IFC_regions"
  
  ##### extracts populations information
  pops=fromXML2_pops(tmp, title_progress = title_progress, display_progress = display_progress)
  class(pops) <- "IFC_pops"
  
  ans = list("spillover"=list(), "graphs"=plots, "pops"=pops, "regions"=regions)
  attr(ans, "class") <- c("IFC_gating")
  return(ans)
}
