################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2021 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2021                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Jean-Pierre Marolleau, Loïc Garçon,       #
#                     CHU Amiens                                               #
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

#' @title ASSIST Database Extraction
#' @description
#' Retrieves ASSIST tests values stored within .cif / .rif files.
#' @param fileName path to file..
#' @param ... other arguments to be passed.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a rif file
#'   file_rif <- system.file("extdata", "example.rif", package = "IFCdata")
#'   ASSIST_db <- getASSIST(fileName = file_rif, from = "analysis")
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return a list of class `IFC_assist` of parsed ASSIST tests database.
#' @keywords internal
getASSIST <- function(fileName, ...) {
  dots = list(...)
  IFD = getIFD(fileName = fileName, offsets = "first", trunc_bytes = 8, force_trunc = FALSE, 
               verbose = FALSE, verbosity = 1, display_progress = FALSE, bypass = TRUE, ...)
  tmp_ass = read_xml(getFullTag(IFD = IFD, which = 1, "33064"), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  return(structure(lapply(as_list(xml_find_first(tmp_ass, "//ASSISTDb")), unlist), class = c("list", "IFC_assist")))
}
