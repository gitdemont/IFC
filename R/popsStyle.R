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

#' @title Population Style
#' @name get_pop_style
#' @description
#' Gets population style
#' @param x type of population.
#' @return a named UTF-8 character
#' @keywords internal
get_pop_style <- function(x) {
  structure(switch(
    x,
    "B"="\u25a0",    #"■"
    "T"="\u29bf",    #"⦿"
    "C"="\u2bba",    #"⮺"
    "line" ="\u21ff",#"⇿"
    "rect" ="\u25ad",#"▭"
    "poly" ="\u23e2",#"⏢"
    "oval" ="\u2b2f",#"⬯"
    "\ufeff"         #""
  ),names=x)
}

#' @title IFC_pops Style
#' @description
#' Helper to extract population style.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @return a named list of UTF-8 characters.
#' @keywords internal
popsStyle <- function(obj) {
  sapply(obj$pops, simplify = FALSE, FUN = function(p) {
    get_pop_style(ifelse(identical(p$type,"G"),ifelse(any(p$region %in% names(obj$regions)), obj$regions[[p$region]]$type,""), p$type))
  })
}
