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

#' @title Object Colorizer
#' @description
#' Colorizes a [0,1] grayscale image.
#' @param mat a [0,1] numeric matrix.
#' @param color a character vector of color(s).
#' A length 1 \code{'color'} will map gray values of \code{'mat'} to hsv \code{'v'} value parameter of the conversion of \code{rgb2hsv(col2rgb('color'))}.
#' Whereas when \code{'color'} is of length > 1, it will be used to bin \code{'mat'} into \code{length('color')} intervals and map them to \code{'color'}.
#' @param as.raster a logical. Default is \code{FALSE}. 
#' @return depending on \code{'as.raster'}:\cr
#' -if \code{FALSE}, a 3D array where 3rd dimension is rgb.
#' -if \code{TRUE}, a `raster` matrix
#' -if \code{NA}, a color matrix
#' @keywords internal
objectColorize <- function(mat, color, as.raster = FALSE) {
  as.raster = as.logical(as.raster); assert(as.raster, len = 1)
  assert(color, typ = "character")
  if(length(color) == 1) {
    col = rgb2hsv(col2rgb(color))
    if(is.na(as.raster)) return(matrix(unclass(grDevices::as.raster(cpp_M_HSV2RGB(mat, h = col[1], s = col[2]))), nrow = nrow(mat), byrow = TRUE))
    if(as.raster) {
      return(grDevices::as.raster(cpp_M_HSV2RGB(mat, h = col[1], s = col[2])))
    } else {
      return(cpp_M_HSV2RGB(mat, h = col[1], s = col[2]))
    }
  }
  lev = level.colors(mat, at = seq(0, 1, length.out = 1 + length(color)), col.regions = color)
  if(is.na(as.raster)) return(matrix(lev, nrow = nrow(mat), byrow = FALSE))
  if(as.raster) return(structure(matrix(t(matrix(lev, nrow = nrow(mat), byrow = FALSE)), nrow = nrow(mat), byrow = FALSE), class = "raster"))
  aperm(array(col2rgb(lev) / 255, dim = c(3, nrow(mat), ncol(mat))), c(2,3,1))
}
