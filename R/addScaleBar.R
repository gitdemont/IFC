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

#' @title Image Scale Bar Incrustation
#' @name addScaleBar
#' @description Adds scale bar to image
#' @param image a [0,1] image.
#' @param size positive integer. Scale's bar size in micro-meter. Default is \code{7}.
#' @param style a character string. Scale's bar style, either \code{"dash"} or \code{"line"}. Default is \code{"dash"}.
#' @param color a character string. Color of the scale. Default is \code{"white"} .
#' @param pix positive numeric. Size of one pixel in micro-meter. Default is \code{0.3} for 60x magnification, use \code{0.5} for 40x and \code{1} for 20x.
#' @param xoff x offset in image to draw scale, starting from bottom left corner. Default is \code{0}.
#' @param yoff y offset in image to draw scale, starting from bottom left corner. Default is \code{0}.
#' @return an image with scale added to the bottom left corner.
#' @keywords internal
addScaleBar <- function(image, size=7, style=c("dash","line")[1], color="white", pix=c(0.3,0.5,1)[1], xoff=0, yoff=0) {
  # several checks
  size = na.omit(as.integer(size)); size = size[size>0]
  assert(size, len = 1, typ = "integer")
  size = as.character(size)
  style = na.omit(as.character(style))
  assert(style, len = 1, alw = c("dash","line"))
  color = na.omit(as.character(color))
  assert(color, len = 1, typ = "character")
  pix = na.omit(as.numeric(pix)); pix = pix[pix>=0]
  assert(pix, len = 1, typ = "numeric")
  xoff = as.integer(na.omit(xoff))
  assert(xoff, len = 1, typ = "integer")
  yoff = as.integer(na.omit(yoff))
  assert(yoff, len = 1, typ = "integer") 
  
  # add text
  d = dim(image)
  bar_w = ceiling(as.numeric(size)/pix)
  if(style == "dash") bar_scheme = c(1,1,1,0,0,0)
  if(style == "line") bar_scheme = c(1,1,1,1,1,1)
  bar_msk = rep_len(bar_scheme,bar_w)
  bar_msk = t(sapply(1:4, FUN=function(x) bar_msk))
  bar_msk[,c(1,bar_w)] <- 1
  bar_msk[c(1,4),] <- 1
  bar_img = bar_msk
  bar_img = objectColorize(bar_img,color)
  ret = array(sapply(1:d[3], FUN=function(x) cpp_mark2(A = image[,,x], B = bar_img[,,x], mask = bar_msk, xoff = 2 + xoff, yoff = d[1] - 12 + 4 - yoff, invert = FALSE)),dim = d)
  yoff = yoff - 3
  xoff = bar_w + 2 + xoff
  ret = addText2(image = ret, text = size, color = color, xoff = xoff, yoff = yoff, anchor = c(0, 0), vjust = "B", hjust = "L")
  xoff = xoff + dim(texttomatrix(size))[2]
  ret = addText2(image = ret, text = "|", color = color, xoff = xoff, yoff = yoff + 2, anchor = c(0, 0), vjust = "B", hjust = "L")
  xoff = xoff + dim(texttomatrix("|"))[2]
  return(addText2(image = ret, text = "m", color = color, xoff = xoff, yoff = yoff, anchor = c(0, 0), vjust = "B", hjust = "L"))
}
