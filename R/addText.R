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

#' @title Image Text Incrustation
#' @name addText
#' @description Adds Text to image.
#' @param image a [0,1] image.
#' @param text a character string.
#' @param color a character string. color of the text.
#' @param xoff positive integer. x offset in image to start writing text. Default is \code{0}.
#' @param yoff positive integer. y offset in image to start writing text. Default is \code{0}.
#' @param corner a character string. where to position text in the image. Allowed are \code{"TL"}, \code{"TR"}, \code{"BL"}, \code{"BR"}, for top-left, top-right, bottom-left, bottom-right, respectively.
#' @details One-lined text will be added so has to be fully contained within image and anchored at desired corner plus x and y offset from it.
#' @return an image with text added.
#' @keywords internal
addText <- function(image, text, color, xoff = 0, yoff = 0, corner = "TL") {
  # several checks
  if(missing(text) || (length(text) == 0)) return(image)
  color = na.omit(as.character(color))
  assert(color, len = 1, typ = "character")
  xoff = na.omit(as.integer(xoff)); xoff = xoff[xoff>=0]
  assert(xoff, len = 1, typ = "integer")
  yoff = na.omit(as.integer(yoff)); yoff = yoff[yoff>=0]
  assert(yoff, len = 1, typ = "integer")
  corner = na.omit(as.character(corner))
  assert(corner, len = 1, alw = c("TL", "TR", "BL", "BR"))
  checkColor(color)
  
  # dim + check
  di = dim(image)
  txt_msk = texttomatrix(text)
  dt = dim(txt_msk)
  if((dt[1] + yoff > di[1]) | dt[2] + xoff > di[2]) {
    new_h = dt[1]
    new_w = dt[2]
    if(dt[1] + yoff > di[1]) new_h = di[1] - yoff
    if(dt[2] + xoff > di[2]) new_w = di[2] - xoff
    txt_msk = with_seed(cpp_resize(mat = txt_msk, new_height = new_h, new_width = new_w, add_noise = FALSE, bg = 0, sd= 0), NULL)
  }
  
  # modify xoff, yoff according to corner anchorage
  switch(corner,
         "TR" = {
           xoff = di[2] - dt[2] - xoff
         }, 
         "BL" = {
           yoff = di[1] - dt[1] - yoff
         }, 
         "BR" = {
           xoff = di[2] - dt[2] - xoff
           yoff = di[1] - dt[1] - yoff
         })
  
  # place text in image
  invert = FALSE
  if(is.na(di[3])) return(cpp_mark(A = image, B = txt_msk, mask = txt_msk, xoff = xoff, yoff = yoff, invert = invert))
  color = tolower(color)
  if(color=="black") invert = TRUE
  txt_img = objectColorize(txt_msk,color)
  return(array(sapply(1:di[3], FUN=function(x) cpp_mark(A = image[,,x], B = txt_img[,,x], mask = txt_msk, xoff = xoff, yoff = yoff, invert = invert)),dim = di))
}

#' @title Image Text Incrustation
#' @name addText2
#' @description Adds Text to image.
#' @param image a [0,1] image.
#' @param text a character string.
#' @param color a character string. color of the text.
#' @param xoff x offset in image to start writing text according to its anchorage and justification. Default is \code{0}.
#' @param yoff y offset in image to start writing text according to its anchorage and justification. Default is \code{0}.
#' @param anchor numeric vector of text anchorage in normalized coordinates starting from bottom-left. Default is \code{c(1,0)}, for top-left anchorage.
#' @param vjust vertical justification of the text. Default is \code{"T"}, allowed is one of \code{"T"}, \code{"M"}, and \code{"B"} for top, middle or bottom justification, respectively.
#' @param hjust horizontal justification of the text. Default is \code{"L"}, allowed is one of \code{"L"}, \code{"C"}, and \code{"R"} for left, center or right justification, respectively.
#' @details One-lined text will be added so has to be fully contained within image.
#' @return an image with text added.
#' @keywords internal
addText2 <- function (image, text, color,
                      xoff = 0, yoff = 0, anchor = c(1,0),
                      vjust = c("T","M","B")[1], hjust = c("L","C","R")[1]) {
  if (missing(text) || (length(text) == 0)) return(image)
  color = na.omit(as.character(color))
  assert(color, len = 1, typ = "character")
  checkColor(color)
  
  anchor = na.omit(as.numeric(anchor))
  anchor = anchor[anchor >= 0 & anchor <= 1]
  assert(anchor, len = 2, typ = "numeric")
  
  xoff = na.omit(as.integer(xoff))
  assert(xoff, len = 1, typ = "integer")
  yoff = na.omit(as.integer(yoff))
  assert(yoff, len = 1, typ = "integer")
  
  hjust = as.character(hjust)
  assert(hjust, len = 1, alw = c("L","C","R"))
  vjust = as.character(vjust)
  assert(vjust, len = 1, alw = c("T","M","B"))
  
  txt_msk = texttomatrix(text)
  dt = dim(txt_msk)
  di = dim(image)
  
  YY = 1 + ceiling((1 - anchor[1]) * (di[1] - 1)) + yoff
  XX = 1 + ceiling(anchor[2] * (di[2] - 1)) + xoff
  #if(vjust == "T") YY = YY
  if(vjust == "M") YY = YY - ceiling(dt[1] / 2)
  if(vjust == "B") YY = YY - dt[1]
  #if(hjust == "L") XX = XX
  if(hjust == "C") XX = XX - ceiling(dt[2] / 2)
  if(hjust == "R") XX = XX - dt[2]
  
  invert = FALSE
  if (is.na(di[3])) return(cpp_mark2(A = image, B = txt_msk, mask = txt_msk, xoff = XX, yoff = YY, invert = invert))
  color = tolower(color)
  txt_img = objectColorize(txt_msk, color)
  return(
    array(
      sapply(1:di[3], FUN = function(x) {
        A = matrix(image[, , x, drop = TRUE], ncol = di[2])
        B = matrix(txt_img[, , x, drop = TRUE], ncol = dt[2])
        cpp_mark2(A = A, B = B, mask = txt_msk, xoff = XX, yoff = YY, invert = invert)
      }), dim = di))
}
