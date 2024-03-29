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

#' @title Object Transformation
#' @description
#' Function to normalize, colorize and add background to images.
#' @param mat a finite numeric matrix.
#' @param msk a finite numeric matrix (mask identifying abnormalities). If missing, the default no cleansing will be done.
#' @param color a color.
#' @param input_range a finite numeric vector of 2 values, sets the range of the input intensity values. Values outside this range are clipped. Default is \code{[0,4095]}.
#' @param mode color mode export. Either \code{"rgb"}, \code{"gray"} or \code{"raw"}.
#' @param type image object type.
#' @param add_noise whether to add normal noise to background or not. Default is \code{TRUE}.
#' @param random_seed a list of elements to pass to \link[base]{set.seed} or a single value, interpreted as an integer, or NULL to be used when \code{'add_noise'} is set to \code{TRUE}. Default is \code{NULL}.
#' Note that \code{NA_integer_} or \code{list(seed = NA_integer_)} can be used to not call \link[base]{set.seed} at all.
#' @param size a length 2 integer vector of final dimensions of the image, height 1st and width 2nd. Default is \code{c(0,0)} for no change.
#' @param bg_mean mean value of the background added if \code{'add_noise'} is \code{TRUE}. Default is \code{0}.
#' @param bg_sd standard deviation of the background added if \code{'add_noise'} is \code{TRUE}. Default is \code{0}.
#' @param full_range only apply when \code{'mode'} is not \code{"raw"}, if \code{'full_range'} is \code{TRUE}, then \code{'input_range'} will be set to \code{[0,4095]} and \code{'gamma'} forced to \code{1}. Default is \code{FALSE}.
#' @param force_range bool, only apply when \code{'mode'} is not \code{"raw"}, if \code{'force_range'} is \code{TRUE}, then \code{'input_range'} will be adjusted to object range in \code{[-4095,+inf]} and \code{'gamma'} forced to \code{1}. Default is \code{FALSE}.\cr
#' Note that this parameter takes the precedence over \code{'input_range'} and \code{'full_range'}.
#' @param gamma correction. Default is \code{1}, for no correction.
#' @details When \code{'add_noise'} is \code{FALSE} and \code{attr(msk, "removal")} is \code{"masked"} or \code{"MC"},
#' background will be automatically set to minimal pixel value.
#' @return the matrix transformed according to input parameters
#' @keywords internal
objectTransform <- function(mat, msk, color, input_range, mode, type, 
                            add_noise = TRUE, random_seed = NULL, size = c(0,0),
                            bg_mean = 0, bg_sd = 0, full_range = FALSE, force_range = FALSE, gamma = 1) {
  foo = mat
  bg_2 = bg_mean
  sd_2 = bg_sd
  if(!missing(msk)) {
    if((attr(msk, "removal") %in% c("masked","raw")) && (!add_noise)) {
      bg_2 = -4096
      sd_2 = 0
    }
    foo = objectCleanse(mat = foo, msk = msk, add_noise = add_noise, random_seed = random_seed, bg = bg_mean, sd = bg_sd)
  }
  foo = objectResize(mat = foo, size = size, add_noise = add_noise, random_seed = random_seed, bg = bg_mean, sd = bg_sd)
  switch(mode,
         "gray" = {
           foo = objectNormalize(mat = foo, input_range = input_range, full_range = full_range, force_range = force_range, gamma = gamma)
         },
         "rgb" = {
           foo = objectColorize(mat = objectNormalize(mat = foo, input_range = input_range, full_range = full_range, force_range = force_range, gamma = gamma), color)
         })
  attr(foo, "input_range") <- input_range
  attr(foo, "full_range") <- full_range
  attr(foo, "force_range") <- force_range
  attr(foo, "gamma") <- gamma
  attr(foo, "color") <- color
  attr(foo, "mode") <- mode
  attr(foo, "raw") <- mat
  attr(foo, "BG_MEAN") <- bg_mean
  attr(foo, "BG_STD") <- bg_sd
  if(type == 2) {
    attr(foo, "class") <- "IFC_img"
  } else {
    attr(foo, "class") <- "IFC_msk"
  }
  return(foo)
}
