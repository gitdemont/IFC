################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2022 Yohann Demont                                             #
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

#' @title IFC Image Coercion
#' @description
#' Helper to build a list image values to allow export
#' @param physicalChannel the channel. Default is 1. Allowed are [1-12].
#' @param ... other arguments to be passed. See members in returned value.
#' @param BF should 'physicalChannel' channel be considered as brightfield. Default is FALSE.
#' @param MODE collection mode (as retrieved by getInfo) determining the range. Default is 1.
#' @return a list whose members are "name", "color", "physicalChannel", "xmin", "xmax",
#' "xmid", "ymid", "scalemin", "scalemax", "tokens", "baseimage", "function", "saturation".
#' @keywords internal
buildImage <- function(physicalChannel = 1, ..., BF = FALSE, MODE = 1) {
  dots = list(...)
  
  # some checks
  physicalChannel = as.integer(physicalChannel); assert(physicalChannel, len = 1, alw = 1L:12L)
  assert(BF, len = 1, alw = c(TRUE,FALSE))
  assert(MODE, len = 1)
  
  # define som variables
  int_names =  c("xmin", "xmax", "xmid", "ymid", "scalemin", "scalemax")
  def_cols = rep(c("DarkOrchid", "Lime", "Yellow", "DarkOrange", "Red", "DeepPink"), 2)
  
  # create default returned object
  m = ifelse(MODE == 0, 1023L, 4095L)
  if(BF) {
    if(MODE == 0) {
      ans = c(100,  300,  200,  127,   95,  305)
    } else {
      ans = c(450, 1000,  725,  127,  445, 1005)
    }
  } else {
    ans = c(0, m, m / 2, 127, 0, m)
  }
  ans = c(list(sprintf("Ch%02i",physicalChannel),def_cols[physicalChannel]),
          as.list(c(physicalChannel,as.integer(ans))),
          list("","","","Cyan"))
  names(ans) = c("name", "color", "physicalChannel",
                 int_names,
                 "tokens", "baseimage", "function", "saturation")
  if(BF) ans$color = "White"
  
  # add / validate values passed from ...
  for(i in names(ans)) {
    if(i %in% names(dots)) {
      if(i %in% int_names) {
        v = na.omit(as.integer(dots[[i]]))
        if(length(v) == 1) {
          if(v < 0) stop("'",i,"' [", v,"] should be >= 0.")
          if(i == "ymid") {
            if(v > 255) stop("'",i,"' [", v,"] should be [0-255].")
          } else {
            if(v > m) stop("'",i,"' [", v,"] should be [0-",m,"], with MODE [",MODE,"].")
          }
          ans[[i]] <- v 
        } else {
          stop("'",i,"' [", dots[[i]],"] should be coercible to a unique integer.")
        }
      } else {
        v = na.omit(as.character(dots[[i]]))
        if(length(v) == 1) {
          if(i %in% c("color","saturation")) {
            alw_cols = unique(unlist(IFC::paletteIFC()[,c("color_R","lightModeColor_R")]))
            if(!(v %in% alw_cols)) stop("'",i,"' [", v,"] color is not allowed. Allowed are: ", paste0(alw_cols, collapse = ", "),".")
          } else {
            if((i != "name") && (v != "")) stop("'",i,"' [", v,"] is not allowed. Allowed is \"\".")
          }
          ans[[i]] <- v
        } else {
          stop("'",i,"' [", dots[[i]],"] should be coercible to a unique non-empty string.")
        }
      }
    }
  }
  ans
}