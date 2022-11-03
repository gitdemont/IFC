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

#' @title IFC_features Raw Conversion
#' @description 
#' Helper to convert features (`IFC_features` object) to raw vector.
#' @param features an `IFC_features` object.
#' @param w_con a connection opened for writing. Default is raw().
#' @param endianness The endian-ness ("big" or "little") of the target system for the file. Default is .Platform$endian.\cr
#' Endianness describes the bytes order of data stored within the files. This parameter may not be modified.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @param ... other arguments to be passed.
#' @return a raw vector of features binaries.
#' @keywords internal
toBIN_features = function(features, w_con = raw(), endianness = .Platform$endian, 
                          verbose = FALSE, display_progress = TRUE, title_progress = "", ...) {
  dots = list(...)
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("exporting features information as binary")
  assert(features, cla = "IFC_features")
  assert(endianness, alw = c("little", "big"))
  assert(display_progress, alw = c(TRUE, FALSE))
  assert(title_progress, len = 1, typ = "character")
  
  L = ncol(features)
  feat_version = cpp_uint32_to_raw(1)
  feat_number = cpp_uint32_to_raw(L)
  obj_number = cpp_uint32_to_raw(nrow(features))
  if(endianness != .Platform$endian) {
    feat_version = rev(feat_version)
    feat_number = rev(feat_number)
    obj_number = rev(obj_number)
  }
  feat_init = writeBin(con = w_con, endian = endianness, c(as.raw(c(0x0a, 0, 30)), feat_version, feat_number, obj_number))
  if(display_progress) {
    pb = newPB(min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
    if(endianness == .Platform$endian) {
      feat = lapply(seq_along(integer(L)), FUN=function(i_feat) {
        setPB(pb, value = i_feat, title = title_progress, label = "converting features values (binary)")
        writeBin(con = w_con, endian = endianness, c(cpp_uint32_to_raw(i_feat-1), writeBin(object=features[[i_feat]], con=raw(), size = 8, endian = endianness, useBytes = TRUE)))
      })
    } else {
      feat = lapply(seq_along(integer(L)), FUN=function(i_feat) {
        setPB(pb, value = i_feat, title = title_progress, label = "converting features values (binary)")
        writeBin(con = w_con, endian = endianness, c(rev(cpp_uint32_to_raw(i_feat-1)), writeBin(object=features[[i_feat]], con=raw(), size = 8, endian = endianness, useBytes = TRUE)))
      })
    }
  } else {
    if(endianness == .Platform$endian) {
      feat = lapply(seq_along(integer(L)), FUN=function(i_feat) {
        writeBin(con = w_con, endian = endianness, c(cpp_uint32_to_raw(i_feat-1), writeBin(object=features[[i_feat]], con=raw(), size = 8, endian = endianness, useBytes = TRUE)))
      })
    } else {
      feat = lapply(seq_along(integer(L)), FUN=function(i_feat) {
        writeBin(con = w_con, endian = endianness, c(rev(cpp_uint32_to_raw(i_feat-1)), writeBin(object=features[[i_feat]], con=raw(), size = 8, endian = endianness, useBytes = TRUE)))
      })
    }
  }
  return(c(feat_init, unlist(feat, recursive = FALSE, use.names = FALSE)))
}
