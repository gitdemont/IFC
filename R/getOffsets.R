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

#' @title RIF/CIF File Image Field Directories Offsets Extraction
#' @description
#' Extracts offsets of the IFDs (Image Field Directories) within a XIF file.
#' Users are highly encouraged to read TIFF specifications to have a better understanding about what offsets and IFDs are.
#' @param fileName path to file.
#' @param fast whether to fast extract 'objects' or not. Default is TRUE.\cr
#' Meaning that offsets will be extracted expecting that 'objects' are stored in ascending order.\cr
#' A message will be thrown since fast extraction method does not ensure correct mapping between objects and offsets.\cr
#' If set to FALSE, all object_ids will be scanned from 'fileName' to ensure extraction of desired offsets.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @source TIFF 6.0 specifications archived from web \url{https://web.archive.org/web/20211209104854/https://www.adobe.io/open/standards/TIFF.html}
#' @details Offsets are byte positions of IFDs found within RIF or CIF file. For more details see TIFF specifications.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   system.time(offsets_fast <- getOffsets(fileName = file_cif, fast = TRUE))
#'   system.time(offsets_slow <- getOffsets(fileName = file_cif, fast = FALSE))
#'   identical(offsets_fast, offsets_slow)   
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return an integer vector of class `IFC_offset` of IFDs offsets found in XIF file.
#' @export
getOffsets <- function(fileName, fast = TRUE, display_progress = TRUE, verbose = FALSE) {
  if(missing(fileName)) stop("'fileName' can't be missing")
  tmp = duplicated(fileName)
  if(any(tmp)) {
    warning(paste0("duplicated files have been removed from 'fileName': ","\n-", paste0(fileName[tmp],collapse="\n-")))
    fileName = fileName[!tmp]
  }
  if(length(fileName) != 1) stop("'fileName' should be of length 1")
  fast = as.logical(fast); assert(fast, len = 1, alw = c(TRUE, FALSE))
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  fileName = enc2native(normalizePath(fileName, winslash = "/", mustWork = FALSE))
  
  XIF_test = testXIF(fileName)
  obj_count = attr(XIF_test, "obj_count")
  obj_estimated = attr(XIF_test, "obj_estimated")
  if((obj_count != obj_estimated) && display_progress) message("offsets extraction is based on an estimation of total number of objects within the file.\nConsequently, it is normal behaviour if progress bar can not reach or exceed 100%.")
  # obj_count = as.integer(getInfo(fileName, warn = FALSE, force_default = TRUE, display_progress = FALSE)$objcount)
  
  type = as.integer(XIF_test == 1) + 1L
  pb = NULL
  if(fast || (XIF_test < 0)) {
    obj_pb = obj_estimated
    if(display_progress) {
      pb = newPB(min = 0, max = ifelse(obj_pb <= 0, 1, obj_pb), initial = 0, style = 3, title = basename(fileName))
      on.exit(endPB(pb))
    }
    offsets = cpp_getoffsets_noid(fileName, obj_count = obj_pb, display_progress = display_progress, pb = pb, verbose = verbose)
    offsets = offsets[offsets != 0]
    offsets_1 = offsets[1]
    offsets = offsets[-1]
    if(obj_count > 0) if(length(offsets) != obj_estimated * type) stop("Number of offsets found is different from expected object count.")
    
    message("Offsets were extracted from XIF file with fast method.\nCorrect mapping between offsets and objects ids is not guaranteed.")
  } else {
    obj_pb = obj_estimated / (as.integer(XIF_test != 1) + 1L)
    if(display_progress) {
      pb = newPB(min = 0, max = ifelse(obj_pb <= 0, 1, obj_pb), initial = 0, style = 3, title = basename(fileName))
      on.exit(endPB(pb))
    }
    offsets = as.data.frame(do.call(what = "cbind", args = cpp_getoffsets_wid(fileName, obj_count = obj_pb, display_progress = display_progress, pb = pb, verbose = verbose)), stringsAsFactors = FALSE)
    offsets = offsets[offsets$OFFSET != 0, ]
    offsets_i = offsets[offsets$TYPE == 2, ]
    offsets_m = offsets[offsets$TYPE == 3, ]
    if(type == 2) if(length(offsets_i$OBJECT_ID) != length(offsets_m$OBJECT_ID)) stop("Offsets contain different numbers of 'img' and 'msk'.")
    
    ORD_i = order(offsets_i$OBJECT_ID)
    offsets_img = offsets_i$OFFSET[ORD_i]
    ORD_m = order(offsets_m$OBJECT_ID)
    offsets_msk = offsets_m$OFFSET[ORD_m]
    
    offsets_1 = offsets[offsets$TYPE == 1, ]
    offsets_1 = offsets_1$OFFSET
    
    if(obj_count > 0) if(length(offsets_img) != obj_count) stop("Number of offsets found is different from expected object count.")
    offsets = as.numeric(apply(cbind(offsets_img, offsets_msk), 1, FUN=function(i) i))
  }
  if(length(offsets) == 0) {
    obj_count = 0
  } else {
    if(obj_count <= 0) obj_count = length(offsets) / type
    N = nchar(sprintf("%1.f",abs(obj_count-1)))
    if(type == 1) {
      names(offsets) = paste0("img_", rep(sprintf(paste0("%0",N,".f"), 0:(obj_count-1))))
    } else {
      names(offsets) = paste0(c("img_", "msk_"), rep(sprintf(paste0("%0",N,".f"), 0:(obj_count-1)), each = 2))
    }
  }
  attr(offsets, "all") = offsets
  attr(offsets, "fileName_image") = fileName
  attr(offsets, "checksum") = checksumXIF(fileName)
  attr(offsets, "obj_count") = obj_count
  attr(offsets, "first") = offsets_1 
  attr(offsets, "test") = XIF_test
  attr(offsets, "class") = c("IFC_offset")
  return(offsets)
}
