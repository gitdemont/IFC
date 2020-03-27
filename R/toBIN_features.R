#' @title IFC_features Raw Conversion
#' @description 
#' Helper to convert features (IFC_features object) to raw vector.
#' @param features an IFC_features object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @return an raw vector of features binaries.
#' @keywords internal
toBIN_features = function(features, endianness = .Platform$endian, 
                          verbose = FALSE, display_progress = TRUE, title_progress = "") {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("exporting features information as binary")
  assert(features, cla = "IFC_features")
  assert(endianness, alw = c("little", "big"))
  assert(display_progress, alw = c(TRUE, FALSE))
  assert(title_progress, len = 1, typ = "character")
  
  L = ncol(features)
  feat_version = packBits(intToBits(1),"raw")
  feat_number = packBits(intToBits(L),"raw")
  obj_number = packBits(intToBits(nrow(features)),"raw")
  if(endianness != .Platform$endian) {
    feat_version = rev(feat_version)
    feat_number = rev(feat_number)
    obj_number = rev(obj_number)
  }
  if(display_progress) {
    pb = newPB(min = 0, max = 1, initial = 0, style = 3)
    on.exit(endPB(pb))
    if(endianness == .Platform$endian) {
      feat = lapply(1:L, FUN=function(i_feat) {
        setPB(pb, value = i_feat/L, title = title_progress, label = "converting features values (binary)")
        c(packBits(intToBits(i_feat-1),"raw"), writeBin(object=features[[i_feat]], con=raw(), size = 8, endian = endianness, useBytes = TRUE))
      })
    } else {
      feat = lapply(1:L, FUN=function(i_feat) {
        setPB(pb, value = i_feat/L, title = title_progress, label = "converting features values (binary)")
        c(rev(packBits(intToBits(i_feat-1),"raw")), writeBin(object=features[[i_feat]], con=raw(), size = 8, endian = endianness, useBytes = TRUE))
      })
    }
  } else {
    if(endianness == .Platform$endian) {
      feat = lapply(1:L, FUN=function(i_feat) {
        c(packBits(intToBits(i_feat-1),"raw"), writeBin(object=features[[i_feat]], con=raw(), size = 8, endian = endianness, useBytes = TRUE))
      })
    } else {
      feat = lapply(1:L, FUN=function(i_feat) {
        c(rev(packBits(intToBits(i_feat-1),"raw")), writeBin(object=features[[i_feat]], con=raw(), size = 8, endian = endianness, useBytes = TRUE))
      })
    }
  }
  return(c(feat_version, feat_number, obj_number, unlist(feat)))
}
