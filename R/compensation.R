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

################################################################################
#             functions described hereunder are experimental                   #
#              inputs and outputs may change in the future                     #
################################################################################

#' @title Shorcut for Batch RIF Images Spatial Correction and Compensation
#' @description
#' Function to shortcut extraction and RIF file spatial correction and compensation 
#' @param ... arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'bypass'(=TRUE).\cr
#' If 'param' is provided 'export'(="matrix"), 'mode'(="raw"), 'removal'(="none"),
#' 'composite'(=""), 'selection'(="all"), 'size'(=c(0, 0)), 'force_width'(=FALSE), 'spatial_correction'(=TRUE), will be overwritten.\cr
#' If 'offsets' are not provided extra arguments can also be passed with ... \code{\link{getOffsets}}.\cr
#' /!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.
#' @param objects integer vector, IDEAS objects ids numbers to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`. 
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.
#' @param spillover spillover matrix to apply. Default is missing to apply no compensation.
#' In such case, you will get images with spatial correction only. However, when set to NULL images 
#' will be compensated using spillover matrix set during acquisition. If not NULL the provided spillover 
#' matrix will be used to compensate images.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @details this function is experimental inputs and outputs may change in the future.
#' @return A list of matrices/arrays of images corresponding to objects extracted.
#' @keywords internal
CompensateFromRIF <- function(...,
                              objects,
                              offsets,
                              spillover,
                              display_progress = TRUE) { 
  dots=list(...)
  
  # check input
  input = whoami(entries = as.list(match.call()))
  if(!any(sapply(input, FUN = function(i) length(i) != 0))) {
    stop("can't determine what to extract with provided parameters.\n try to input at least one of: 'fileName', 'info', 'param' or 'offsets'")
  }
  
  # reattribute needed param
  offsets = input[["offsets"]]
  param = input[["param"]]
  if(length(offsets) == 0) {
    fileName = input[["fileName"]]
  } else {
    fileName = attr(offsets, "fileName_image")
  }
  
  # check mandatory param
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  
  # process extra parameters
  if(length(dots[["verbose"]]) == 0) { 
    verbose = FALSE
  } else {
    verbose = dots[["verbose"]]
  }
  if(length(dots[["verbosity"]]) == 0) { 
    verbosity = 1
  } else {
    verbosity = dots[["verbosity"]]
  }
  if(length(dots[["fast"]]) == 0) { 
    fast = TRUE
  } else {
    fast = dots[["fast"]]
  }
  
  fast = as.logical(fast); assert(fast, len = 1, alw = c(TRUE, FALSE))
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  verbosity = as.integer(verbosity); assert(verbosity, len = 1, alw = c(1, 2))
  param_extra = names(dots) %in% c("ifd","param","export","bypass","verbose",
                                   "mode", "removal","composite",
                                   "selection","size","force_width","spatial_correction")
  dots = dots[!param_extra] # remove not allowed param
  param_param = names(dots) %in% c("write_to","base64_id","base64_att","overwrite",
                                   "random_seed","add_noise","full_range","force_range")
  dots_param = dots[param_param] # keep param_param for objectParam
  dots = dots[!param_param]
  
  # compute object param
  # 1: prefer using 'param' if found,
  # 2: otherwise use 'info' if found,
  # 3: finally look at fileName
  filename = input$fileName
  if(length(param) == 0) {
    if(length(input$info) == 0) {
      fileName = param$fileName_image 
    } else {
      fileName = input$info$fileName_image
    }
  } else {
    fileName = param$fileName_image
  }
  if(length(fileName) == 0) fileName = attr(input$offsets, "fileName_image")
  
  # compute object param
  info = getInfo(fileName = fileName, from = "acquisition", warn = FALSE)
  param = objectParam(info = info, selection = "all", 
                      add_noise = FALSE, removal = "none",
                      mode = "raw", export = "matrix", 
                      size = c(0,0), force_width = FALSE, spatial_correction = TRUE,
                      force_range = FALSE)
  
  if(!missing(spillover)) {
    if(length(spillover) == 0) {
      spillover = info$CrossTalkMatrix[which(info$in_use), which(info$in_use)]
    } else {
      if((ncol(spillover) != sum(info$in_use)) || (nrow(spillover) != sum(info$in_use))) stop("'spillover' differs from number of acquired channels")
    }
    spill = spillover
  } else {
    spill = NULL
  }
  
  fileName = param$fileName_image
  title_progress = basename(fileName)
  
  # check input offsets if any
  compute_offsets = TRUE
  if(missing(offsets)) offsets = NULL
  if(length(offsets) != 0) {
    if(!("IFC_offset" %in% class(offsets))) {
      warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
    } else {
      if(attr(offsets, "checksum") != checksumXIF(param$fileName_image)) {
        warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
      } else {
        compute_offsets = FALSE
      }
    }
  }
  if(compute_offsets) {
    offsets = suppressMessages(getOffsets(fileName = param$fileName_image, fast = fast, display_progress = display_progress, verbose = verbose))
  }
  
  # check objects to extract
  nobj = as.integer(attr(x = offsets, which = "obj_count"))
  
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }
  
  # extract objects
  sel = subsetOffsets(offsets = offsets, objects = objects, image_type = "img")
  sel = split(sel, ceiling(seq_along(sel)/20))
  L=length(sel)
  if(L == 0) {
    warning("CompensateFromRIF: No objects to extract, check the objects you provided.", immediate. = TRUE, call. = FALSE)
    return(NULL)
  }
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
    ans = lapply(1:L, FUN=function(i) {
      setPB(pb, value = i, title = title_progress, label = "exporting images to matrix")
      pre = do.call(what = "objectExtract", args = c(list(ifd = lapply(sel[[i]],
                                                                 FUN = function(off) cpp_getTAGS(fname = param$fileName_image,
                                                                                                 offset = off,
                                                                                                 trunc_bytes = 1, 
                                                                                                 force_trunc = TRUE, 
                                                                                                 verbose = verbose)),
                                                    param = param,
                                                    verbose = verbose,
                                                    bypass = TRUE),
                                               dots))
      lapply(1:length(pre), FUN = function(i_pre) {
        d = dim(pre[[i_pre]][[1]]) - c(3, 2)
        if(!is.null(spill)) foo = compensate(c(pre[[i_pre]][[i_ch]]), spill)
        foo = lapply(1:length(pre[[i_pre]]), FUN = function(i_ch) {
          img = matrix(foo[, i_ch], nrow = d[1], ncol = d[2])
          attr(img, which = "mask") <- attr(pre[[i_pre]][[i_ch]], which = "mask")
          attr(img, "input_range") <- attr(pre[[i_pre]][[i_ch]], "input_range")
          attr(img, "full_range") <- attr(pre[[i_pre]][[i_ch]], "full_range")
          attr(img, "force_range") <- attr(pre[[i_pre]][[i_ch]], "force_range")
          attr(img, "gamma") <- attr(pre[[i_pre]][[i_ch]], "gamma")
          attr(img, "color") <- attr(pre[[i_pre]][[i_ch]], "color")
          attr(img, "removal") <- attr(pre[[i_pre]][[i_ch]], "removal")
          attr(img, "BG_MEAN") <- attr(pre[[i_pre]][[i_ch]], "BG_MEAN")
          attr(img, "BG_STD") <- attr(pre[[i_pre]][[i_ch]], "BG_STD")
          return(img)
        })
        attr(foo, "names") <- attr(pre[[i_pre]], "names")
        attr(foo, "object_id") <- attr(pre[[i_pre]], "object_id")
        attr(foo, "offset_id") <- attr(pre[[i_pre]], "offset_id")
        attr(foo, "channel_id") <- attr(pre[[i_pre]], "channel_id")
        attr(foo, "removal") <- attr(pre[[i_pre]], "removal")
        return(foo)
      })
    })
  } else {
    ans = lapply(1:L, FUN=function(i) { 
      pre = do.call(what = "objectExtract", args = c(list(ifd = lapply(sel[[i]],
                                                                 FUN = function(off) cpp_getTAGS(fname = param$fileName_image,
                                                                                                 offset = off,
                                                                                                 trunc_bytes = 1, 
                                                                                                 force_trunc = TRUE, 
                                                                                                 verbose = verbose)),
                                                    param = param,
                                                    verbose = verbose,
                                                    bypass = TRUE),
                                               dots))
      lapply(1:length(pre), FUN = function(i_pre) {
        d = dim(pre[[i_pre]][[1]]) - c(3, 2)
        if(!is.null(spill)) foo = compensate(c(pre[[i_pre]][[i_ch]]), spill)
        foo = lapply(1:length(pre[[i_pre]]), FUN = function(i_ch) {
          img = matrix(foo[, i_ch], nrow = d[1], ncol = d[2])
          attr(img, which = "mask") <- attr(pre[[i_pre]][[i_ch]], which = "mask")
          attr(img, "input_range") <- attr(pre[[i_pre]][[i_ch]], "input_range")
          attr(img, "full_range") <- attr(pre[[i_pre]][[i_ch]], "full_range")
          attr(img, "force_range") <- attr(pre[[i_pre]][[i_ch]], "force_range")
          attr(img, "gamma") <- attr(pre[[i_pre]][[i_ch]], "gamma")
          attr(img, "color") <- attr(pre[[i_pre]][[i_ch]], "color")
          attr(img, "removal") <- attr(pre[[i_pre]][[i_ch]], "removal")
          attr(img, "BG_MEAN") <- attr(pre[[i_pre]][[i_ch]], "BG_MEAN")
          attr(img, "BG_STD") <- attr(pre[[i_pre]][[i_ch]], "BG_STD")
          return(img)
        })
        attr(foo, "names") <- attr(pre[[i_pre]], "names")
        attr(foo, "object_id") <- attr(pre[[i_pre]], "object_id")
        attr(foo, "offset_id") <- attr(pre[[i_pre]], "offset_id")
        attr(foo, "channel_id") <- attr(pre[[i_pre]], "channel_id")
        attr(foo, "removal") <- attr(pre[[i_pre]], "removal")
        return(foo)
      })
    }) 
  }
  channel_id = attr(ans[[1]][[1]], "channel_id")
  if(L>1) {
    ans = do.call(what="c", args=ans)
  } else {
    ans = ans[[1]]
  }
  ids = sapply(ans, attr, which="object_id")
  if(!all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  names(ans) = objects
  attr(ans, "channel_id") <- channel_id
  return(ans)
}

#' @title Shorcut for Batch CIF Images Compensation
#' @description
#' Function to shortcut extraction and CIF file compensation / decompensation 
#' @param ... arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'bypass'(=TRUE).\cr
#' If 'param' is provided 'export'(="matrix"), 'mode'(="raw"), 'removal'(="none"),
#' 'composite'(=""), 'selection'(="all"), 'size'(=c(0, 0)), 'force_width'(=FALSE), 'spatial_correction'(=FALSE) will be overwritten.\cr
#' If 'offsets' are not provided extra arguments can also be passed with ... \code{\link{getOffsets}}.\cr
#' /!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.
#' @param objects integer vector, IDEAS objects ids numbers to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`. 
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.
#' @param spillover spillover matrix to apply. Default is missing to do nothing.
#' When set to NULL images will be decompensated using spillover matrix set within .cif file. Otherwise,
#' the images will be decompensated with spillover matrix set within the file and recompensated with the one provided.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @details this function is experimental inputs and outputs may change in the future.
#' @return A list of matrices/arrays of images corresponding to objects extracted.
#' @keywords internal
CompensateFromCIF <- function(...,
                              objects,
                              offsets,
                              spillover,
                              display_progress = TRUE) { 
  dots=list(...)
  
  # check input
  input = whoami(entries = as.list(match.call()))
  if(!any(sapply(input, FUN = function(i) length(i) != 0))) {
    stop("can't determine what to extract with provided parameters.\n try to input at least one of: 'fileName', 'info', 'param' or 'offsets'")
  }
  
  # reattribute needed param
  offsets = input[["offsets"]]
  param = input[["param"]]
  if(length(offsets) == 0) {
    fileName = input[["fileName"]]
  } else {
    fileName = attr(offsets, "fileName_image")
  }
  
  # check mandatory param
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  
  # process extra parameters
  if(length(dots[["verbose"]]) == 0) { 
    verbose = FALSE
  } else {
    verbose = dots[["verbose"]]
  }
  if(length(dots[["verbosity"]]) == 0) { 
    verbosity = 1
  } else {
    verbosity = dots[["verbosity"]]
  }
  if(length(dots[["fast"]]) == 0) { 
    fast = TRUE
  } else {
    fast = dots[["fast"]]
  }
  
  fast = as.logical(fast); assert(fast, len = 1, alw = c(TRUE, FALSE))
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  verbosity = as.integer(verbosity); assert(verbosity, len = 1, alw = c(1, 2))
  param_extra = names(dots) %in% c("ifd","param","export","bypass","verbose",
                                   "mode", "removal","composite",
                                   "selection","size","force_width","spatial_correction")
  dots = dots[!param_extra] # remove not allowed param
  param_param = names(dots) %in% c("write_to","base64_id","base64_att","overwrite",
                                   "random_seed","add_noise","full_range","force_range")
  dots_param = dots[param_param] # keep param_param for objectParam
  dots = dots[!param_param]
  
  # compute object param
  # 1: prefer using 'param' if found,
  # 2: otherwise use 'info' if found,
  # 3: finally look at fileName
  filename = input$fileName
  if(length(param) == 0) {
    if(length(input$info) == 0) {
      fileName = param$fileName_image 
    } else {
      fileName = input$info$fileName_image
    }
  } else {
    fileName = param$fileName_image
  }
  if(length(fileName) == 0) fileName = attr(input$offsets, "fileName_image")
  
  # compute object param
  info = getInfo(fileName = fileName, from = "analysis", warn = FALSE)
  param = objectParam(info = info, selection = "all", 
                      add_noise = FALSE, removal = "none",
                      mode = "raw", export = "matrix", 
                      size = c(0,0), force_width = FALSE, spatial_correction = FALSE,
                      force_range = FALSE)
  
  if(missing(spillover)) {
    spill = 0
  } else {
    if(length(spillover) == 0) {
      spill = NULL
    } else {
      if((ncol(spillover) != sum(info$in_use)) || (nrow(spillover) != sum(info$in_use))) stop("'spillover' differs from number of acquired channels")
      spill = spillover 
    }
  }
  decomp = info$CrossTalkMatrix[which(info$in_use), which(info$in_use)]
  
  # check input offsets if any
  compute_offsets = TRUE
  if(missing(offsets)) offsets = NULL
  if(length(offsets) != 0) {
    if(!("IFC_offset" %in% class(offsets))) {
      warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
    } else {
      if(attr(offsets, "checksum") != checksumXIF(param$fileName_image)) {
        warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
      } else {
        compute_offsets = FALSE
      }
    }
  }
  if(compute_offsets) {
    offsets = suppressMessages(getOffsets(fileName = param$fileName_image, fast = fast, display_progress = display_progress, verbose = verbose))
  }
  
  # check objects to extract
  nobj = as.integer(attr(x = offsets, which = "obj_count"))
  
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }
  
  # extract objects
  sel = subsetOffsets(offsets = offsets, objects = objects, image_type = "img")
  sel = split(sel, ceiling(seq_along(sel)/20))
  L=length(sel)
  if(L == 0) {
    warning("CompensateFromCIF: No objects to extract, check the objects you provided.", immediate. = TRUE, call. = FALSE)
    return(NULL)
  }
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
    ans = lapply(1:L, FUN=function(i) {
      setPB(pb, value = i, title = title_progress, label = "exporting images to matrix")
      pre = do.call(what = "objectExtract", args = c(list(ifd = lapply(sel[[i]],
                                                                 FUN = function(off) cpp_getTAGS(fname = param$fileName_image,
                                                                                                 offset = off,
                                                                                                 trunc_bytes = 1, 
                                                                                                 force_trunc = TRUE, 
                                                                                                 verbose = verbose)),
                                                    param = param,
                                                    verbose = verbose,
                                                    bypass = TRUE),
                                               dots))
      lapply(1:length(pre), FUN = function(i_pre) {
        d = dim(pre[[i_pre]][[1]])
        foo = sapply(1:length(pre[[i_pre]]), FUN = function(i_ch) {
          return(c(attr(pre[[i_pre]][[i_ch]], which = "raw")))
        })
        if(is.null(spill)) {
          foo = decompensate(foo, decomp)
        } else {
          if(length(spill) != 1) {
            foo = recompensate(foo, decomp, spill)
          }
        }
        foo = lapply(1:length(pre[[i_pre]]), FUN = function(i_ch) {
          img = matrix(foo[, i_ch], nrow = d[1], ncol = d[2])
          attr(img, "mask") <- attr(pre[[i_pre]][[i_ch]], "mask")
          attr(img, "input_range") <- attr(pre[[i_pre]][[i_ch]], "input_range")
          attr(img, "full_range") <- attr(pre[[i_pre]][[i_ch]], "full_range")
          attr(img, "force_range") <- attr(pre[[i_pre]][[i_ch]], "force_range")
          attr(img, "gamma") <- attr(pre[[i_pre]][[i_ch]], "gamma")
          attr(img, "color") <- attr(pre[[i_pre]][[i_ch]], "color")
          attr(img, "removal") <- attr(pre[[i_pre]][[i_ch]], "removal")
          attr(img, "BG_MEAN") <- attr(pre[[i_pre]][[i_ch]], "BG_MEAN")
          attr(img, "BG_STD") <- attr(pre[[i_pre]][[i_ch]], "BG_STD")
          return(img)
        })
        attr(foo, "names") <- attr(pre[[i_pre]], "names")
        attr(foo, "object_id") <- attr(pre[[i_pre]], "object_id")
        attr(foo, "offset_id") <- attr(pre[[i_pre]], "offset_id")
        attr(foo, "channel_id") <- attr(pre[[i_pre]], "channel_id")
        attr(foo, "removal") <- attr(pre[[i_pre]], "removal")
        return(foo)
      })
    })
  } else {
    ans = lapply(1:L, FUN=function(i) { 
      pre = do.call(what = "objectExtract", args = c(list(ifd = lapply(sel[[i]],
                                                                 FUN = function(off) cpp_getTAGS(fname = param$fileName_image,
                                                                                                 offset = off,
                                                                                                 trunc_bytes = 1, 
                                                                                                 force_trunc = TRUE, 
                                                                                                 verbose = verbose)),
                                                    param = param,
                                                    verbose = verbose,
                                                    bypass = TRUE),
                                               dots))
      lapply(1:length(pre), FUN = function(i_pre) {
        d = dim(pre[[i_pre]][[1]])
        foo = sapply(1:length(pre[[i_pre]]), FUN = function(i_ch) {
          return(c(attr(pre[[i_pre]][[i_ch]], which = "raw")))
        })
        if(is.null(spill)) {
          foo = decompensate(foo, decomp)
        } else {
          if(length(spill) != 1) {
            foo = recompensate(foo, decomp, spill)
          }
        }
        foo = lapply(1:length(pre[[i_pre]]), FUN = function(i_ch) {
          img = matrix(foo[, i_ch], nrow = d[1], ncol = d[2])
          attr(img, "mask") <- attr(pre[[i_pre]][[i_ch]], "mask")
          attr(img, "input_range") <- attr(pre[[i_pre]][[i_ch]], "input_range")
          attr(img, "full_range") <- attr(pre[[i_pre]][[i_ch]], "full_range")
          attr(img, "force_range") <- attr(pre[[i_pre]][[i_ch]], "force_range")
          attr(img, "gamma") <- attr(pre[[i_pre]][[i_ch]], "gamma")
          attr(img, "color") <- attr(pre[[i_pre]][[i_ch]], "color")
          attr(img, "removal") <- attr(pre[[i_pre]][[i_ch]], "removal")
          attr(img, "BG_MEAN") <- attr(pre[[i_pre]][[i_ch]], "BG_MEAN")
          attr(img, "BG_STD") <- attr(pre[[i_pre]][[i_ch]], "BG_STD")
          return(img)
        })
        attr(foo, "names") <- attr(pre[[i_pre]], "names")
        attr(foo, "object_id") <- attr(pre[[i_pre]], "object_id")
        attr(foo, "offset_id") <- attr(pre[[i_pre]], "offset_id")
        attr(foo, "channel_id") <- attr(pre[[i_pre]], "channel_id")
        attr(foo, "removal") <- attr(pre[[i_pre]], "removal")
        return(foo)
      })
    }) 
  }
  channel_id = attr(ans[[1]][[1]], "channel_id")
  if(L>1) {
    ans = do.call(what="c", args=ans)
  } else {
    ans = ans[[1]]
  }
  ids = sapply(ans, attr, which="object_id")
  if(!all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  names(ans) = objects
  attr(ans, "channel_id") <- channel_id
  return(ans)
}

decompensate <- function(comp, spillover) {
  # compared with flowcore spillover is t()
  # N = colnames(comp)
  # if(any(grepl("^Raw", N))) stop("'comp' is already decompensated")
  ans = t(solve(a = solve(spillover), b = t(as.matrix(comp))))
  # colnames(ans) = paste0("Raw ", N)
  return(ans)
}

compensate <- function(raw, spillover) {
  # compared with flowcore spillover is t()
  # N = colnames(raw)
  # if(!all(grepl("^Raw", N))) stop("'raw' is not raw")
  ans = as.matrix(raw) %*% solve(t(spillover))
  # colnames(ans) = gsub("^Raw ", "", N)
  return(ans)
}

recompensate <- function(val, spill_dec, spill_comp) {
  # compared with flowcore spillover is t()
  ans = compensate(as.matrix(val), solve(spill_dec, spill_comp))
  return(ans)
}