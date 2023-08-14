################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2023 Yohann Demont                                             #
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

#' @title Shortcut for Batch Images or Masks Extraction
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of images or masks to various format.
#' @param \dots arguments to be passed to \code{\link{objectExtract}} with the exception of \code{'ifd'} and \code{'bypass'}(=\strong{TRUE}).\cr
#' \strong{/!\\} If not any of \code{'fileName'}, \code{'info'} and \code{'param'} can be found in \code{'...'} then \code{attr(offsets, "fileName_image")} will be used as \code{'fileName'} input parameter to pass to \code{\link{objectParam}}.
#' @param objects integer vector, IDEAS objects ids numbers to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`.
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.\cr
#' If \code{'offsets'} are not provided, extra arguments can also be passed with \code{'...'} to \code{\link{getOffsets}}.
#' @param display_progress whether to display a progress bar. Default is \code{TRUE}.
#' @param image_type (\code{\link{subsetOffsets}} argument) type of desired object offsets. Default is \code{"img"}. Allowed are \code{"img"} or \code{"msk"}.
#' @param export (\code{\link{objectParam}} argument) format mode export. Either \code{"file"}, \code{"matrix"}, \code{"base64"}. Default is \code{"matrix"}.
#' @param mode (\code{\link{objectParam}} argument) color mode export. Either \code{"raw"}, \code{"rgb"}, \code{"gray"}. Default is \code{"raw"}.
#' @note Arguments of \code{\link{objectExtract}} will be deduced from input arguments.
#' @details If \code{'param'} is provided in \code{'...'}, \code{'param$export'<-'export'} and \code{'param$mode'<-'mode'} \strong{only} will be overwritten.
#' @inherit objectExtract return
#' @keywords internal
polyExtractTo <- function(...,
                          objects,
                          offsets,
                          display_progress = TRUE,
                          image_type = "img",
                          export = "matrix",
                          mode = "raw") {
  # check mandatory parameters
  assert(image_type, len = 1, alw = c("img", "msk"))
  assert(export, len = 1, alw = c("file", "matrix", "base64"))
  
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
    fileName = enc2native(input[["fileName"]])
  } else {
    fileName = enc2native(attr(offsets, "fileName_image"))
  }
  
  # process extra parameters
  if(length(dots[["display_progress"]]) == 0) { 
    display_progress = TRUE
  } else {
    display_progress = dots[["display_progress"]]
  }
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
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  
  # force 'mode' and 'write_to' parameters depending on 'export'
  if(export == "matrix") {
    if(length(param) != 0) mode = param$mode
    assert(mode, len = 1, alw = c("rgb", "gray", "raw"))
    write_to = "%o_%c.bmp"
  } else {
    assert(mode, len = 1, alw = c("rgb", "gray"))
    write_to = NULL
    if(export == "base64") write_to = "%o_%c.bmp"
    if(length(param) != 0) write_to = param$write_to
    if(length(dots[["write_to"]]) != 0) write_to = dots[["write_to"]]
    if(length(write_to) == 0) stop("'write_to' can't be missing")
    write_to = na.omit(as.character(write_to))
    assert(write_to, len = 1, typ = "character")
  }
  
  # define arguments input to param
  param_extra_a = setdiff(names(formals(objectParam)), "...")
  param_extra_n = c("ifd","param","bypass","verbose", # arguments to objectExtract
                    "export","mode","write_to")       # arguments to objectParam
  param_extra = names(dots) %in% param_extra_n
  dots = dots[!param_extra] # remove not allowed param
  param_extra_a = setdiff(param_extra_a, param_extra_n)
  param_param = names(dots) %in% param_extra_a
  dots_param = dots[param_param] # keep param_param for objectParam
  dots = dots[!param_param]
  
  # compute object param
  # 1: prefer using 'param' if found,
  # 2: otherwise use 'info' if found,
  # 3: finally look at fileName
  if(length(param) == 0) {
    if(length(input$info) == 0) { 
      param = do.call(what = "objectParam",
                      args = c(list(fileName = fileName,
                                    export = export,
                                    write_to = write_to,
                                    mode = mode), dots_param))
    } else {
      param = do.call(what = "objectParam",
                      args = c(list(info = quote(input$info),
                                    export = export,
                                    write_to = write_to,
                                    mode = mode), dots_param))
    }
  } else {
    param = input$param
    param$export = export
    param$mode = mode
    if(export != "matrix") {
      param$write_to <- write_to
      param$splitp_obj = splitp(param$write_to)
      param$dir_name <- dirname(formatn(splitp_obj = splitp(param$write_to), splitf_obj = param$splitf_obj))
      type = getFileExt(param$write_to)
      switch(type,
             "jpg" = {type <- "jpeg"},
             "tif" = {type <- "tiff"})
      ##### check type
      assert(type, len = 1, alw = c("bmp", "jpeg", "png", "tiff"))
      param$type <- type
    }
  }
  param$fileName_image = enc2native(param$fileName_image)
  fileName = param$fileName_image
  title_progress = basename(fileName)
  
  # check input offsets if any
  compute_offsets = TRUE
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
      warning("Some objects that are not in ", fileName, "\nhave been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "), immediate. = TRUE, call. = FALSE)
      objects = objects[tokeep]
    }
  }
  
  # extract objects
  sel = subsetOffsets(offsets = offsets, objects = objects, image_type = image_type)
  sel = split(sel, ceiling(seq_along(sel)/20))
  L=length(sel)
  if(L == 0) {
    warning("No objects to extract, check the objects you provided.", immediate. = TRUE, call. = FALSE)
    return(invisible(NULL))
  }
  if(display_progress) {
    pb = newPB(min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
    ans = lapply(1:L, FUN=function(i) {
      setPB(pb, value = i, title = title_progress, label = "exporting images to file")
      do.call(what = "objectExtract", args = c(list(ifd = lapply(sel[[i]],
                                                                 FUN = function(off) cpp_getTAGS(fname = param$fileName_image,
                                                                                                 offset = off,
                                                                                                 trunc_bytes = 1, 
                                                                                                 force_trunc = TRUE, 
                                                                                                 verbose = verbose)),
                                                    param = param,
                                                    verbose = verbose,
                                                    bypass = TRUE),
                                               dots))
    })
  } else{
    ans = lapply(1:L, FUN=function(i) {
      do.call(what = "objectExtract", args = c(list(ifd = lapply(sel[[i]],
                                                                 FUN = function(off) cpp_getTAGS(fname = param$fileName_image,
                                                                                                 offset = off,
                                                                                                 trunc_bytes = 1, 
                                                                                                 force_trunc = TRUE, 
                                                                                                 verbose = verbose)),
                                                    param = param,
                                                    verbose = verbose,
                                                    bypass = TRUE),
                                               dots))
    })
  }
  channel_id = attr(ans[[1]][[1]], "channel_id")
  if(L>1) {
    ans = do.call(what="c", args=ans)
  } else {
    ans = ans[[1]]
  }
  if(image_type == "img") {
    ids = sapply(ans, attr, which="object_id")
  } else {
    ids = as.integer(gsub("msk_", "", gsub("img_", "", sapply(ans, attr, which = "offset_id"), fixed = TRUE), fixed = TRUE))
  }
  if((length(ids) != 0) && !all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  names(ans) = num_to_string(objects)
  if(export != "file") {
    attr(ans, "channel_id") <- channel_id
    return(ans)
  }
  return(invisible(ans))
}

#' @title Shortcut for Batch Images Extraction to Matrices/Arrays
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of images to matrix ! excludes mask.
#' @param \dots arguments to be passed to \code{\link{objectExtract}} with the exception of \code{'ifd'} and \code{'bypass'}(=\strong{TRUE}).\cr
#' \strong{/!\\} If not any of \code{'fileName'}, \code{'info'} and \code{'param'} can be found in \code{'...'} then \code{attr(offsets, "fileName_image")} will be used as \code{'fileName'} input parameter to pass to \code{\link{objectParam}}.
#' @param objects integer vector, IDEAS objects ids numbers to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`.
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.\cr
#' If \code{'offsets'} are not provided, extra arguments can also be passed with \code{'...'} to \code{\link{getOffsets}}.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @note Arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractImages_toMatrix}} input arguments.
#' @details If \code{'param'} is provided in \code{'...'}, \code{'param$export'<-"matrix"} \strong{only} will be overwritten.
#' @return A list of matrices/arrays of images corresponding to objects extracted.
#' @export
ExtractImages_toMatrix <- function(...,
                                   objects,
                                   offsets,
                                   display_progress = TRUE) { 
  dots = list(...)
  dots = dots[!(names(dots) %in% c("image_type", "export"))]
  args = list(image_type = "img", export = "matrix")
  if(!missing(objects)) args = c(args, list(objects = objects))
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  args = c(args, list(display_progress = display_progress))
  do.call(what = polyExtractTo, args = c(dots, args))
}

#' @title Shortcut for Batch Masks Extraction to Matrices/Arrays
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of masks to matrix ! excludes image.
#' @inheritParams ExtractImages_toMatrix
#' @note Arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractMasks_toMatrix}} input arguments.
#' @inherit ExtractImages_toMatrix details return
#' @export
ExtractMasks_toMatrix <- function(...,
                                  objects,
                                  offsets,
                                  display_progress = TRUE) {
  dots = list(...)
  dots = dots[!(names(dots) %in% c("image_type", "export"))]
  args = list(image_type = "msk", export = "matrix")
  if(!missing(objects)) args = c(args, list(objects = objects))
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  args = c(args, list(display_progress = display_progress))
  do.call(what = polyExtractTo, args = c(dots, args))
}

#' @title Shortcut for Batch Images Extraction to Base64
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of images to base64 ! excludes mask.
#' @inheritParams ExtractImages_toMatrix
#' @param mode (\code{\link{objectParam}} argument) color mode export. Either \code{"rgb"}, \code{"gray"}. Default is \code{"rgb"}. 
#' @note Arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractImages_toBase64}} input arguments.\cr
#' @details If \code{'param'} is provided 'in \code{'...'}, \code{'param$export'<-"base64"} and \code{'param$mode'<-'mode'} \strong{only} will be overwritten.
#' @return A list of base64 encoded images corresponding to objects extracted.
#' @export
ExtractImages_toBase64 <- function(...,
                                   objects,
                                   offsets,
                                   display_progress = TRUE,
                                   mode = c("rgb","gray")[1]) {
  dots = list(...)
  dots = dots[!(names(dots) %in% c("image_type", "export"))]
  args = list(image_type = "img", export = "base64")
  if(!missing(objects)) args = c(args, list(objects = objects))
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  args = c(args, list(display_progress = display_progress, mode = mode))
  do.call(what = polyExtractTo, args = c(dots, args))
}


#' @title Shortcut for Batch Masks Extraction to Base64
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of masks to base64 ! excludes image.
#' @inheritParams ExtractImages_toBase64
#' @note Arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractMasks_toBase64}} input arguments.
#' @inherit ExtractImages_toBase64 details return
#' @export
ExtractMasks_toBase64 <- function(...,
                                  objects,
                                  offsets,
                                  display_progress = TRUE,
                                  mode = c("rgb","gray")[1]) {
  dots = list(...)
  dots = dots[!(names(dots) %in% c("image_type", "export"))]
  args = list(image_type = "msk", export = "base64")
  if(!missing(objects)) args = c(args, list(objects = objects))
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  args = c(args, list(display_progress = display_progress, mode = mode))
  do.call(what = polyExtractTo, args = c(dots, args))
}

#' @title Shortcut for Batch Images Extraction to Files
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of images to file ! excludes mask.
#' @param mode (\code{\link{objectParam}} argument) color mode export. Either \code{"rgb"}, \code{"gray"}. Default is \code{"rgb"}.
#' @inheritParams ExtractImages_toBase64
#' @param write_to used to compute respectively exported file name.\cr
#' Exported \code{"file"} extension will be deduced from this pattern. Allowed export are \code{".bmp"}, \code{".jpg"}, \code{".jpeg"}, \code{".png"}, \code{".tif"}, \code{".tiff"}.
#' Note that \code{".bmp"} is faster but files are not compressed producing bigger data.\cr
#' Placeholders, if found, will be substituted:\cr
#' -\code{\%d}: with full path directory\cr
#' -\code{\%p}: with first parent directory\cr
#' -\code{\%e}: with extension (without leading .)\cr
#' -\code{\%s}: with shortname (i.e. basename without extension)\cr
#' -\code{\%o}: with object_id\cr
#' -\code{\%c}: with channel_id\cr
#' A good trick is to use:\cr
#' -\code{"\%d/\%s/\%s_\%o_\%c.tiff"}, when \code{'export'} is \code{"file"}.
#' @note Arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractImages_toFile}} input arguments.
#' @details If \code{'param'} is provided in \code{'...'}, \code{'param$export'<-"file"}, \code{'param$write_to'<-'write_to'} and \code{'param$mode'<-'mode'} \strong{only} will be overwritten.
#' @return It invisibly returns a list of exported file path of corresponding to objects extracted.
#' @export
ExtractImages_toFile <- function(...,
                                 objects,
                                 offsets,
                                 display_progress = TRUE,
                                 mode = c("rgb","gray")[1], 
                                 write_to) {
  dots = list(...)
  dots = dots[!(names(dots) %in% c("image_type", "export"))]
  args = list(image_type = "img", export = "file")
  if(!missing(objects)) args = c(args, list(objects = objects))
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  args = c(args, list(display_progress = display_progress, mode = mode))
  if(!missing(write_to)) args = c(args, list(write_to = write_to))
  do.call(what = polyExtractTo, args = c(dots, args))
}

#' @title Shortcut for Batch Masks Extraction to Files
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of masks to file ! excludes image.
#' @inheritParams ExtractImages_toFile
#' @note Arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractMasks_toFile}} input arguments.
#' @inherit ExtractImages_toFile details return
#' @export
ExtractMasks_toFile <- function(...,
                                objects,
                                offsets,
                                display_progress = TRUE,
                                mode = c("rgb","gray")[1], 
                                write_to) {
  dots = list(...)
  dots = dots[!(names(dots) %in% c("image_type", "export"))]
  args = list(image_type = "msk", export = "file")
  if(!missing(objects)) args = c(args, list(objects = objects))
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  args = c(args, list(display_progress = display_progress, mode = mode))
  if(!missing(write_to)) args = c(args, list(write_to = write_to))
  do.call(what = polyExtractTo, args = c(dots, args))
}
