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

#' @title Numpy Export
#' @description
#' Exports IFC objects to numpy array [objects,height,width,channels]
#' @param \dots arguments to be passed to \code{\link{objectExtract}} with the exception of \code{'ifd'} and \code{'bypass'}(=\strong{TRUE}).\cr
#' \strong{/!\\} If not any of \code{'fileName'}, \code{'info'} and \code{'param'} can be found in \code{'...'} then \code{attr(offsets, "fileName_image")} will be used as \code{'fileName'} input parameter to pass to \code{\link{objectParam}}.
#' @param objects integer vector, IDEAS objects ids numbers to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`.
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.\cr
#' If \code{'offsets'} are not provided, extra arguments can also be passed with \code{'...'} to \code{\link{getOffsets}}.
#' @param image_type type of desired object offsets. Either \code{"img"} or \code{"msk"}. Default is \code{"img"}.
#' @param size a length 2 integer vector of final dimensions of the image, height 1st and width 2nd. Default is \code{c(64,64)}.
#' @param force_width whether to use information in \code{'info'} to fill \code{'size'}. Default is \code{FALSE}.
#' When set to \code{TRUE}, width of \code{'size'} argument will be overwritten.
#' @param display_progress whether to display a progress bar. Default is \code{TRUE}.
#' @param python path to python. Default is \code{Sys.getenv("RETICULATE_PYTHON")}.\cr
#' Note that \code{numpy} should be available in this \code{python} to allow export to \code{".npy"} file, otherwise \code{'export'} will be forced to \code{"matrix"}.
#' @param dtype desired arrays data-type. Default is \code{"double"}. Allowed are \code{"uint8"}, \code{"int16"}, \code{"uint16"} or \code{"double"}. If \code{'mode'} is \code{"raw"}, this parameter will be forced to \code{"int16"}.
#' @param mode (\code{\link{objectParam}} argument) color mode export. Either \code{"raw"}, \code{"gray"}. Default is \code{"raw"}.
#' @param export export format. Either \code{"file"}, \code{"matrix"}. Default is \code{"matrix"}.\cr
#' Note that you will need \code{reticulate} package installed to be able to export to \code{".npy"} file, otherwise \code{'export'} will be forced to \code{"matrix"}.
#' @param write_to used when \code{'export'} is \code{"file"} to compute exported file name and type.
#' Exported type will be deduced from this pattern. Allowed exported file extension is \code{".npy"}.\cr
#' Placeholders, if found, will be substituted:\cr
#' -\code{\%d}: with full path directory\cr
#' -\code{\%p}: with first parent directory\cr
#' -\code{\%e}: with extension of (without leading .)\cr
#' -\code{\%s}: with shortname (i.e. basename without extension)\cr
#' -\code{\%o}: with objects (at most 10, will be collapse with "_", if more than one).\cr
#' -\code{\%c}: with channel_id (will be collapse with "_", if more than one, composite if any will be bracketed).
#' A good trick is to use:\cr
#' -\code{"\%d/\%s_Obj[\%o]_Ch[\%c].npy"}, when \code{'export'} is \code{"file"}.
#' @param overwrite whether to overwrite file or not. Default is \code{FALSE}.
#' @note Arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExportToNumpy}} input arguments.
#' @details Please note that \code{'size'} parameter has to be supplied and could not be set to (0,) when \code{'object'} length is not equal to one\cr
#' \code{\link{ExportToNumpy}} requires \code{reticulate} package, \code{python} and \code{numpy} installed to create \code{".npy"} file.\cr
#' If one of these is missing, \code{'export'} will be set to \code{"matrix"}.
#' If \code{'param'} is provided in \code{'...'}, \code{param$export <- "matrix"}, \code{param$mode <- 'mode'} and \code{param$size <- 'size'} and will be overwritten.
#' @return Depending on \code{'export'}:\cr
#' -\code{"matrix"}, an array whose dimensions are [object, height, width, channel].\cr
#' -\code{"file"}, it invisibly returns path of \code{".npy"} exported file. 
#' @export
ExportToNumpy <- function(...,
                          objects,
                          offsets, 
                          image_type = "img", 
                          size = c(64,64),
                          force_width = FALSE,
                          display_progress = TRUE,
                          python = Sys.getenv("RETICULATE_PYTHON"),
                          dtype = c("uint8", "int16", "uint16", "double")[3],
                          mode = c("raw", "gray")[1], 
                          export = c("file", "matrix")[2],
                          write_to, 
                          overwrite = FALSE) {
  dots=list(...)
  
  # check mandatory param
  assert(image_type, len = 1, alw = c("img", "msk"))
  assert(mode, len = 1, alw = c("raw", "gray"))
  if(mode == "raw") dtype = "int16"
  assert(dtype, len = 1, alw = c("uint8", "int16", "uint16", "double"))
  python_back = Sys.getenv("RETICULATE_PYTHON")
  on.exit(Sys.setenv("RETICULATE_PYTHON" = python_back), add = TRUE)
  Sys.setenv("RETICULATE_PYTHON" = python)
  assert(export, len = 1, alw = c("file", "matrix"))
  if(export == "file") {
    if(!requireNamespace("reticulate")) {
      warning("ExportToNumpy: Please install 'reticulate' to export to numpy array file. 'export' has been forced to \"matrix\"")
      export = "matrix"
    } else {
      if(reticulate::py_numpy_available(initialize = TRUE)) {
        np <- reticulate::import("numpy", convert = FALSE)
      } else {
        warning("ExportToNumpy: Can't find numpy in your python installation. 'export' has been forced to \"matrix\"")
        export = "matrix"
      }
    }
  }
  size = na.omit(as.integer(size[1:2]))
  assert(size, len=2, typ="integer")
  force_width = as.logical(force_width); assert(force_width, len = 1, alw = c(TRUE,FALSE)) 
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE,FALSE))
  overwrite = as.logical(overwrite); assert(overwrite, len = 1, alw = c(TRUE,FALSE))
  assert(python, len = 1, typ = "character")
  
  # precompute param
  args=list(mode = mode,
            size = size,
            force_width = force_width,
            export = "matrix",
            overwrite = overwrite)
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  param = do.call(what = dotsParam, args = c(dots, args))
  H = param$size[1]
  W = param$size[2]
  C = length(param$chan_to_keep) + length(param$composite_desc)
  channel_id = c(param$chan_to_keep, param$composite)
  channel_names = c(sapply(param$chan_to_keep, USE.NAMES = FALSE, FUN = function(x) param$channels$name[param$channels$physicalChannel == as.integer(x)]),param$composite)
  
  fileName = param$fileName_image
  title_progress = basename(fileName)
  file_extension = getFileExt(fileName)
  
  # check objects
  nobj = as.integer(param$objcount)
  N = nchar(sprintf("%1.f",abs(nobj-1)))
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

  if(length(objects)!=1) {
    if(W == 0) stop("'size' width should be provided when 'object' length not equal to one")
    if(H == 0) stop("'size' height should be provided when 'object' length not equal to one")
  }
  
  # check export/write_to
  overwritten = FALSE
  if(export != "matrix") {
    if(missing(write_to)) {
      if(export == "file") stop("'write_to' can't be missing")
      write_to = "%s_numpy.npy"
    }
    assert(write_to, len = 1, typ = "character")
    type = getFileExt(write_to)
    assert(type, len = 1, alw = "npy")
    splitf_obj = splitf(param$fileName_image)
    splitp_obj = splitp(write_to)
    if(length(objects) > 10) {
      obj_text = paste0(sprintf(paste0("%0",N,".f"), objects[1:10]), collapse="_")
      obj_text = paste0(obj_text, "_...")
    } else {
      obj_text = paste0(sprintf(paste0("%0",N,".f"), objects), collapse="_")
    }
    chan_text = paste0(sprintf(paste0("%0",2,".f"), as.integer(param$chan_to_keep)), collapse="_")
    if(length(param$composite) != 0) {
      comp_text = paste0("_(",gsub("/",",",param$composite),")", collapse="_")
    } else {
      comp_text = ""
    }
    write_to = formatn(splitp_obj, 
                       splitf_obj, 
                       object = obj_text,
                       channel = paste0(chan_text,comp_text))
    if(export == "file") {
      dir_name = dirname(write_to)
      if(!dir.exists(dir_name)) if(!dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dir_name))
      if(file.exists(write_to)) {
        if(!overwrite) stop(paste0("file ", write_to, " already exists"))
        write_to = normalizePath(write_to, winslash = "/")
        overwritten = TRUE
      }
      message(paste0("file will be exported in :\n", normalizePath(dirname(write_to), winslash = "/")))
    }
  }
  
  # extract images/masks
  dots=dots[setdiff(names(dots), c("param","mode","objects","display_progress"))]
  args = list(param = param, mode = mode, objects = objects, display_progress = display_progress)
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  if(image_type == "img") { fun = ExtractImages_toMatrix } else { fun = ExtractMasks_toMatrix }
  ans = do.call(what = fun, args = c(dots, args))
  L = length(ans)
  if(L == 0) {
    if(export == "file") {
      warning(paste0("ExportToNumpy: No objects to export, check the objects you provided.\n",
                     "Can't create 'write_to' =", write_to, " from file.\n", param$fileName_image),
              immediate. = TRUE, call. = FALSE)
      return(invisible(structure(character(), object_id = structure(numeric(), names = character()))))
    } else {
      warning("ExportToNumpy: No objects to export, check the objects you provided.\n", immediate. = TRUE, call. = FALSE)
      return(structure(array(NA_real_, dim = c(0,H,W,C), list("object" = character(0),
                                                              "height" = NULL,
                                                              "width" = NULL,
                                                              "channel" = channel_id)),
                       object_id = structure(numeric(), names = character()),
                       offset_id = structure(numeric(), names = character()),
                       channel_id = channel_id,
                       channel_names = channel_names))
    }
  }
  
  # check object_ids
  if(image_type == "img") { 
    ids = sapply(ans, attr, which = "object_id")
  } else {
    ids = as.integer(gsub("msk_", "", gsub("img_", "", sapply(ans, attr, which = "offset_id"), fixed = TRUE), fixed = TRUE))
  }
  if(!all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  if(length(objects) >= 1) {
    H = nrow(ans[[1]][[1]])
    W = ncol(ans[[1]][[1]])
    C = length(ans[[1]])
  }
  
  # create array [obj,height,width,channel]
  switch(dtype,
         "uint8" = {
           ret = aperm(array(as.integer(unlist(ans, use.names = FALSE) * 255), dim = c(H, W, C, length(objects))), perm = c(4,1,2,3))
         },
         "int16" = {
           if(mode == "raw") {
             ret = aperm(array(as.integer(unlist(ans)), dim = c(H, W, C, length(objects))), perm = c(4,1,2,3))
           } else {
             ret = aperm(array(as.integer(unlist(ans, use.names = FALSE) * 32767), dim = c(H, W, C, length(objects))), perm = c(4,1,2,3))
           }
         },
         "uint16" = {
           ret = aperm(array(as.integer(unlist(ans, use.names = FALSE) * 65535), dim = c(H, W, C, length(objects))), perm = c(4,1,2,3))
         },
         {
           ret = aperm(array(unlist(ans, use.names = FALSE), dim = c(H, W, C, length(objects))), perm = c(4,1,2,3))
         })
  dimnames(ret) = list("object" = num_to_string(ids),
                       "height" = NULL,
                       "width" = NULL,
                       "channel" = channel_id)
  attr(ret, "fileName_image") <- param$fileName_image
  attr(ret, "object_id") <- ids
  attr(ret, "offset_id") <- sapply(ans, attr, which = "offset_id")
  attr(ret, "channel_id") <- channel_id
  attr(ret, "channel_names") <- channel_names
  tryCatch({
    if(export == "file") {
      np$save(file = write_to, arr = np$array(ret, dtype=np[[dtype]], order='C'))
      message(paste0("\n######################\n", normalizePath(write_to, winslash = "/", mustWork = FALSE), "\nhas been successfully ", ifelse(overwritten, "overwritten", "exported"), "\n"))
      attr(write_to, "fileName_image") <- param$fileName_image
      attr(write_to, "object_id") <- ids
      attr(write_to, "offset_id") <- sapply(ans, attr, which = "offset_id")
      attr(write_to, "channel_id") <- channel_id
      attr(write_to, "channel_names") <- channel_names
      return(invisible(write_to))
    }
    return(ret)
  }, error = function(e) {
    stop(ifelse(export == "file", paste0(normalizePath(write_to, winslash = "/", mustWork = FALSE), " has been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n", e$message), e$message), call. = FALSE)
  })
}
