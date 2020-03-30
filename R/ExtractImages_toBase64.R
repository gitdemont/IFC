#' @title Shorcut for Batch Images Extraction to Base64
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of images to matrix ! excludes mask.
#' @param fileName path of file to read data from. This parameter can't be missing.
#' @param display object of class IFC_display, rich information extracted by \code{\link{getDisplayInfo}}. If missing, the default, 'display' will be extracted from fileName.\cr
#' This param is not mandatory but it may allow to pass custom display info by manually modifying display$Images for example.
#' @param offsets object of class IFC_offsets. If missing, the default, offsets will be extracted from fileName.\cr
#' This param is not mandatory but it may allow to save time for repeated image export on same file.
#' @param objects integers, indices of objects to use.
#' This param is not mandatory, if missing, the default, all objects will be used.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param mode (\code{\link{objectExtract}} argument) color mode export. Either "rgb", "gray". Default is "rgb".
#' @param ... other arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd', 'export'(="base64") and 'mode'.\cr
#' If 'display' and/or 'offsets' are not provided arguments can also be passed to \code{\link{getDisplayInfo}} and/or \code{\link{getOffsets}}.
#' @details arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractImages_toBase64}} input arguments.
#' @return A list of base64 encoded images corresponding to objects extracted.
#' @export
ExtractImages_toBase64 <- function(fileName, display, offsets, objects, display_progress = TRUE,
                                   mode = c("rgb","gray")[1], ...) {
  dots=list(...)
  
  # check mandatory param
  if(missing(fileName)) stop("'fileName' can't be missing")
  if(!file.exists(fileName)) stop(paste("can't find",fileName,sep=" "))
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf", "rif", "cif"))
  title_progress = basename(fileName)
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  assert(mode, len = 1, alw = c("rgb", "gray"))
  
  # process extra parameters
  param_infos = names(dots) %in% c("from","verbose","verbosity","warn","force_default")
  param_extra = names(dots) %in% c("ifd","display","export","mode","verbose")
  param_param = names(dots) %in% c("composite","selection","size","force_width","random_seed",
                                   "removal","add_noise","full_range","force_range")
  if(length(dots[["verbose"]]) == 0) { # param for objectExtract, getDisplayInfo, getIFD, getOffsets
    verbose = FALSE
  } else {
    verbose = dots[["verbose"]]
  }
  if(length(dots[["verbosity"]]) == 0) { # param for getDisplayInfo, getIFD
    verbosity = 1
  } else {
    verbosity = dots[["verbosity"]]
  }
  if(length(dots[["fast"]]) == 0) { # param for getOffsets
    fast = TRUE
  } else {
    fast = dots[["fast"]]
  }
  dots_infos = dots[param_infos] # keep param_infos fo getDisplayInfo
  dots = dots[!param_extra] # remove not allowed param
  dots_param = dots[param_param] # keep param_param for objectParam
  dots = dots[!param_param]
  
  # check input display if any
  infos = do.call(what = "getDisplayInfo", args = c(list(fileName=fileName), dots_infos))
  if(!missing(display)) {
    if(!("IFC_display" %in% class(display))) {
      warning("provided 'display' do not match with expected one, 'display' will not be used", immediate. = TRUE, call. = FALSE)
    } else {
      if(display$checksum == checksumXIF(infos$fileName_image)) {
        infos = display
      } else {
        warning("provided 'display' do not match with expected one, 'display' will not be used", immediate. = TRUE, call. = FALSE)
      }
    }
  }
  
  # check input offsets if any
  compute_offsets = TRUE
  if(!missing(offsets)) {
    if(!("IFC_offsets" %in% class(offsets))) {
      warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
    } else {
      if(attr(offsets, "checksum") != checksumXIF(infos$fileName_image)) {
        warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
      } else {
        compute_offsets = FALSE
      }
    }
  }
  if(compute_offsets) {
    offsets = suppressMessages(getOffsets(fileName = infos$fileName_image, fast = fast, display_progress = display_progress, verbose = verbose))
  }
  
  # check objects to extract
  nobj = as.numeric(infos$objcount)
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(length(tokeep) == 0) {
      warning("ExtractImages_toBase64: No objects to extract, check the objects you provided.", immediate. = TRUE, call. = FALSE)
      return(NULL)
    }
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }
  sel = split(objects, ceiling(seq_along(objects)/20))
  L=length(sel)
  
  # compute object param
  is_param = names(dots) %in% "param"
  if(any(is_param)) {
    param = dots$param
    dots = dots[!is_param]
  } else {
    param = do.call(what = "objectParam", args = c(list(display = infos), dots))
  }
  
  # extract objects
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
    ans = lapply(1:L, FUN=function(i) {
      setPB(pb, value = i, title = title_progress, label = "exporting images to base64")
      do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = infos$fileName_image, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], objects_type = "img"), trunc_bytes = 8, force_trunc = TRUE, verbose = verbose, verbosity = verbosity),
                                                    display = infos,
                                                    param = param,
                                                    export = "base64",
                                                    mode = mode,
                                                    verbose = verbose),
                                               dots))
    })
  } else{
    ans = lapply(1:L, FUN=function(i) {
      do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = infos$fileName_image, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], objects_type = "img"), trunc_bytes = 8, force_trunc = TRUE, verbose = verbose, verbosity = verbosity),
                                                    display = infos,
                                                    param = param,
                                                    export = "base64",
                                                    mode = mode,
                                                    verbose = verbose),
                                               dots))
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
