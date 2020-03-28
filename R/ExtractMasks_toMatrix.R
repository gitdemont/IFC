#' @title Shorcut for Batch Masks Extraction to Matrices/Arrays
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of masks to matrix ! excludes image.
#' @param fileName path of file to read data from. This parameter can't be missing.
#' @param display object of class IFC_display, rich information extracted by \code{\link{getDisplayInfo}}. If missing, the default, 'display' will be extracted from fileName.\cr
#' This param is not mandatory but it may allow to pass custom display info by manually modifying display$Images for example.
#' @param offsets object of class IFC_offsets. If missing, the default, 'offsets' will be extracted from fileName.\cr
#' This param is not mandatory but it may allow to save time for repeated image export on same file.
#' @param objects integers, indices of objects to use.
#' This param is not mandatory, if missing, the default, all objects will be used.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd' and 'export'(="matrix").\cr
#' If 'display' and/or 'offsets' are not provided arguments can also be passed to \code{\link{getDisplayInfo}} and/or \code{\link{getOffsets}}.
#' @details arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractMasks_toMatrix}} input arguments.
#' @return A list of matrices/arrays of masks corresponding to objects extracted.
#' @export
ExtractMasks_toMatrix <- function(fileName, display, offsets, objects, display_progress = TRUE, 
                                  ...) {
  dots=list(...)
  
  # check mandatory param
  if(missing(fileName)) stop("'fileName' can't be missing")
  if(!file.exists(fileName)) stop(paste("can't find",fileName,sep=" "))
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf", "rif", "cif"))
  title_progress = basename(fileName)
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  
  # process extra parameters
  param_infos = names(dots) %in% c("from","verbose","verbosity","warn","force_default")
  param_extra = names(dots) %in% c("ifd","display","export","verbose")
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
  dots_infos = dots[param_infos]
  dots = dots[!param_extra]
  
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
    offsets = suppressMessages(getOffsets(fileName = infos$fileName_image, fast = fast, display_progress = display_progress))
  }
  
  # check objects to extract
  nobj = as.numeric(infos$objcount)
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(length(tokeep) == 0) {
      warning("ExtractMasks_toMatrix: No objects to extract, check the objects you provided.", immediate. = TRUE, call. = FALSE)
      return(NULL)
    }
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }
  sel = split(objects, ceiling(seq_along(objects)/20))
  L=length(sel)
  
  # extract objects
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = 1, initial = 0, style = 3)
    on.exit(endPB(pb))
    ans = lapply(1:L, FUN=function(i) {
      setPB(pb, value = i/L, title = title_progress, label = "exporting masks to matrix")
      do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = infos$fileName_image, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], objects_type = "msk"), trunc_bytes = 8, force_trunc = TRUE, verbose = verbose, verbosity = verbosity),
                                                    display = infos,
                                                    export = "matrix",
                                                    verbose = verbose),
                                               dots))
    })
  } else {
    ans = lapply(1:L, FUN=function(i) {
      do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = infos$fileName_image, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], objects_type = "msk"), trunc_bytes = 8, force_trunc = TRUE, verbose = verbose, verbosity = verbosity),
                                                    display = infos,
                                                    export = "matrix",
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
  ids = as.integer(gsub("^.*_(.*)$", "\\1", sapply(ans, attr, which="offset_id")))
  if(length(ids) != 0) if(!all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  names(ans) = objects
  attr(ans, "channel_id") <- channel_id
  return(ans)
}
