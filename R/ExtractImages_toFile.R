#' @title Shorcut for Batch Images Extraction to Files
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
#' @param mode (\code{\link{objectExtract}} argument) color mode export. Either "rgb", "gray" . Default is "rgb".
#' @param export_to (\code{\link{objectExtract}} argument) used to compute exported file name.\cr
#' Exported "file" extension will be deduced from this pattern. Allowed export are '.bmp', '.jpg', '.jpeg', '.png', '.tif', '.tiff'.
#' Note that '.bmp' are faster but not compressed producing bigger data.\cr
#' Placeholders, if found, will be substituted:\cr
#' -\%d: with full path directory of 'fileName'\cr
#' -\%p: with first parent directory of 'fileName'\cr
#' -\%e: with extension of 'fileName' (without leading .)\cr
#' -\%s: with shortname from 'fileName' (i.e. basename without extension)\cr
#' -\%o: with object_id\cr
#' -\%c: with channel_id\cr
#' A good trick is to use "\%d/\%s/\%s_\%o_\%c.tiff".
#' @param ... other arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd', 'export'(="file"), 'export_to' and 'mode'.\cr
#' If 'display' and/or 'offsets' are not provided arguments can also be passed to \code{\link{getDisplayInfo}} and/or \code{\link{getOffsets}}.
#' @details arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractImages_toFile}} input arguments.
#' @return It invisibly returns a list of exported file path of corresponding to objects extracted.
#' @export
ExtractImages_toFile <- function(fileName, display, offsets, objects, display_progress = TRUE,
                                 mode = c("rgb","gray")[1], export_to, ...) {
  dots=list(...)
  
  # check mandatory param
  if(missing(fileName)) stop("'fileName' can't be missing")
  if(!file.exists(fileName)) stop(paste("can't find",fileName,sep=" "))
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf", "rif", "cif"))
  title_progress = basename(fileName)
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  if(missing(export_to)) stop("'export_to' can't be missing")
  export_to = na.omit(as.character(export_to))
  assert(export_to, len = 1, typ = "character")
  assert(mode, len = 1, alw = c("rgb", "gray"))
  
  # process extra parameters
  param_infos = names(dots) %in% c("from","verbose","verbosity","warn","force_default")
  param_extra = names(dots) %in% c("ifd","display","export","export_to","mode","verbose")
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
      warning("ExtractImages_toFile: No objects to extract, check the objects you provided.", immediate. = TRUE, call. = FALSE)
      return(invisible(NULL))
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
    pb = newPB(min = 0, max = 1, initial = 0, style = 3)
    on.exit(endPB(pb))
    ans = lapply(1:L, FUN=function(i) {
      setPB(pb, value = i/L, title = title_progress, label = "exporting images to file")
      do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = infos$fileName_image, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], objects_type = "img"), trunc_bytes = 8, force_trunc = TRUE, verbose = verbose, verbosity = verbosity),
                                                    display = infos,
                                                    export = "file",
                                                    export_to = export_to,
                                                    mode = mode,
                                                    verbose = verbose),
                                               dots))
    })
  } else{
    ans = lapply(1:L, FUN=function(i) {
      do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = infos$fileName_image, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], objects_type = "img"), trunc_bytes = 8, force_trunc = TRUE, verbose = verbose, verbosity = verbosity),
                                                    display = infos,
                                                    export = "file",
                                                    export_to = export_to,
                                                    mode = mode,
                                                    verbose = verbose),
                                               dots))
    })
  }
  if(L>1) {
    ans = do.call(what="c", args=ans)
  } else {
    ans = ans[[1]]
  }
  ids = sapply(ans, attr, which="object_id")
  if(!all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  return(invisible(ans))
}
