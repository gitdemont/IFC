#' @title Shorcut for Batch Masks Extraction to Matrices/Arrays
#' @description
#' Function to shortcut extraction, normalization and eventually colorization of masks to matrix ! excludes image.
#' @param offsets object of class `IFC_offset`. 
#' This param is not mandatory but it may allow to save time for repeated image export on same file.
#' @param objects integers, indices of objects to use.
#' This param is not mandatory, if missing, the default, all objects will be used.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd', 'export'(="matrix") and bypass(=TRUE).\cr
#' If 'offsets' are not provided arguments can also be passed to \code{\link{getOffsets}}.\cr
#' /!\ If not any of 'fileName', 'info' and 'param' can be found in ... then attr(offsets, "fileName_image") will be used as 'fileName' input parameter to pass to \code{\link{objectParam}}.
#' @details arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExtractMasks_toMatrix}} input arguments.
#' @return A list of matrices/arrays of masks corresponding to objects extracted.
#' @export
ExtractMasks_toMatrix <- function(offsets,
                                  objects, 
                                  display_progress = TRUE,
                                  ...) {
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
  param_extra = names(dots) %in% c("ifd","export")
  dots = dots[!param_extra] # remove not allowed param
  param_param = names(dots) %in% c("write_to","mode","base64_id","base64_att","overwrite",
                                   "composite","selection","random_seed","size","force_width",
                                   "removal","add_noise","full_range","force_range")
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
                                    export = "matrix"), dots_param))
    } else {
      param = do.call(what = "objectParam",
                      args = c(list(info = input$info,
                                    export = "matrix"), dots_param))
    }
  } else {
    param = dots$param
  }
  fileName = param$fileName
  title_progress = basename(fileName)

  # check objects to extract
  nobj = as.numeric(param$objcount)
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
  
  # check input offsets if any
  compute_offsets = TRUE
  if(!missing(offsets)) {
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
  
  # extract objects
  sel = split(objects, ceiling(seq_along(objects)/20))
  L=length(sel)
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
    ans = lapply(1:L, FUN=function(i) {
      setPB(pb, value = i, title = title_progress, label = "exporting masks to matrix")
      do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = param$fileName_image,
                                                                 offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], image_type = "msk"),
                                                                 trunc_bytes = 8, 
                                                                 force_trunc = FALSE, 
                                                                 verbose = verbose, 
                                                                 verbosity = verbosity,
                                                                 bypass = TRUE),
                                                    param = param,
                                                    verbose = verbose,
                                                    bypass = TRUE),
                                               dots))
    })
  } else {
    ans = lapply(1:L, FUN=function(i) {
      do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = param$fileName_image,
                                                                 offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], image_type = "msk"),
                                                                 trunc_bytes = 8, 
                                                                 force_trunc = FALSE, 
                                                                 verbose = verbose, 
                                                                 verbosity = verbosity,
                                                                 bypass = TRUE), 
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
  ids = as.integer(gsub("^.*_(.*)$", "\\1", sapply(ans, attr, which="offset_id")))
  if(length(ids) != 0) if(!all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  names(ans) = objects
  attr(ans, "channel_id") <- channel_id
  return(ans)
}
