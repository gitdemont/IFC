#' @title Numpy Export
#' @description
#' Exports IFC objects to numpy array [objects,height,width,channels]
#' @param display object of class IFC_display, rich information extracted by \code{\link{getDisplayInfo}}. This parameter can't be missing.
#' @param offsets object of class IFC_offsets. If missing, the default, offsets will be extracted from display$fileName.\cr
#' This param is not mandatory but it may allow to save time for repeated image export on same file.
#' @param objects integers, indices of objects to use.
#' This param is not mandatory, if missing, the default, all objects will be used.
#' @param objects_type objects_type of desired offsets. Either "img" or "msk". Default is "img".
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param python path to python. Default is Sys.getenv("RETICULATE_PYTHON").\cr
#' Note that this numpy should be available in this python to be able to export to numpy array file, otherwise 'export' will be forced to "matrix".
#' @param dtype desired array’s data-type. Default is "double". Allowed are "uint8", "int16", "uint16" or "double". If 'mode' is "raw", this parameter will be forced to "int16".
#' @param mode (\code{\link{objectExtract}} argument) color mode export. Either "raw", "gray" . Default is "gray".
#' @param export export format. Either "file", "matrix", "base64". Default is "matrix".\cr
#' Note that you will need 'reticulate' package installed to be able to export to numpy array file, otherwise 'export' will be forced to "matrix".
#' @param export_to used when 'export' is "file" to compute respectively filename.
#' Exported type will be deduced from this pattern. Allowed export are '.npy'.\cr
#' Placeholders, if found, will be substituted:\cr
#' -\%d: with full path directory of 'display$fileName_image'\cr
#' -\%p: with first parent directory of ''display$fileName_image'\cr
#' -\%e: with extension of 'display$fileName_image' (without leading .)\cr
#' -\%s: with shortname from 'display$fileName_image' (i.e. basename without extension)\cr
#' -\%o: with objects (at most 10, will be collapse with "_", if more than one).\cr
#' -\%c: with channel_id (will be collapse with "_", if more than one, composite in any will be bracketed).
#' A good trick is to use:\cr
#' -"\%d/\%s_Obj[\%o]_Ch[\%c].npy", when 'export' is "file"\cr
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' @param ... other arguments to be passed to \code{\link{objectExtract}} with the exception of 'ifd'.
#' If 'offsets' are not provided arguments can also be passed to \code{\link{getOffsets}}.
#' @details arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExportToNumpy}} input arguments.\cr
#' \code{\link{ExportToNumpy}} requires reticulate package, python and numpy installed. to create npy file.\cr
#' If one of these is missing, 'export' will be set to "matrix".
#' @return Depending on 'export':\cr
#' -"matrix", a numpy array,\cr
#' -"file", an invisible vector of ids corresponding to the objects exported. 
#' @export
ExportToNumpy <- function(display, offsets, objects, objects_type = "img", display_progress = TRUE,
                          python = Sys.getenv("RETICULATE_PYTHON"),
                          dtype = c("uint8", "int16", "uint16", "double")[3],
                          mode = c("raw", "gray")[1], 
                          export = c("file", "matrix")[2],
                          export_to, overwrite = FALSE,
                          ...) {
  dots=list(...)
  # check mandatory param
  if(missing(display)) stop("'display' can't be missing")
  assert(display, cla = "IFC_display")
  assert(objects_type, len = 1, alw = c("img", "msk"))
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
      if(!reticulate::py_module_available("numpy")) {
        warning("ExportToNumpy: Can't find numpy in your python installation. 'export' has been forced to \"matrix\"")
        export = "matrix"
      } else {
        np <- reticulate::import("numpy", convert  = FALSE)
      }
    }
  }
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE,FALSE))
  overwrite = as.logical(overwrite); assert(overwrite, len = 1, alw = c(TRUE,FALSE))
  assert(python, len = 1, typ = "character")
  
  fileName = display$fileName
  title_progress = basename(fileName)
  file_extension = getFileExt(fileName)
  channels = display$Images[display$Images$physicalChannel %in% which(display$in_use), ]
  
  # process extra parameters
  param_extra = names(dots) %in% c("ifd","export","mode","size","force_width","verbose")
  param_param = names(dots) %in% c("composite","selection","random_seed",
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
  
  # should be checked before being passed to objectParam/objectExtract
  if(length(dots[["size"]]) == 0) {
    size = c(0,0)
  } else {
    size = dots[["size"]]
  }
  # should be checked before being passed to objectParam/objectExtract
  if(length(dots[["force_width"]]) == 0) {
    force_width = TRUE
  } else {
    force_width = dots[["force_width"]]
  }
  dots = dots[!param_extra] # remove not allowed param
  dots_param = dots[param_param] # keep param_param for objectParam
  dots = dots[!param_param]
  
  # check objects to extract
  nobj = as.integer(display$objcount)
  N = nchar(sprintf("%1.f",abs(nobj-1)))
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(length(tokeep) == 0) {
      if(export == "file") {
        warning(paste0("ExportToNumpy: No objects to export, check the objects you provided.\n",
                       "Can't create 'export_to' =", export_to, " from file.\n", display$fileName_image),
                immediate. = TRUE, call. = FALSE)
        return(invisible(NULL))
      } else {
        warning("ExportToNumpy: No objects to export, check the objects you provided.\n", immediate. = TRUE, call. = FALSE)
        return(NULL)
      }
    }
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }
  
  # check size
  force_width = as.logical(force_width); assert(force_width, len = 1, alw = c(TRUE,FALSE)) 
  size = na.omit(as.integer(size[1:2]))
  assert(size, len=2, typ="integer")
  if(!force_width) {
    if(length(objects)!=1) if(size[2] == 0) stop("'size' width should be provided when 'force_width' is set to FALSE and 'objects' length not equal to one")
  } else {
    size = c(size[1], as.integer(display$channelwidth))
  }
  
  # check input offsets if any
  compute_offsets = TRUE
  if(!missing(offsets)) {
    if(!("IFC_offsets" %in% class(offsets))) {
      warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
    } else {
      if(attr(offsets, "checksum") != display$checksum) {
        warning("provided 'offsets' do not match with expected ones, 'offsets' will be recomputed", immediate. = TRUE, call. = FALSE)
      } else {
        compute_offsets = FALSE
      }
    }
  }
  if(compute_offsets) {
    offsets = suppressMessages(getOffsets(fileName = display$fileName_image, fast = fast, display_progress = display_progress, verbose = verbose))
  }
  
  # compute object param
  is_param = names(dots) %in% "param"
  if(any(is_param)) {
    param = dots$param
    param$size = size
    dots = dots[!is_param]
  } else {
    param = do.call(what = "objectParam",
                    args = c(list(display = display, 
                                  size = size, 
                                  force_width = force_width), dots))
  }
  
  # check export/export_to
  overwritten = FALSE
  if(export != "matrix") {
    if(missing(export_to)) {
      if(export == "file") stop("'export_to' can't be missing")
      export_to = "%s_numpy.npy"
    }
    assert(export_to, len = 1, typ = "character")
    type = getFileExt(export_to)
    assert(type, len = 1, alw = "npy")
    splitf_obj = splitf(display$fileName_image)
    splitp_obj = splitp(export_to)
    if(length(objects) > 10) {
      obj_text = paste0(sprintf(paste0("%0",N,".f"), objects[1:10]), collapse="_")
      obj_text = paste0(obj_text, "_...")
    } else {
      obj_text = paste0(sprintf(paste0("%0",N,".f"), objects), collapse="_")
    }
    export_to = formatn(splitp_obj, splitf_obj, object = obj_text,
                        channel = paste0(paste0(sprintf(paste0("%0",2,".f"), param$chan_to_keep), collapse="_"),"_",paste0("(",gsub("/",",",param$composite),")", collapse="_")))
    if(export == "file") {
      dir_name = dirname(export_to)
      if(!dir.exists(dir_name)) if(!dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dir_name))
      if(file.exists(export_to)) {
        if(!overwrite) stop(paste0("file ", export_to, " already exists"))
        export_to = normalizePath(export_to, winslash = "/")
        overwritten = TRUE
      }
      message(paste0("file will be exported in :\n", normalizePath(dirname(export_to), winslash = "/")))
    }
  }
  
  # extract objects
  sel = split(objects, ceiling(seq_along(objects)/20))
  L = length(sel)
  tryCatch({
    if(display_progress) {
      pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
      ans = lapply(1:L, FUN = function(i) {
        setPB(pb, value = i, title = title_progress, label = "exporting objects to numpy")
        do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = display$fileName_image, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], objects_type = objects_type), trunc_bytes = 8, force_trunc = FALSE, verbose = verbose, verbosity = verbosity), 
                                                      display = display, 
                                                      param = param,
                                                      export = "matrix", 
                                                      mode = mode,
                                                      verbose = verbose),
                                                 dots))
      })
    } else {
      ans = lapply(1:L, FUN = function(i) {
        do.call(what = "objectExtract", args = c(list(ifd = getIFD(fileName = display$fileName_image, offsets = subsetOffsets(offsets = offsets, objects = sel[[i]], objects_type = objects_type), trunc_bytes = 8, force_trunc = FALSE, verbose = verbose, verbosity = verbosity), 
                                                      display = display, 
                                                      param = param,
                                                      export = "matrix", 
                                                      mode = mode,
                                                      verbose = verbose),
                                                 dots))
      })
    }
  }, error = function(e) {
    stop(e$message, call. = FALSE)
  }, finally = {
    if(display_progress) endPB(pb)
  })
  channel_id = attr(ans[[1]][[1]], "channel_id")
  if (L > 1) {
    ans = do.call(what = "c", args = ans)
  } else {
    ans = ans[[1]]
  }
  channel_names = names(ans[[1]])
  
  # check object_ids
  if(objects_type == "img") { 
    ids = sapply(ans, attr, which = "object_id")
  } else {
    ids = as.integer(gsub("^.*_(.*)$", "\\1", sapply(ans, attr, which = "offset_id")))
  }
  if(!all(objects == ids)) warning("Extracted object_ids differ from expected ones. Concider running with 'fast' = FALSE", call. = FALSE, immediate. = TRUE)
  
  # create array [obj,height,width,channel]
  switch(dtype,
         "uint8" = {
           ret = aperm(array(as.integer(unlist(ans) * 255), dim = c(nrow(ans[[1]][[1]]), ncol(ans[[1]][[1]]), length(ans[[1]]), length(objects))), perm = c(4,1,2,3))
         },
         "uint16" = {
           ret = aperm(array(as.integer(unlist(ans) * 65535), dim = c(nrow(ans[[1]][[1]]), ncol(ans[[1]][[1]]), length(ans[[1]]), length(objects))), perm = c(4,1,2,3))
         },
         {
           ret = aperm(array(unlist(ans), dim = c(nrow(ans[[1]][[1]]), ncol(ans[[1]][[1]]), length(ans[[1]]), length(objects))), perm = c(4,1,2,3))
         })
  attr(ret, "object_id") <- ids
  attr(ret, "channel_id") <- channel_id
  attr(ret, "channel_names") <- channel_names
  tryCatch({
    if(export == "file") {
      np$save(file = export_to, arr = np$array(ret, dtype=np[[dtype]], order='C'))
      message(paste0("\n######################\n", normalizePath(export_to, winslash = "/", mustWork = FALSE), "\nhas been successfully ", ifelse(overwritten, "overwritten", "exported"), "\n"))
      attr(ids, "filename") <- export_to
      attr(ids, "channel_id") <- channel_id
      attr(ids, "channel_names") <- channel_names
      return(invisible(ids))
    }
    return(ret)
  }, error = function(e) {
    stop(ifelse(export == "file", paste0(normalizePath(export_to, winslash = "/", mustWork = FALSE), " has been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n", e$message), e$message), call. = FALSE)
  })
}
