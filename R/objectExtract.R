#' @title IFC_object Extraction
#' @description
#' Extracts / Decompress objects stored in RIF or CIF Files.
#' @param ifd list of sub elements of IFD data information extracted by \code{\link{getIFD}}. This parameter can't be missing.
#' @param display object of class IFC_display, information extracted by \code{\link{getDisplayInfo}}. This parameter can't be missing.
#' @param param object of class IFC_param, containing extraction parameters defined by \code{\link{objectParam}}.\cr
#' If this parameter is missing, \code{\link{objectExtract}} will use extra ... to pass arguments to \code{\link{objectParam}} to control object extraction.\cr
#' However, if provided, ... will be ignored.
#' @param mode color mode export. Either "rgb", "gray" or "raw". Default is "raw".
#' Note that "raw" is only possible when 'export' is "matrix".
#' @param export format mode export. Either "file", "matrix", "base64". Default is "matrix".
#' @param export_to Used when export is "file" or "base64" to compute respectively exported file name or base64 id attribute.\cr
#' Exported "file" extension and "base64" MIME type will be deduced from this pattern. Allowed export are ".bmp", ".jpg", ".jpeg", ".png", ".tif", ".tiff".
#' Note that '.bmp' are faster but not compressed producing bigger data.\cr
#' Placeholders, if found, will be substituted:\cr
#' -\%d: with full path directory of 'display$fileName_image'\cr
#' -\%p: with first parent directory of 'display$fileName_image'\cr
#' -\%e: with extension of 'display$fileName_image' (without leading .)\cr
#' -\%s: with shortname from 'display$fileName_image' (i.e. basename without extension)\cr
#' -\%o: with object_id\cr
#' -\%c: with channel_id\cr
#' A good trick is to use:\cr
#' -"\%d/\%s/\%s_\%o_\%c.tiff", when 'export' is "file"\cr
#' -"\%o_\%c.bmp", when 'export' is "base64".\cr
#' Note that if missing and 'export' is not "file", 'export_to' will be set to "\%o_\%c.bmp".
#' @param base64_id whether to add id attribute to base64 exported object. Default is FALSE.\cr
#' Only applied when export is "base64".
#' @param base64_att attributes to add to base64 exported object. Default is "".\cr
#' Only applied when export is "base64". For example, use "class=draggable".\cr
#' Note that id (if base64_id is TRUE) and width and height are already used.
#' @param overwrite only apply when 'export' is "file" whether to overwrite file or not. Default is FALSE.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @param bypass logical to avoid several checking. Default is FALSE.
#' @param ... other arguments to be passed \code{\link{objectParam}}.
#' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/openmicroscopy/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
#' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}\cr\cr
#' BSD implementations of Bio-Formats readers and writers\cr
#' %%
#' Copyright (C) 2005 - 2017 Open Microscopy Environment:\cr
#'   - Board of Regents of the University of Wisconsin-Madison\cr
#'   - Glencoe Software, Inc.\cr
#'   - University of Dundee\cr
#' %%
#' Redistribution and use in source and binary forms, with or without
#' modification, are permitted provided that the following conditions are met:\cr
#' 1. Redistributions of source code must retain the above copyright notice,
#'    this list of conditions and the following disclaimer.\cr
#' 2. Redistributions in binary form must reproduce the above copyright notice,
#'     this list of conditions and the following disclaimer in the documentation
#'     and/or other materials provided with the distribution.\cr
#' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
#' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#' POSSIBILITY OF SUCH DAMAGE.
#' @details when a mask is detected, add_noise and force_range are set to FALSE.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   cif_offs <- getOffsets(fileName = file_cif, fast = TRUE)
#'   ## extract display infomation
#'   disp <- getDisplayInfo(fileName = file_cif, from = "analysis")
#'   ## retrieve number of objects stored
#'   nobj <- as.integer(disp$objcount)
#'   ## randomly subset the offsets of at most 5 "img" objects
#'   sel = sample(0:(nobj-1), min(5, nobj))
#'   sub_offs <- subsetOffsets(cif_offs, objects = sel, objects_type = "img")
#'   ## read IFDs from these "img" objects
#'   IFDs <- getIFD(fileName = file_cif, offsets = sub_offs)
#'   ## extract raw data of these"img" objects to matrix
#'   raw = objectExtract(ifd = IFDs, display = disp, mode = "raw", 
#'                       export = "matrix")
#'   ## extract base64 "rgb" colorized version of these "img" objects to base64
#'   b64 = objectExtract(ifd = IFDs, display = disp, mode = "rgb", 
#'                       export = "base64", base64_id = TRUE,
#'                       export_to = "example_%o_%c.png")
#'   ## use DisplayGallery to show the first "img" objects and play with ... extra parameters
#'   ## force_range, add_noise, selection, composite, see objectParam
#'   DisplayGallery(display = disp, offsets = cif_offs, objects = sel,
#'                  base64_id = TRUE, export_to = "example_%o_%c.png",
#'                  force_range = c(FALSE,TRUE,FALSE,TRUE), add_noise = FALSE,
#'                  selection = c(1,2,4,6), composite = "1.7/4.3")
#' } else {
#'   message(sprintf('Please type `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return A list (for every extracted objects) of list (for every exported channels) depending on "export" parameter:\cr
#' -"matrix", a matrix when 'mode' is set to "raw" or "gray" OR an array when 'mode' == "rgb",\cr
#' -"base64", a data-uri string,\cr
#' -"file", an invisible file path corresponding to the location of exported file(s). 
#' @export
objectExtract <- function(ifd, display, param,
                          mode = c("rgb", "gray", "raw")[3],
                          export = c("file", "matrix", "base64")[2],
                          export_to, base64_id = FALSE, base64_att = "",
                          overwrite = FALSE, verbose = FALSE, bypass = FALSE, ...) {
  dots=list(...)
  bypass = as.logical(bypass); assert(bypass, len = 1, alw = c(TRUE,FALSE))
  # if(!bypass) {
    if(missing(ifd)) stop("'ifd' can't be missing") 
    if(missing(display)) stop("'display' can't be missing") 
    assert(ifd, cla = "IFC_ifd_list")
    if(length(ifd) == 0) stop("can't extract object when 'ifd is empty")
    if("IFC_first_ifd" %in% class(ifd)) stop("can't extract object from 'ifd' of class \"IFC_first_ifd\"")
    assert(display, cla = "IFC_display") 
    if(attr(ifd, "checksum") != display$checksum) stop("'ifd' and 'display' do not match, please ensure that they originate from same file")
  # }
  if(missing(param)) {
    param = do.call(what = "objectParam", args = c(list(display = display, bypass = bypass), dots))
  } else {
    assert(param, cla = "IFC_param")
    if(display$checksum != param$checksum) stop("'param' and 'display' do not match, please ensure that they originate from same file")
  }
  assert(export, len = 1, alw = c("file", "matrix", "base64")) 
  assert(mode, len = 1, alw = c("rgb", "gray", "raw")) 

  fname = display$fileName_image
  if(export != "matrix") { # file or base64
    if(mode == "raw") stop("can't export as \"raw\" when '", export, "' is choosen")
    if(export == "file") { # file
      if(missing(export_to)) stop("'export_to' can't be missing when 'export' is \"file\"")
      # determines type from export_to
      type = getFileExt(export_to)
      switch(type,
             "jpg" = {type <- "jpeg"},
             "tif" = {type <- "tiff"}
      )
      splitf_obj = splitf(fname)
      splitp_obj = splitp(export_to)
      # create folder for file export
      dir_name = dirname(formatn(splitp_obj, splitf_obj))
      if(!dir.exists(dir_name)) if(!dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dir_name))
    } else { # base64
      base64_id = as.logical(base64_id)
      assert(base64_id, len = 1, alw = c(TRUE,FALSE))
      base64_att = na.omit(as.character(base64_att))
      assert(base64_att, len = 1, typ = "character")
      if(missing(export_to)) {
        export_to = "%o_%c.bmp"
        type = "bmp"
      } else {
        type = getFileExt(export_to)
        switch(type,
               "jpg" = {type <- "jpeg"},
               "tif" = {type <- "tiff"}
        )
      }
      if(base64_id) {
        splitf_obj = splitf(fname)
        splitp_obj = splitp(export_to)
      }
    }
  }
  
  # shortcut
  chan_to_extract = param$chan_to_extract 
  chan_to_keep = param$chan_to_keep 
  channels = param$channels
  composite = param$composite
  l_ifd = length(ifd)
  n_ifd = names(ifd)
  
  # set seed if any
  if(param$add_noise) {
    set.seed(param$random_seed)
    on.exit(set.seed(NULL))
  }

  # extract
  foo = lapply(1:l_ifd, FUN=function(i_ifd) {
    img = cpp_extract(fname = fname, ifd = ifd[[i_ifd]], 
                      colors = param$colors, channels = channels, chan_to_extract = param$chan_to_extract - 1, # index start 0 in C, 1 in R, 
                      extract_msk = param$extract_msk, mode = mode, size = param$size, verbose = verbose)
    names(img) = channels[chan_to_extract,"physicalChannel"]
    
    ##### only keep selected channels + channels needed for composite and create composite
    img = c(img[chan_to_keep], lapply(param$composite_desc, FUN=function(i) {
      tmp = img[[as.character(i[1,"int"])]]*i[1,"dec"]
      R = nrow(i)
      if(R>1) for(r in 2:R) tmp = tmp + img[[as.character(i[r,"int"])]]*i[r,"dec"]
      return(tmp)
    }))
    ##### export image
    switch(export,
           "file" = {
             img = lapply(1:length(img), FUN = function(i) {
               export_name = formatn(splitp_obj = splitp_obj,
                                     splitf_obj = splitf_obj,
                                     channel = c(sprintf("Ch%02.f",channels[chan_to_extract,"physicalChannel"]),composite)[i],
                                     object = n_ifd[i_ifd],
                                     bypass = bypass)
               if(file.exists(export_name)) {
                 if(overwrite) {
                   objectWrite(x = img[[i]], type = type, export_name)
                 } else {
                   warning(paste0("file ", export_name, " already exists and will not be overwritten"), call. = FALSE, immediate. = TRUE)
                 }
               } else {
                 objectWrite(x = img[[i]], type = type, export_name)
               }
               return(normalizePath(export_name, winslash = "/", mustWork = FALSE))
             })
           },
           "base64" = {
             if(base64_id) {
               img = lapply(1:length(img), FUN=function(i) {
                 sprintf("<img id=%s %s width='%s' height='%s' src='data:image/%s;base64,%s'>",
                         formatn(splitp_obj = splitp_obj,
                                 splitf_obj = splitf_obj,
                                 channel = c(sprintf("Ch%02.f",channels[chan_to_extract,"physicalChannel"]),composite)[i],
                                 object = n_ifd[i_ifd],
                                 bypass = bypass),
                         base64_att,
                         ncol(img[[i]]),
                         nrow(img[[i]]),
                         type,
                         base64_encode(objectWrite(x = img[[i]], type = type, raw())))
               })
             } else {
               img = lapply(1:length(img), FUN=function(i) {
                 sprintf("<img %s width='%s' height='%s' src='data:image/%s;base64,%s'>",
                         base64_att,
                         ncol(img[[i]]),
                         nrow(img[[i]]),
                         type,
                         base64_encode(objectWrite(x = img[[i]], type = type, raw())))
               })
             }
           })
    
    names(img) <- c(channels$name[channels$physicalChannel %in% chan_to_keep],composite)
    attr(img, "object_id") <- ifd[[i_ifd]]$infos$OBJECT_ID # adds object_id number so as to further check that extracted image is expected one
    attr(img, "offset_id") <- n_ifd[i_ifd] # adds offset_id number further check that extracted mask is expected one
    attr(img, "channel_id") <- c(chan_to_keep, composite) # adds channel_id (physical's one) number so as to be able to create a Gallery
    attr(img, "removal") <- param$removal
    return(img)
  })
  if(export == "file") {
    return(invisible(foo))
  }
  return(foo)
}
