#' @title IFC_object Object Extraction Parameters Definition
#' @description
#' Defines IFC_object object extraction parameters.
#' @param display object of class IFC_display, rich information extracted by \code{\link{getDisplayInfo}}. THis parameter can't be missing.
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
#' @param composite character vector of image composite. Default is "", for no image composite.\cr
#' Should be like "1.05/2.4/4.55" for a composition of 5 perc. of channel 1, 40 perc. of channel 2 and 50 perc. of channel 55.\cr
#' Note that channels should have been acquired and final image composition should be 100 perc., otherwise an error is thrown.\cr
#' Note that each composite will be appended after 'selection'.
#' @param selection physical channels to extract.\cr
#' Note that this parameter will be ordered.\cr
#' Default is "all" to extract all acquired channels.\cr
#' Use "none" to only extract composite.
#' @param random_seed a single value, interpreted as an integer when 'add_noise' is set to TRUE. Default is NULL.
#' @param size a length 2 integer vector of final dimensions of the image, height 1st and width 2nd. Default is c(0,0) for no change.
#' @param force_width whether to use information in display to fill size. Default is TRUE.
#' When set to TRUE, width of 'size' argument will be overwritten.
#' @param removal removal method: Either "none", "clipped", "masked", "MC".\cr
#' -"none", to keep image as is\cr
#' -"clipped", to remove clipped object from image.\cr
#' -"masked", to only keep masked object from image.\cr
#' -"MC", to only keep MC masked object from image.
#' This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to 'selection' and 'composite' parameters.
#' @param add_noise if TRUE adds normal noise to background using rnorm(), from \pkg{Rcpp}. Default is TRUE.\cr
#' Note that it is better to set it to FALSE when 'removal' is "masked" or "MC". Doing so will allow to place masked object in a zero filled background,
#' otherwise background will still be filled with noise.
#' This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to 'selection' and 'composite' parameters.
#' @param full_range only apply when mode is not "raw", if full_range is TRUE, then [0,4095] range will be kept. Default is FALSE.\cr
#' It is like "raw" mode but allowing normalization to [0,1].
#' This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to 'selection' and 'composite' parameters.\cr
#' @param force_range only apply when mode is not "raw", if force_range is TRUE, then image display range will be adjusted for each object resulting in normalization. Default is FALSE.\cr
#' Note that this parameter takes the precedence over 'full_range'.\cr
#' This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to 'selection' and 'composite' parameters.
#' @param ... other arguments to be passed.
#' @return A object of class IFC_param. 
#' @export
objectParam <- function(display, 
                        mode = c("rgb", "gray", "raw")[3],
                        export = c("file", "matrix", "base64")[2],
                        export_to, base64_id = FALSE, base64_att = "",
                        overwrite = FALSE, 
                        composite = "", selection = "all",
                        size = c(0,0), force_width = TRUE, 
                        random_seed = NULL, removal = "none", add_noise = TRUE, 
                        full_range = FALSE, force_range = FALSE,
                        ...) {
  dots=list(...)
  # dots=list(...)
  if(missing(display)) stop("'display' can't be missing") 
  assert(display, cla = "IFC_display")
  
  ##### check size
  size = na.omit(as.integer(size[1:2]))
  assert(size, len=2, typ="integer")
  force_width = as.logical(force_width); assert(force_width, len = 1, alw = c(TRUE,FALSE)) 
  if(force_width) size = c(size[1], as.integer(display$channelwidth))
  
  ##### check input param
  removal = as.character(removal); assert(removal, alw = c("none", "clipped", "masked", "MC"))
  add_noise = as.logical(add_noise); assert(add_noise, alw = c(TRUE,FALSE))
  full_range = as.logical(full_range); assert(full_range, alw = c(TRUE,FALSE))
  force_range = as.logical(force_range); assert(force_range, alw = c(TRUE,FALSE))
  assert(export, len = 1, alw = c("file", "matrix", "base64")) 
  assert(mode, len = 1, alw = c("rgb", "gray", "raw")) 
  if(!missing(random_seed)) {
    random_seed = na.omit(as.integer(random_seed[is.finite(random_seed)]))
    assert(random_seed, len = 1, typ = "integer")
  }
  
  ##### retrieve channels
  channels = display$Images[display$Images$physicalChannel %in% which(display$in_use), ]

  ##### checks selection
  if(any(c("all","none") %in% selection)) {
    if(length(selection)!=1) stop("when 'selection' is \"all\" or \"none\" no other terms are accepted")
    if("all" %in% selection) {
      selection = channels$physicalChannel
    } else {
      selection = NULL
    }
  } else {
    selection = as.integer(selection); mess = assert(selection, alw = channels$physicalChannel)
  }
  
  ##### ensures composite is well formatted
  composite = as.character(composite); assert(composite, typ = "character")
  if(!all(gsub("\\.|\\/|[[:digit:]]","",composite) %in% "")) stop("'composite' is not well formatted")
  
  ##### treats composite
  if(any(composite != "")) {
    composite = setdiff(composite, "")
    composite_desc = lapply(strsplit(composite, split = "/", fixed = TRUE), FUN=function(comp) {
      tmp = do.call(rbind, args = strsplit(comp, split = ".", fixed = TRUE))
      return(cbind("int"=as.numeric(tmp[,1]), "dec"=as.numeric(tmp[,2])/10^nchar(tmp[,2])))
    })
    composite_chan = unique(unlist(lapply(composite_desc, FUN=function(x) {
      if(sum(x[, "dec"]) != 1) stop("'composite' final composition is not 100%")
      return(x[, "int"])
    })))
    tmp = !(composite_chan %in% channels$physicalChannel)
    if(any(tmp)) stop(paste0("'composite' requires channel",ifelse(sum(tmp)>1,"s","")," [",
                             paste0(composite_chan[tmp],collapse=","),"] which ",ifelse(sum(tmp)>1,"are","is"),
                             " not part of physical channels acquired [",paste0(channels$physicalChannel, collapse=","),"]"))
    composite = gsub("\\/", ",", composite)
    } else {
    composite = NULL
    composite_chan = NULL
    composite_desc = list()
  }
  
  ##### define chan to extract, to keep, ...
  chan_to_extract = which(channels$physicalChannel %in% c(selection, composite_chan))
  chan_to_keep = as.character(channels$physicalChannel[channels$physicalChannel %in% selection])
  if(sum(chan_to_extract) == 0) stop("can't export object with 'selection'/'composite' parameters provided")
  
  ##### fill removal, add_noise, full_range, force_range for every extracted channels
  channels[,"removal"] <- "none"
  channels[chan_to_extract,"removal"] <- rep_len(removal, length.out = length(chan_to_extract))
  channels[,"removal"] = as.integer(factor(x = channels[,"removal"], levels = c("none", "clipped", "masked", "MC"))) - 1
  channels[,"add_noise"] <- F
  channels[chan_to_extract,"add_noise"] <- rep_len(add_noise, length.out = length(chan_to_extract))
  channels[,"full_range"] <- F
  channels[chan_to_extract,"full_range"] <- rep_len(full_range, length.out = length(chan_to_extract))
  channels[,"force_range"] <- F
  channels[chan_to_extract,"force_range"] <- rep_len(force_range, length.out = length(chan_to_extract))
  
  ##### checks gamma
  gamma = channels[,c("gamma")]
  # recomputes gamma to check if c("xmin", "xmax", "xmid", "ymid") were not modified
  gamma_c = apply(channels[,c("xmin", "xmax", "xmid", "ymid")], 1, cpp_computeGamma)
  # if modified gamma is reset to 1 to use linear visualization
  gamma[gamma != gamma_c] <- 1
  channels[,"gamma"] <- gamma
  
  # minimal ans for export == "matrix"
  ans = list(mode = mode,
             export = export,
             base64_id = FALSE,
             base64_att = "",
             colors = sapply(channels[,"color"], simplify = FALSE, FUN=function(x) {tmp = c(rgb2hsv(col2rgb(x))); names(tmp) = x; tmp }), 
             channels = channels, 
             chan_to_extract = chan_to_extract,
             chan_to_keep = chan_to_keep,
             removal = channels[channels$physicalChannel[channels$physicalChannel %in% selection], "removal"],
             add_noise = any(channels$add_noise),
             composite = composite,
             composite_desc = composite_desc,
             extract_msk = max(channels[, "removal"]), 
             size = size,
             random_seed = random_seed,
             checksum = display$checksum,
             fileName_image = display$fileName_image)
  
  # compute extre param for export == "file" or ""base64"
  if(export != "matrix") { # file or base64
    if(mode == "raw") stop("can't export as \"raw\" when '", export, "' is choosen")
    if(export == "file") { # file
      assert(overwrite, len = 1, alw = c(TRUE, FALSE))
      if(missing(export_to)) stop("'export_to' can't be missing when 'export' is \"file\"")
      # determines type from export_to
      type = getFileExt(export_to)
      switch(type,
             "jpg" = {type <- "jpeg"},
             "tif" = {type <- "tiff"}
      )
      # check type
      assert(type, len = 1, alw = c("bmp", "jpg", "png", "tif"))
      splitf_obj = splitf(display$fileName_image)
      splitp_obj = splitp(export_to)
      # create folder for file export
      dir_name = dirname(formatn(splitp_obj, splitf_obj))
      ans = c(list(dir_name = dir_name,
                   export_to = export_to,
                   type = type,
                   overwrite = overwrite,
                   splitf_obj = splitf_obj,
                   splitp_obj = splitp_obj), ans)
    } else { # base64
      base64_id = as.logical(base64_id)
      assert(base64_id, len = 1, alw = c(TRUE,FALSE))
      base64_att = na.omit(as.character(base64_att))
      assert(base64_att, len = 1, typ = "character")
      ans$base64_id = base64_id
      ans$base64_att = base64_att
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
      # check type
      assert(type, len = 1, alw = c("bmp", "jpg", "png", "tif"))
      ans = c(list(export_to = export_to,
                   type = type), ans)
      if(base64_id) {
        ans = c(list(splitf_obj = splitf(display$fileName_image),
                     splitp_obj = splitp(export_to)), ans)
      }
    }
  }
  class(ans) <- "IFC_param"
  return(ans)
}
