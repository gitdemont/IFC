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

#' @title Object Extraction Parameters Definition
#' @description
#' Defines `IFC_object` object extraction parameters.
#' @param ... arguments to be passed to \code{\link{getInfo}}, \strong{only} if \code{'info'} is \strong{not} provided.
#' @param info object of class `IFC_info`, rich information extracted by \code{\link{getInfo}}. 
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.
#' If missing, the default, \code{'info'} will be extracted thanks to \code{'...'}.
#' @param mode color mode export. Either \code{"rgb"}, \code{"gray"} or \code{"raw"}. Default is \code{"raw"}.
#' Note that \code{"raw"} is only possible when \code{'export'} is \code{"matrix"} or \code{"multi"}.
#' @param export format mode export. Either \code{"file"}, \code{"matrix"}, \code{"base64"}, or \code{"multi"}. Default is \code{"matrix"}.
#' @param write_to used when export is not \code{"matrix"} to compute exported file name or base64 id attribute.\cr
#' Exported file extension and base64 MIME type will be deduced from this pattern. Allowed export are \code{".bmp"}, \code{".jpg"}, \code{".jpeg"}, \code{".png"}, \code{".tif"}, \code{".tiff"}.
#' Note that \code{".bmp"} are faster but not compressed producing bigger data.\cr
#' Placeholders, if found, will be substituted:\cr
#' -\code{\%d}: with full path directory,\cr
#' -\code{\%p}: with first parent directory,\cr
#' -\code{\%e}: with extension (without leading .),\cr
#' -\code{\%s}: with shortname (i.e. basename without extension),\cr
#' -\code{\%o}: with object_id,\cr
#' -\code{\%c}: with channel_id (not possible when \code{'export'} is \code{"multi"}).\cr
#' A good trick is to use:\cr
#' -\code{"\%d/\%s/\%s_\%o_\%c.tiff"}, when \code{'export'} is \code{"file"},\cr
#' -\code{"\%d/\%s/\%s_\%o.tiff"}, when \code{'export'} is \code{"multi"},\cr
#' -\code{"\%o_\%c.bmp"}, when \code{'export'} is \code{"base64"}.\cr
#' Note that if missing and \code{'export'} is not \code{"file"}, \code{'write_to'} will be set to \code{"\%o_\%c.bmp"}.
#' @param base64_id whether to add id attribute to base64 exported object. Default is \code{FALSE}.\cr
#' Only applied when export is \code{"base64"}.
#' @param base64_att attributes to add to base64 exported object. Default is \code{""}.\cr
#' Only applied when export is \code{"base64"}. For example, use \code{"class='draggable'"}.\cr
#' Note that \code{id} (if \code{'base64_id'} is \code{TRUE}) and \code{width} and \code{height} are already used.
#' @param overwrite only apply when \code{'export'} is \code{"file"} whether to overwrite file or not. Default is \code{FALSE}.
#' @param composite character vector of image composite. Default is \code{""}, for no image composite.\cr
#' Should be like \code{"1.05/2.4/4.55"} for a composition of 5 perc. of channel 1, 40 perc. of channel 2 and 50 perc. of channel 55.\cr
#' Note that channels should have been acquired and final image composition should be 100 perc., otherwise an error is thrown.\cr
#' Note that each composite will be appended after \code{'selection'}.\cr
#' Note that composite will be forced to \code{""} when \code{'export'} is \code{"multi"}.
#' @param selection physical channels to extract.\cr
#' Note that this parameter will be ordered.\cr
#' Default is \code{"all"} to extract all acquired channels.\cr
#' Use \code{"none"} to only extract composite.
#' @param random_seed a list of elements to pass to \link[base]{set.seed} or a single value, interpreted as an integer, or NULL to be used when \code{'add_noise'} is set to \code{TRUE}. Default is \code{NULL}.
#' Note that \code{NA_integer_} or \code{list(seed = NA_integer_)} can be used to not call \link[base]{set.seed} at all.
#' @param size a length 2 integer vector of final dimensions of the image, height 1st and width 2nd. Default is \code{c(0,0)} for no change.
#' @param force_width whether to use information in \code{'info'} to fill size. Default is \code{TRUE}.
#' When set to \code{TRUE}, width of \code{'size'} argument will be overwritten.
#' @param removal removal method: Either \code{"none"}, \code{"raw"}, \code{"clipped"}, \code{"masked"}, \code{"MC"}.\cr
#' -\code{"none"}, to keep image as is, no mask extraction will be performed resulting in faster extraction,\cr
#' -\code{"raw"}, to keep image as is, it provides a convenient way to retrieve \code{"raw"} value for the mask,\cr
#' -\code{"clipped"}, to remove clipped object from image,\cr
#' -\code{"masked"}, to only keep masked object from image,\cr
#' -\code{"MC"}, to keep MC masked object from image.
#' This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to \code{'selection'} and \code{'composite'} parameters.
#' @param add_noise if \code{TRUE} adds normal noise to background using \pkg{Rcpp}. Default is \code{TRUE}.\cr
#' Note that it is better to set it to \code{FALSE} when \code{'removal'} is \code{"masked"} or \code{"MC"}. Doing so will allow to place masked object in a zero filled background,
#' otherwise background will still be filled with noise.
#' This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to \code{'selection'} and \code{'composite'} parameters.
#' @param full_range only apply when \code{'mode'} is not \code{"raw"}, if \code{'full_range'} is \code{TRUE}, then object range will be considered as 0 to 4095, it is like \code{"raw"} \code{'mode'} but resulting in [0,4095] normalization to [0,1]. Default is \code{FALSE}.\cr
#' This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to \code{'selection'} and \code{'composite'} parameters.
#' @param force_range only apply when \code{'mode'} is not \code{"raw"}, if \code{'force_range'} is \code{TRUE}, then range will be adjusted to object range in \code{[-4095,+inf]} resulting in normalization to [0,1]. Default is \code{FALSE}.\cr
#' This parameter will be repeated with rep_len() from \pkg{base} for every physical channel that needs to be extracted according to \code{'selection'} and \code{'composite'} parameters.\cr
#' Note that this parameter takes the precedence over \code{'full_range'}.
#' @param spatial_correction only apply on RIF file, whether to apply spatial correction. Default is \code{FALSE}.
#' @details When a mask is detected, \code{'add_noise'}, \code{'full_range'} and \code{'force_range'} are set to \code{FALSE} and range used will be forced to \code{[0,3]}.\cr\cr
#' Range of image is controlled by \code{'Images'} information from supplied \code{'info'} or as extracted by \code{\link{getInfo}} and will be returned as \code{'channels'} by \code{\link{objectParam}}.
#' In case \code{'mode'} is not \code{"raw"}, '\code{channels$xmin}', '\code{channels$xmax}', '\code{channels$gamma}' will be used for object extraction by \code{\link{objectExtract}} unless any of \code{'force_range'} or \code{'full_range'} is \code{TRUE}.\cr\cr
#' Experimental (as of v0.2.0.501): once returned by \code{\link{objectParam}}, those '\code{channels$xmin}' and '\code{channels$xmax}' can be manually adjusted to \code{]0,1[} so as to be used as \code{'probs'} argument to \link[stats]{quantile} to allow quantile normalization during object extraction (\code{\link{objectExtract}}) afterwards.
#' @return an object of class `IFC_param`. 
#' @export
objectParam <- function(...,
                        info,
                        mode = c("rgb", "gray", "raw")[3],
                        export = c("file", "matrix", "base64", "multi")[2],
                        write_to,
                        base64_id = FALSE,
                        base64_att = "",
                        overwrite = FALSE, 
                        composite = "",
                        selection = "all",
                        size = c(0,0),
                        force_width = TRUE, 
                        random_seed = NULL,
                        removal = "none",
                        add_noise = TRUE, 
                        full_range = FALSE,
                        force_range = FALSE,
                        spatial_correction = FALSE) {
  dots=list(...)
  
  #### check input
  input = whoami(entries = as.list(match.call()), search = list(info = "IFC_info"))
  info = input$info
  if(length(info) == 0) { # info was not found. use extra param to get info
    param_info = names(dots) %in% c("fileName","from","verbose",
                                    "verbosity","warn","force_default",
                                    "cifdir","ntry")
    info = do.call(what = "getInfo", args = dots[param_info])  
  }
  
  # TODO add the folowing lines
  # provided = names(as.list(match.call())[-(unique(attr(input, "was"))+1)])
  # expected = setdiff(names(formals(objectParam)), c("...", "info",
  #                                                   "fileName","from","verbose",
  #                                                   "verbosity","warn","force_default",
  #                                                   "cifdir","ntry"))
  # matches = charmatch(provided, expected)
  # untreated = provided[is.na(matches)]
  # multiple = provided[!is.na(matches) & matches == 0]
  # if(length(untreated) != 0) warning(paste0("objectParam: provided argument", ifelse(length(untreated) == 1, "", "s"), " [", paste0("'", untreated, "'", collapse = ","), "] will not be used"))
  # if(length(multiple) != 0) warning(paste0("objectParam: provided argument", ifelse(length(multiple) == 1, "", "s"), " [", paste0("'", multiple, "'", collapse = ","), "] match with several parameters"))
  #
  
  ##### check size
  size = na.omit(as.integer(size[1:2]))
  assert(size, len=2, typ="integer")
  force_width = as.logical(force_width); assert(force_width, len = 1, alw = c(TRUE,FALSE)) 
  if(force_width) size = c(size[1], as.integer(info$channelwidth))
  
  ##### check input param
  removal = as.character(removal); assert(removal, alw = c("none", "raw", "clipped", "masked", "MC"))
  add_noise = as.logical(add_noise); assert(add_noise, alw = c(TRUE,FALSE))
  full_range = as.logical(full_range); assert(full_range, alw = c(TRUE,FALSE))
  force_range = as.logical(force_range); assert(force_range, alw = c(TRUE,FALSE))
  spatial_correction = as.logical(spatial_correction); assert(spatial_correction, alw = c(TRUE, FALSE))
  assert(export, len = 1, alw = c("file", "matrix", "base64", "multi")) 
  assert(mode, len = 1, alw = c("rgb", "gray", "raw")) 
  overwrite = as.logical(overwrite); assert(overwrite, alw = c(TRUE, FALSE))
  
  ##### retrieve channels
  channels = info$Images[info$Images$physicalChannel %in% which(info$in_use), ]
  
  ##### check selection
  if(any(c("all","none") %in% selection)) {
    if(length(selection)!=1) stop("when 'selection' is \"all\" or \"none\" no other terms are accepted")
    if("all" %in% selection) {
      sel_int = channels$physicalChannel
    } else {
      sel_int = NULL
    }
  } else {
    sel_int = as.integer(selection); mess = assert(sel_int, alw = channels$physicalChannel)
  }
  
  ##### ensure composite is well formatted
  if((export == "multi") && !identical(composite, "")) {
    warning("'composite' has been forced to \"\" for 'export'=\"multi\"")
    composite = ""
  }
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
  chan_to_extract = which(channels$physicalChannel %in% c(sel_int, composite_chan))
  chan_to_keep = as.character(channels$physicalChannel[channels$physicalChannel %in% sel_int])
  if(sum(chan_to_extract) == 0) stop("can't export object with 'selection'/'composite' parameters provided")
  
  ##### fill removal, add_noise, full_range, force_range for every extracted channels
  channels[,"string_removal"] <- "none"
  XIF_test = ifelse(length(info$XIF_test) == 0, testXIF(info$fileName_image), info$XIF_test) 
  
  if(XIF_test > 0) {
    channels[chan_to_extract,"string_removal"] <- rep_len(removal, length.out = length(chan_to_extract))
  } else {
    if(any(removal != "none")) warning("'removal' was forced to \"none\" because no mask can be found in the file", call. = FALSE, immediate. = TRUE)
  }
  channels[,"removal"] = as.integer(factor(x = channels[,"string_removal"], levels = c("none", "raw", "clipped", "masked", "MC"))) - 1
  channels[,"add_noise"] <- FALSE
  channels[chan_to_extract,"add_noise"] <- rep_len(add_noise, length.out = length(chan_to_extract))
  channels[,"full_range"] <- FALSE
  channels[chan_to_extract,"full_range"] <- rep_len(full_range, length.out = length(chan_to_extract))
  channels[,"force_range"] <- FALSE
  channels[chan_to_extract,"force_range"] <- rep_len(force_range, length.out = length(chan_to_extract))
  
  ##### check gamma
  gamma = channels[,c("gamma")]
  # recompute gamma to check if c("xmin", "xmax", "xmid", "ymid") were not modified
  gamma_c = apply(channels[,c("xmin", "xmax", "xmid", "ymid")], 1, cpp_computeGamma)
  # if modified gamma is reset to 1 to use linear visualization
  gamma[gamma != gamma_c] <- 1
  channels[,"gamma"] <- gamma
  
  # add spatial correction information
  if(spatial_correction) {
    if(getFileExt(info$fileName_image) == "rif") {
      mag = switch(as.character(info$magnification), "20" = "20x", "40" = "", "60" = "60x")
      ASSISTDb = getASSIST(info$fileName_image)
      spa_off = sapply(c("X","Y"), FUN = function(off) {
        as.numeric(strsplit(ASSISTDb[[paste0(off, "Offsets", mag, "_Gen2_0_11")]], split = " ", fixed = TRUE)[[1]])
      })
      channels[, "spatial_X"] <- spa_off[, "X"][chan_to_extract]
      channels[, "spatial_Y"] <- spa_off[, "Y"][chan_to_extract]
    } else {
      warning("'spatial_correction' can only be applied on .rif file")
    }
  }
  
  ##### build ans
  ans = list(mode = mode,
             export = export,
             write_to = NULL,
             type = "",
             splitf_obj = NULL,
             splitp_obj = NULL,
             base64_id = FALSE,
             base64_att = "",
             overwrite = overwrite,
             colors = sapply(channels[,"color"], simplify = FALSE, FUN=function(x) {tmp = c(rgb2hsv(col2rgb(x))); names(tmp) = x; tmp }), 
             channels = channels, 
             chan_to_extract = chan_to_extract,
             chan_to_keep = chan_to_keep,
             removal = channels[chan_to_extract, "string_removal"],
             add_noise = any(channels$add_noise),
             composite = composite,
             selection = selection,
             composite_desc = composite_desc,
             extract_msk = max(channels[, "removal"]), 
             size = size,
             random_seed = fetch_seed(random_seed),
             objcount = info$objcount,
             channelwidth = info$channelwidth,
             in_use = info$in_use,
             brightfield = info$brightfield,
             coremode = info$coremode,
             magnification = info$magnification,
             checksum = info$checksum,
             XIF_test = XIF_test,
             fileName_image = info$fileName_image)
  
  ##### compute extre param for export == "file" or ""base64"
  ans$splitf_obj <- splitf(ans$fileName_image)
  if(export != "matrix") { 
    if((mode == "raw" && export != "multi") ||
       (mode != "raw" && export == "multi")) stop("can't 'export' to \"",export,"\" when 'mode' \"", mode, "\" is choosen")
    if(export == "file" || export == "multi") { # file
      # not allowed to write file without user input
      if(missing(write_to)) stop("'write_to' can't be missing when 'export' is \"file\" or \"multi\"")
    } else { # base64
      base64_id = as.logical(base64_id)
      assert(base64_id, len = 1, alw = c(TRUE,FALSE))
      ans$base64_id = base64_id
      base64_att = na.omit(as.character(base64_att))
      assert(base64_att, len = 1, typ = "character")
      ans$base64_att = base64_att
    }
  }
  if(missing(write_to)) {
    ans$write_to = "%o_%c.bmp"
  } else {
    ans$write_to <- write_to
  }
  ans$splitp_obj = splitp(ans$write_to)
  ans$dir_name <- dirname(formatn(splitp_obj = splitp(ans$write_to), splitf_obj = ans$splitf_obj))
  type = getFileExt(ans$write_to)
  switch(type,
         "jpg" = {type <- "jpeg"},
         "tif" = {type <- "tiff"} )
  ##### check type
  assert(type, len = 1, alw = c("bmp", "jpeg", "png", "tiff"))
  if(export == "multi" && type != "tiff") stop("when 'export' is \"multi\", file extension has to be \"tiff\" not \"", type, "\"")
  ans$type <- type
  
  class(ans) <- "IFC_param"
  return(ans)
}
