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

#' @title Object Extraction
#' @description
#' Extracts / Decompress objects stored in RIF or CIF Files.
#' @param ifd list of sub elements of IFD data information extracted by \code{\link{getIFD}}. This argument can't be missing.
#' @param param object of class `IFC_param`, containing extraction parameters defined by \code{\link{objectParam}}.\cr
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.
#' If this argument is missing, \code{\link{objectExtract}} will use named arguments in \code{'...'} as input to \code{\link{objectParam}} to control object extraction.\cr
#' However, if \code{'param'} is provided, \code{'...'} will be \strong{ignored}.
#' @param verbose whether to display information (use for debugging purpose). Default is \code{FALSE}.
#' @param bypass whether to bypass checks on \code{'ifd'} and \code{'param'}. Default is \code{FALSE}.
#' @inheritDotParams objectParam
#' @note \code{'...'} will be used to compute \code{'param'} \strong{only} if \code{'param'} is \strong{not} provided.\cr
#' \strong{/!\\} When \code{'param'} is not provided, if not any of \code{'fileName'}, \code{'info'} can be found in \code{'...'}, then \code{attr(ifd, "fileName_image")} will be used as \code{'fileName'} input parameter to pass to \code{\link{objectParam}}.
#' @inherit objectParam details
#' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
#' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}\cr
#' \verb{
#' BSD implementations of Bio-Formats readers and writers
#' %%
#' Copyright (C) 2005 - 2017 Open Microscopy Environment:
#'   - Board of Regents of the University of Wisconsin-Madison
#'   - Glencoe Software, Inc.
#'   - University of Dundee
#' %%
#' Redistribution and use in source and binary forms, with or without
#' modification, are permitted provided that the following conditions are met:
#' 
#' 1. Redistributions of source code must retain the above copyright notice,
#'    this list of conditions and the following disclaimer.
#' 2. Redistributions in binary form must reproduce the above copyright notice,
#'    this list of conditions and the following disclaimer in the documentation
#'    and/or other materials provided with the distribution.
#'
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
#' }
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a cif file
#'   file_cif <- system.file("extdata", "example.cif", package = "IFCdata")
#'   cif_offs <- getOffsets(fileName = file_cif, fast = TRUE)
#'   ## extract infomation
#'   info <- getInfo(fileName = file_cif, from = "analysis")
#'   ## retrieve number of objects stored
#'   nobj <- as.integer(info$objcount)
#'   ## randomly subset the offsets of at most 5 "img" objects
#'   sel = sample(0:(nobj-1), min(5, nobj))
#'   sub_offs <- subsetOffsets(cif_offs, objects = sel, image_type = "img")
#'   ## read IFDs from these "img" objects
#'   IFDs <- getIFD(fileName = file_cif, offsets = sub_offs)
#'   ## extract raw data of these"img" objects to matrix
#'   raw = objectExtract(ifd = IFDs, info = info, mode = "raw", 
#'                       export = "matrix")
#'   ## extract base64 "rgb" colorized version of these "img" objects to base64
#'   b64 = objectExtract(ifd = IFDs, info = info, mode = "rgb", 
#'                       export = "base64", base64_id = TRUE,
#'                       write_to = "example_%o_%c.bmp")
#'   ## use DisplayGallery to show the first "img" objects and play with ... extra parameters
#'   ## force_range, add_noise, selection, composite, see objectParam
#'   DisplayGallery(info = info, offsets = cif_offs, objects = sel,
#'                  base64_id = TRUE, write_to = "example_%o_%c.bmp",
#'                  force_range = c(FALSE,TRUE,FALSE,TRUE), add_noise = FALSE,
#'                  selection = c(1,2,4,6), composite = "1.7/4.3")
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return A list (for every extracted objects) of list (for every exported channels) depending on '\code{param$export}' parameter:\cr
#' -\code{"matrix"}, a matrix when '\code{param$mode}' is set to \code{"raw"} or \code{"gray"} OR an array when '\code{param$mode}' is \code{"rgb"},\cr
#' -\code{"base64"}, a data-uri string,\cr
#' -\code{"file"}, an invisible file path corresponding to the location of exported file(s).\cr
#' -\code{"multi"}, an invisible file path corresponding to the location of exported file(s).
#' @export
objectExtract <- function(ifd, 
                          param, 
                          verbose = FALSE, 
                          bypass = FALSE, 
                          ...) {
  dots=list(...)
  assert(verbose, len = 1, alw = c(TRUE, FALSE))
  assert(bypass, len = 1, alw = c(TRUE, FALSE))
  # bypass ifd / param checking
  if(!bypass) {
    #### check input
    input = whoami(entries = as.list(match.call()), search = list(ifd = "IFC_ifd_list",
                                                                  param = "IFC_param",
                                                                  info = "IFC_info"))
    ifd = input$ifd
    param = input$param
    info = input$info
    fileName = input$fileName
    if(length(ifd) == 0) stop("'ifd' can't be missing")
    if("IFC_first_ifd" %in% class(ifd)) stop("can't extract object from 'ifd' of class `IFC_first_ifd`")
    if(length(param) == 0) {
      dots = dots[!(names(dots) %in% c("fileName", "verbose"))]
      if((length(fileName) == 0) && (length(info) == 0)) {
        param = do.call(what = "objectParam", args = c(list(fileName = attr(ifd, "fileName_image"), verbose = verbose), dots))
      } else {
        dots = dots[!(names(dots) %in% c("info"))]
        param = do.call(what = "objectParam", args = c(list(verbose = verbose, info = quote(info), fileName = fileName), dots))
      }
    } else {
      if(attr(ifd, "checksum") != param$checksum) stop("'ifd' and 'param' do not match, please ensure that they originate from same file")
    }
  }
  # create dir to export files
  if(param$export == "file" || param$export == "multi") if(!dir.exists(param$dir_name)) if(!dir.create(param$dir_name, recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", param$dir_name))
  
  # shortcut
  chan_to_extract = param$chan_to_extract
  chan_to_keep = param$chan_to_keep 
  channels = param$channels
  composite = param$composite
  
  # retrieve length + names of IFDs
  l_ifd = length(ifd)
  n_ifd = names(ifd)
  
  # set seed if any
  f = function(x) { x }
  if(param$add_noise) {
    SEED = param$random_seed
    f = function(x) {
      with_seed(x, SEED$seed, SEED$kind, SEED$normal.kind, SEED$sample.kind)
    }
  }
  if(sum(names(param$channels) %in% c("spatial_X","spatial_Y")) == 2) {
    spatialX = param$channels$spatial_X
    spatialY = param$channels$spatial_Y
  } else {
    spatialX = rep(NA_real_, nrow(param$channels))
    spatialY = spatialX
  }

  # extracts
  foo = f(lapply(1:l_ifd, FUN=function(i_ifd) {
    img = cpp_extract(fname = enc2native(param$fileName_image), 
                      ifd = ifd[[i_ifd]], 
                      colors = param$colors,  
                      physicalChannel = channels$physicalChannel,
                      xmin = channels$xmin,
                      xmax = channels$xmax,
                      spatialX = spatialX,
                      spatialY = spatialY,
                      removal = channels$removal,
                      add_noise = channels$add_noise,
                      full_range = channels$full_range,
                      force_range = channels$force_range,
                      gamma = channels$gamma,
                      chan_to_extract = param$chan_to_extract - 1, # index start 0 in C, 1 in R, 
                      extract_msk = param$extract_msk, 
                      mode = param$mode, 
                      size = param$size, 
                      verbose = verbose)
    names(img) = channels[chan_to_extract,"physicalChannel"]
    
    ##### only keep selected channels + channels needed for composite and create composite
    img = c(img[chan_to_keep], lapply(param$composite_desc, FUN=function(i) {
      tmp = img[[as.character(i[1,"int"])]]*i[1,"dec"]
      R = nrow(i)
      if(R>1) for(r in 2:R) tmp = tmp + img[[as.character(i[r,"int"])]]*i[r,"dec"]
      return(tmp)
    }))
    
    ##### export image
    switch(param$export,
           "multi" = { #TODO add code to fill ImageDescription TAG 270
             meta = list("269" = list(tag = 269, typ = 2, map = n_ifd[i_ifd]),
                         "270" = list(tag = 270, typ = 2, map = "N/A"),
                         "306" = list(tag = 306, typ = 2, map = formatdate(suppressWarnings(getFullTag(ifd, i_ifd, "306")))),
                         "33003" = list(tag = 33003, typ = 4, map = ifd[[i_ifd]]$infos$OBJECT_ID))
             export_name = formatn(splitp_obj = param$splitp_obj,
                                   splitf_obj = param$splitf_obj,
                                   object = n_ifd[i_ifd])
             if(file.exists(export_name)) {
               if(param$overwrite) {
                 objectWrite(x = array(unlist(img, use.names = FALSE, recursive = TRUE), dim = c(dim(img[[1]]), length(img), 1)),
                             type = "multi", export_name, tags = meta)
               } else {
                 warning(paste0("multi ", export_name, " already exists and will not be overwritten"), call. = FALSE, immediate. = TRUE)
               }
             } else {
               if(!dir.exists(dirname(export_name))) if(!dir.create(dirname(export_name), recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dirname(export_name)))
               objectWrite(x = array(unlist(img, use.names = FALSE, recursive = TRUE), dim = c(dim(img[[1]]), length(img), 1)),
                           type = "multi", export_name, tags = meta)
             }
             img = normalizePath(export_name, winslash = "/", mustWork = FALSE)
           }, 
           "file" = {
             meta = list("269" = list(tag = 269, typ = 2, map = n_ifd[i_ifd]),
                         "270" = list(tag = 270, typ = 2, map = "N/A"),
                         "306" = list(tag = 306, typ = 2, map = formatdate(suppressWarnings(getFullTag(ifd, i_ifd, "306")))),
                         "33003" = list(tag = 33003, typ = 4, map = ifd[[i_ifd]]$infos$OBJECT_ID))
             img = lapply(1:length(img), FUN = function(i) {
               export_name = formatn(splitp_obj = param$splitp_obj,
                                     splitf_obj = param$splitf_obj,
                                     channel = c(sprintf("Ch%02.f",channels[chan_to_extract,"physicalChannel"]),composite)[i],
                                     object = n_ifd[i_ifd])
               if(file.exists(export_name)) {
                 if(param$overwrite) {
                   objectWrite(x = img[[i]], type = param$type, export_name, tags = meta)
                 } else {
                   warning(paste0("file ", export_name, " already exists and will not be overwritten"), call. = FALSE, immediate. = TRUE)
                 }
               } else {
                 if(!dir.exists(dirname(export_name))) if(!dir.create(dirname(export_name), recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dirname(export_name)))
                 objectWrite(x = img[[i]], type = param$type, export_name, tags = meta)
               }
               return(normalizePath(export_name, winslash = "/", mustWork = FALSE))
             })
           },
           "base64" = {
             meta = list("269" = list(tag = 269, typ = 2, map = n_ifd[i_ifd]),
                         "270" = list(tag = 270, typ = 2, map = "N/A"),
                         "306" = list(tag = 306, typ = 2, map = formatdate(suppressWarnings(getFullTag(ifd, i_ifd, "306")))),
                         "33003" = list(tag = 33003, typ = 4, map = ifd[[i_ifd]]$infos$OBJECT_ID))
             if(param$base64_id) {
               img = lapply(1:length(img), FUN=function(i) {
                 sprintf("<img id='%s' %s width='%i' height='%i' src='data:image/%s;base64,%s'>",
                         formatn(splitp_obj = param$splitp_obj,
                                 splitf_obj = param$splitf_obj,
                                 channel = c(sprintf("Ch%02.f",channels[chan_to_extract,"physicalChannel"]),composite)[i],
                                 object = n_ifd[i_ifd]),
                         param$base64_att,
                         ncol(img[[i]]),
                         nrow(img[[i]]),
                         param$type,
                         cpp_base64_encode(objectWrite(x = img[[i]], type = param$type, raw())))
               })
             } else {
               img = lapply(1:length(img), FUN=function(i) {
                 sprintf("<img %s width='%i' height='%i' src='data:image/%s;base64,%s'>",
                         param$base64_att,
                         ncol(img[[i]]),
                         nrow(img[[i]]),
                         param$type,
                         cpp_base64_encode(objectWrite(x = img[[i]], type = param$type, raw())))
               })
             }
           })
    if(param$export != "multi") names(img) <- c(channels$name[channels$physicalChannel %in% chan_to_keep],composite)
    attr(img, "object_id") <- ifd[[i_ifd]]$infos$OBJECT_ID # adds object_id number so as to further check that extracted image is expected one
    attr(img, "offset_id") <- n_ifd[i_ifd] # adds offset_id number further check that extracted mask is expected one
    attr(img, "channel_id") <- c(chan_to_keep, composite) # adds channel_id (physical's one) number so as to be able to create a Gallery
    attr(img, "removal") <- param$removal
    return(img)
  }))
  if(param$export == "file" || param$export == "multi") {
    return(invisible(foo))
  }
  return(foo)
}
