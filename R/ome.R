################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2024 Yohann Demont                                             #
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

################################################################################
#             functions described hereunder are experimental                   #
#              inputs and outputs may change in the future                     #
################################################################################

#' @title Image Field Directory Tag Editor
#' @description
#' Edits tag value in file IFDs (Image Field Directory).
#' @param IFD an object of class `IFC_ifd_list` extracted by \code{\link{getIFD}}.
#' @param which scalar, integer (index) or the name of 'IFD' sub-element to edit 'tag' to.
#' @param tag scalar, integer (index) or the name of the desired 'tag' in IFD[[which]].
#' @param value value used to fill the 'tag'. For the moment only string is supported.
#' @param force whether to force edition. Needed when 'tag' value already exceed 4 bytes. Default is \code{FALSE}.
#' @details /!\ file will be modified /!\
#' /!\ When \code{'force'} needs to be \code{TRUE}, 'value' will be appended to file \strong{but} former tag 'value' will still be present in file \strong{and not} erased.\cr
#' Will only work when 'tag' is present in TIFF file.\cr
#' Will only work for TIFF files produced by \pkg{IFC}
#' @return It invisibly returns path of edited file.
#' @keywords internal
editTag <- function(IFD, which, tag, value, force = FALSE) {
  if(!("IFC_ifd_list" %in% class(IFD))) stop("'IFD' object is not of class `IFC_ifd_list`")
  fileName = attr(x = IFD, which = "fileName_image")
  endian = ifelse(readBin(fileName, what = "raw", n = 1) == 0x49, "little", "big")
  if(!missing(value)) {
    if(length(value) == 0) {
      warning("empty 'value' is not supported")
      return(invisible(fileName)) 
    }
    stopifnot(length(value) == 1, typeof(value) == "character")
  } else {
    return(invisible(fileName))
  }
  force = na.omit(force)
  stopifnot(length(which) == 1, length(tag) == 1, length(force) == 1, typeof(force) == "logical")
  i_ifd = 0
  if(typeof(which) == "character") {
    tmp = names(IFD) %in% which
    if(any(tmp)) i_ifd = head(na.omit(which(tmp)))
  } else {
    which = na.omit(as.integer(which))
    if(which >= 1 && which <= length(IFD)) i_ifd = which
  }
  if(i_ifd == 0) {
    warning("can't find 'which'[", which, "] in 'IFD'")
    return(invisible(fileName))
  }
  if(!identical(substr(suppressWarnings(IFC::getFullTag(IFD, which = which, tag = "305")),0,3), "IFC")) stop("'fileName' should have been produced by `IFC` package")
  ifd = IFD[[which]]
  i_tag = 0
  if(typeof(tag) == "character") {
    tmp = names(ifd$tags) == tag
    if(any(tmp)) i_tag = head(na.omit(which(tmp)))
  } else {
    tag = na.omit(as.integer(tag))
    if((tag >= 1) && (tag <= length(ifd$tags))) i_tag = tag
  }
  if(i_tag == 0) {
    warning("can't find 'tag'[", tag, "] in IFD[[", which, "]] tags")
    return(invisible(fileName))
  }
  if(ifd$tags[[i_tag]]$typ != 2) {
    warning("'tag'[", tag, "] in IFD[[", which, "]] of type=",ifd[[i_tag]]$typ," is not supported (should be 2)")
    return(invisible(fileName))
  }
  if(ifd$tags[[i_tag]]$byt > 4 && !force) {
    warning("'tag'[", tag, "] in IFD[[", which, "]] can't be modified unless 'force'=TRUE")
    return(invisible(fileName))
  }
  tag_offset = ifd$curr_IFD_offset + 2 + (i_tag - 1) * 12
  raw = unlist(c(iconv(value, from = "UTF8", to = "UTF8", toRaw = TRUE), raw(1)), use.names = FALSE, recursive = FALSE)
  l = length(raw)
  fs = file.size(fileName)
  if(l + fs > (2^32-1)) stop("'value' is too big to be stored")
  if(l > 2^31-1) l = l - 2^32
  towrite = file(description = enc2native(fileName), open = "r+b")
  on.exit(close(towrite))
  # modify tag count
  seek(towrite, tag_offset + 4, origin = "start", rw = "write")
  writeBin(object = as.integer(l), size = 4, con = towrite, endian = endian)
  if(l > 4) {
    # set tag value to point to end of file
    seek(towrite, tag_offset + 8, origin = "start", rw = "write")
    if(fs > 2^31-1) fs = fs - 2^32
    writeBin(object = as.integer(fs), size = 4, con = towrite, endian = endian)
    # append value to end of file
    seek(towrite, 0, origin = "end", rw = "write")
    writeBin(object = raw, con = towrite)
  } else {
    # set tag value to new value
    seek(towrite, tag_offset + 8, origin = "start", rw = "write")
    writeBin(object = raw, con = towrite)
    while(l < 4) raw = c(raw, raw(1))
  }
  return(invisible(fileName))
}

#' @title Metadata Editor
#' @description
#' Edits metadata in TIFF file
#' @param fileName path to file.
#' @param metadata string used to fill tag "270" of first IFD of 'fileName'. Default is \code{NULL}.
#' @details /!\ file will be modified /!\
#' It will only work for 1st time edition of TIFF files produced by \pkg{IFC}.
#' @return It invisibly returns path of edited file.
#' @keywords internal
editMetadata <- function(fileName, metadata = NULL) {
  IFD = IFC::getIFD(fileName)
  if(!identical(suppressWarnings(IFC::getFullTag(IFD, which = 1, tag = "270")), "N/A")) stop("'metadata' edition is not allowed")
  editTag(IFD = IFD, which = 1, tag = "270", value = metadata, force = FALSE)
}

#' @title OME XML Writer
#' @name createOME
#' @description
#' Creates OME xml
#' @param dims image dimensions [h,w,c,f]\cr
#' with h=height, w=width, c=channel, f=frame.
#' @param param object of class `IFC_param`, containing extraction parameters defined by \code{\link{objectParam}}.
#' @param name name of the image. When missing, the default it will not be incorporated.
#' @param date date of the image. Should be formatted as \verb{"\%Y:\%m:\%d \%H:\%M:\%S"}, see \code{\link[base]{format.Date}}. When missing, the default it will not be incorporated.
#' @param what bits mode used to store image. Default is \code{"int16"}. Allowed are \code{"uint8"}, \code{"int8"}, \code{"uint16"}, \code{"int16"}, \code{"uint32"}, \code{"int32"}, \code{"float"} and \code{"double"}.
#' @param endianness The endian-ness ("big" or "little") of the return object. Default is .Platform$endian.
#' @return An OME xml.
#' @keywords internal
createOME <- function(dims, param, name, date, what = "int16", endianness = .Platform$endian) {
  # FIXME ... 1
  # Should uuid be used ?
  # this would involve dependency to uuid library
  # in addition AFAIU uuid are only needed when ref to other file(s) is required which is not the case here
  
  # FIXME ... 2
  # there may be some bad interplay between <Pixels Size[XYCZT], Interleaved /> and <Channel SamplePerPixel />
  # writemulti stores [height,width,channels,frame] as frame IFD(s) with SamplePerPixel = c
  
  # FIXME ... 3
  # Is color value OK ?
  
  # FIXME ... 4
  # What should be the correct <Channel AcquisitionMode /> for ImageStream / FlowSight?
  
  # FIXME ... 5
  # pass illumination parameters ? but it is not part of 'param'
  # <Channel ExcitationWavelength /> applies to channel but with Imagestream all channels are excited by all lasers (for same camera)
  # <Channel EmissionWavelength /> but with ImageStream light beam is splitted from by filter stack to camera 
  # <LightSourceSettings ID="LightSourceID:n", PercentFraction="power ?", Wavelength="375/405/488/561/592/642/730/785" /> but there is nothing like power value in OME specification, besides same problem applies with channels and excitation as described above
  
  while(length(dims) < 4) dims = c(dims, 1)
  stopifnot(dims[3] == length(param$chan_to_keep))
  stopifnot(inherits(param, "IFC_param"))
  fun <- getFromNamespace("xml_new_node", "IFC")
  what = match.arg(arg = what, choices = c("uint8","int8","uint16","int16","uint32","int32","float","double"), several.ok = FALSE)
  SignificantBits = z = switch(what, "int8" = 8, "uint8" = 8, "int16" = 16, "uint16" = 16, "int32" = 32, "uint32" = 32, "float" = 32, "double" = 64)
  img_attr = list(ID="Image:0")
  if(!missing(name)) {
    stopifnot(length(name) == 1, typeof(name) == "character")
    img_attr = c(img_attr, list(Name=name))
  }
  date_xml = NULL
  if(!missing(date)) {
    stopifnot(length(date) == 1, typeof(date) == "character")
    date_xml = list(fun("AcquisitionDate", text = formatdate(date, iso = TRUE)))
  }
  
  # instrument is not returned by objectParam but is extracted by getInfo
  # so before running createOME one can run
  # info <- getInfo(file, ...)
  # param <- objectParam(info, ...)
  # param$instrument <- info$instrument
  inst = try({
    v = strsplit(param$instrument, "-", fixed = TRUE)[[1]]
    list(
      fun("Instrument",
          attrs = list(ID="Instrument:0"),
          .children = fun("Microscope", attrs = list("Model"=v[1], "SerialNumber"=v[2]))))
  }, silent = TRUE)
  if(inherits(inst, "try-error")) inst = NULL
  
  PhysicalSizeX=switch(as.character(param$magnification),"20"=1,"40"=0.5,"60"=0.3,1)
  PhysicalSizeY=ifelse(PhysicalSizeX == 1 && param$coremode == 0, 1, PhysicalSizeX*(param$coremode + 1))
  ome = fun(
    root = TRUE,
    "OME", 
    attrs = list(
      `xmlns`="http://www.openmicroscopy.org/Schemas/OME/2016-06",
      `xmlns:xsi`="http://www.w3.org/2001/XMLSchema-instance",
      `xsi:schemaLocation`="http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd",
      # UUID = "uuid value",                # FIXME 1
      Creator = paste0("IFC ", asNamespace("IFC")$.pkgenv$version)
    ),
    .children = c(
      inst,
      list(
        fun(
          "Image",
          attrs = img_attr,
          .children = c(
            date_xml,
            list(fun(
              "Pixels", 
              attrs = list(
                ID="Pixels:0",
                DimensionOrder="XYZTC",      # FIXME 2
                SizeX=dims[2],               # FIXME 2
                SizeY=dims[1],               # FIXME 2
                SizeC=dims[3],               # FIXME 2
                SizeT=dims[4],               # FIXME 2
                SizeZ="1",                   # FIXME 2
                PhysicalSizeX=PhysicalSizeX, # FIXME 2
                PhysicalSizeY=PhysicalSizeY, # FIXME 2
                Interleaved = "false",       # FIXME 2
                BigEndian=ifelse(endianness == "big", "true", "false"),
                SignificantBits=SignificantBits,
                Type=what
              ),
              .children = c(
                lapply(seq_along(param$chan_to_keep), FUN = function(i) {
                  chn = as.integer(param$chan_to_keep[i])
                  col = sum(as.integer(c(param$colors[[i]], 1) * 255) * c(1,256,65536,16777216)) # FIXME 3
                  if(col > 2147483647) col = col - 4294967296                                    # FIXME 3
                  fun(
                    "Channel",
                    attrs = list(
                      ID=paste0("Channel:0:",i-1),
                      Name=param$channels[chn, "name"],
                      SamplesPerPixel=length(param$chan_to_keep), # FIXME 2
                      AcquisitionMode="SpectralImaging",          # FIXME 4
                      ContrastMethod=ifelse(param$brightfield$channel[chn], "Brightfield", ifelse(any(chn %in% c(6, 12)), "Darkfield", "Fluorescence")),
                      Color=col
                    ) #,                    # FIXME 5
                    # .children = fun("LightSourceSettings", attrs = list())
                  )
                }),
                list(
                  fun("TiffData" #, lines below are optional (and not desirable since they require uuid)
                      # attrs = list(FirstC="0", FirstT="0", FirstZ="0", IFD="0", PlaneCount="1"),
                      # .children = fun("UUID", attrs = list(Filename = basename(f)), text = "uuid value") # FIXME 1
                  )
                )))))))))
  return(ome)
}
