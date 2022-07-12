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

#' @title Image Field Directory Builder
#' @description Builds Image Field Directory (IFD)
#' @param val the value of the IFD
#' @param typ desired IFD type
#' @param tag the desired IFD 'tag'
#' @param endianness the desired endian-ness ("big" or "little"). Default is .Platform$endian.\cr
#' Endianness describes the bytes order of data stored within the files. This parameter may not be modified.
#' @details if 'val' if of type "character", 'tag' is automatically set to 2.\cr
#' if 'val' is of length 0 NULL is returned.
#' @return NULL or a list of 2 members:\cr
#' -min_content: the minimal IFD content,\cr
#' -add_content: the additional IFD content if 'val' converted to raw does not fit in 4 bytes.
#' @keywords internal
buildIFD <- function(val, typ, tag, endianness = .Platform$endian) {
  sizes = c(1,1,2,4,4,1,1,2,4,4,4,8)
  multi = c(1,1,1,1,2,1,1,1,1,2,1,1)
  switch(typeof(val),
         "character" = { 
           typ <- 2
           val = strsplit(x = unname(val), split = character())
           if(length(val) == 1) val = val[[1]]
         })
  val_raw = lapply(unname(val), FUN = function(x) {
    switch(typ,
           { x # 1 BYTE
           },
           { 
             if(typeof(val) == "raw") {
               x
             } else {
               charToRaw(as.character(x)) # 2 ASCII
           }
           },
           { cpp_uint32_to_raw(x)[1:2] # 3 SHORT 2 bytes, what happen when endianness is swapped ?
           },
           { cpp_uint32_to_raw(x) # 4 LONG, 4 bytes
           },
           { cpp_uint32_to_raw(x) # 5 RATIONAL = 2 LONG
           },
           { x # 6 SBYTE
           },
           { x # 7 UNDEFINED, 1 Byte
           },
           { packBits(intToBits(x),type="raw")[1:2] # 8 SSHORT, 2 bytes, what happen when endianness is swapped ?
           },
           { packBits(intToBits(x),type="raw") # 9 SLONG, 4 bytes
           },
           { packBits(intToBits(x),type="raw") # 10 SRATIONAL, 2 SLONG
           },
           { writeBin(x, raw(), size = 4) # 11 FLOAT, 4 bytes
           },
           { writeBin(x, raw(), size = 8) # 12 DOUBLE, 8 bytes
           })
  })
  bytes = length(unlist(val_raw))
  count = bytes / (sizes[typ] * multi[typ])
  ifd = list(cpp_uint32_to_raw(tag)[1:2], #tag
             cpp_uint32_to_raw(typ)[1:2], #typ
             cpp_uint32_to_raw(count)) #count
  if(endianness != .Platform$endian) {
    ifd = lapply(1:length(ifd), FUN = function(i_tag) rev(ifd[[i_tag]]))
    val_raw = lapply(1:length(val_raw), FUN = function(i_tag) rev(val_raw[[i_tag]]))
  }
  if(bytes > 4) {
    ifd = c(ifd, as.raw(c(0x00, 0x00, 0x00, 0x00))) #val/offsets
    add = val_raw
  } else {
    ifd = c(ifd, sapply(1:4, FUN=function(i) { ifelse(i <= bytes, unlist(val_raw)[i], as.raw(0x00)) }))
    add = raw()
  }
  ifd = list(list(min_content = as.raw(unlist(ifd)), add_content = unlist(add), bytes = bytes))
  names(ifd) = tag
  return(ifd)
}

#' @title Image Field Directory Writer
#' @description Writes Image Field Directory (IFD)
#' @param ifd an ifd extracted by cpp_fastTAGS
#' @param r_con a connection opened for reading
#' @param w_con a connection opened for writing
#' @param pos current position within 'w_con'. Default is 0.
#' @param extra extra entries to add to 'ifd'. Default is NULL
#' @param endianness the desired endian-ness ("big" or "little"). Default is .Platform$endian.\cr
#' Endianness describes the bytes order of data stored within the files. This parameter may not be modified.
#' @return the position within 'w_con' after 'IFD' and 'extra' content have been written\cr
#' @keywords internal
writeIFD <- function(ifd, r_con, w_con, pos = 0, extra = NULL, endianness = .Platform$endian, ...) {
  swap = endianness != .Platform$endian
  
  # extract image bytes
  if(all(c("273","279") %in% names(ifd)) && !("add_content" %in% names(ifd[["273"]]))) {
    extra = c(extra, list("273" = list(min_content = ifd[["273"]]$raw, 
                                       add_content = list(typ = 1, val = ifd[["273"]]$val, len = ifd[["279"]]$val), bytes = ifd[["279"]]$val)))
    ifd = ifd[names(ifd) != "273"]
  }
  ifd = sapply(ifd, simplify = FALSE, FUN = function(i) list(min_content = i$raw,
                                                             add_content = i[c("typ","val","len")],
                                                             bytes = i$byt))
  
  # add extra content to ifd
  ifd = c(ifd, extra)
  
  # stop if duplicated names is found
  if(any(duplicated(names(ifd)))) stop("found duplicated ifd names [",paste0(names(ifd)[duplicated(names(ifd))], collpase=","),"] in ",
                                       showConnections(all = FALSE)[as.character(r_con), "description"])
  
  # reorder ifd
  ifd = ifd[order(as.integer(names(ifd)))]
  N = names(ifd)
  # define new offset position of current ifd
  l_min = cumsum(sapply(1:length(ifd), FUN=function(i_tag) length(ifd[[i_tag]]$min_content)))
  names(l_min) <- N
  l_add = cumsum(sapply(1:length(ifd), FUN=function(i_tag) ifelse(ifd[[i_tag]]$bytes > 4, ifd[[i_tag]]$bytes, 0)))
  names(l_add) <- N
  
  # modify 33080 feature offset to point to features in 33083
  if(all(c("33080","33083") %in% names(ifd))) {
    tmp = cpp_uint32_to_raw(l_add["33083"] - length(ifd[["33083"]]$add_content) + 8)
    if(swap) tmp = rev(tmp)
    ifd[["33080"]]$min_content[9:12] <- tmp
  }
  
  # write this new offset
  pos = pos + 4
  tmp = cpp_uint32_to_raw(pos + l_add[length(l_add)])
  if(swap) tmp = rev(tmp)
  writeBin(object = tmp, con = w_con, endian = endianness)
  
  # write all additional extra content
  lapply(1:length(ifd), FUN=function(i_tag) {
    if(ifd[[i_tag]]$bytes <= 4) return(NULL)
    if(inherits(x = ifd[[i_tag]]$add_content, what = "list")) {
      seek(r_con, ifd[[i_tag]]$add_content$val)
      # if(swap && (ifd[[i_tag]]$len != ifd[[i_tag]]$bytes)) {
      #   foo = readBin(con = r_con, what = "raw", n = ifd[[i_tag]]$bytes)
      #   lapply(split(foo, ceiling(seq_along(foo)/ifd[[i_tag]]$add_content$len)), FUN = function(x) {
      #     writeBin(object = rev(x), what = "raw", con = w_con, endian = endianness)
      #   })
      # } else {
        writeBin(readBin(con = r_con, what = "raw", n = ifd[[i_tag]]$byt), con = w_con, endian = endianness)
      # }
    } else {
      writeBin(ifd[[i_tag]]$add_content, con = w_con, endian = endianness)
    }
    # modify ifd val/offset of minimal content 
    tmp = cpp_uint32_to_raw(pos)
    if(swap) tmp = rev(tmp)
    ifd[[i_tag]]$min_content[9:12] <<- tmp
    pos <<- pos + ifd[[i_tag]]$bytes
  })
  # modify number of directory entries
  n_entries = length(ifd)
  tmp = cpp_uint32_to_raw(n_entries)
  if(swap) tmp = rev(tmp)
  
  # write modified number of entries
  writeBin(object = tmp[1:2], con = w_con, endian = endianness)
  
  # write ifd
  lapply(1:length(ifd), FUN=function(i_tag) {
    writeBin(object = ifd[[i_tag]]$min_content,
             con = w_con, endian = endianness)
  })
  
  # compute pos
  pos = unname(pos + l_min[length(l_min)] + 2)
  if(pos > 4294967295) stop("file is too big")
  pos
}

#' @title RIF/CIF Image Order Test
#' @description Tests order of IFD within RIF and XIF file
#' @param fileName path of file.
#' @return an integer\cr
#' -1: not a XIF file\cr
#' 0: non regular XIF file, i.e. no mask found after 1st Image itself after 1st IFD\cr
#' +1: regular XIF file, i.e. a mask is found after 1st Image itself after 1st IFD.
#' @keywords internal
testXIF <- function(fileName) {
  ans = -1L
  fileName = enc2native(fileName)
  IFD_first = getIFD(fileName = fileName, 
                     offsets = "first", 
                     trunc_bytes = 8, 
                     force_trunc = TRUE, 
                     verbose = FALSE, 
                     verbosity = 1, 
                     display_progress = FALSE,
                     bypass = TRUE)

  obj_count = suppressWarnings(as.integer(getFullTag(IFD_first, 1, "33018")))
  IFD_second = list(next_IFD_offset = 0, curr_IFD_offset = 0)
  IFD_third = list(next_IFD_offset = 0, curr_IFD_offset = 0)
  
  if(!((length(IFD_first[[1]]$infos$TYPE) == 0) || (IFD_first[[1]]$infos$TYPE != 1) || (IFD_first[[1]]$next_IFD_offset == 0))) {
    IFD_second = cpp_getTAGS(fileName, IFD_first[[1]]$next_IFD_offset, FALSE, 8, TRUE)
    if(!((length(IFD_second$infos$TYPE) == 0) || (IFD_second$infos$TYPE != 2))) {
      ans = +0L
      if(IFD_second$next_IFD_offset != 0) {
        IFD_third = cpp_getTAGS(fileName, IFD_second$next_IFD_offset, FALSE, 8, TRUE)
        if((length(IFD_third$infos$TYPE) != 0) && (IFD_third$infos$TYPE == 3)) {
          ans = +1L
        }
      }
    }
  }
  attr(ans, "obj_count") <- obj_count
  attr(ans, "obj_estimated") <- obj_count
  
  # try to evaluate number of objects in file when it can not be retrieved from tag 33018
  if((length(obj_count) == 0) || (obj_count == 0)) {
    delta_second = ifelse(IFD_second$next_IFD_offset == 0, 0, abs(IFD_second$next_IFD_offset - IFD_second$curr_IFD_offset))
    delta_third = ifelse(IFD_third$next_IFD_offset == 0, 0, abs(IFD_third$next_IFD_offset - IFD_third$curr_IFD_offset))
    delta = abs(delta_second + delta_third) / (as.integer((ans == 0)) + 1L)
    obj_estimated = ceiling(file.size(fileName) / delta)
    obj_estimated = obj_estimated[is.finite(obj_estimated)]
    if(length(obj_estimated) == 0) obj_estimated = 0
    attr(ans, "obj_count") <- 0
    attr(ans, "obj_estimated") <- obj_estimated
  }
  return(ans)
}

#' @title Raw Vectors Collapse
#' @description Collapses raw vectors together
#' @param x a list of raw vectors.
#' @collapse a raw vector used to collapse. Default is as.raw(0x7c)
#' @return a collapsed raw vector
#' @keywords internal
collapse_raw <- function(x, collapse = as.raw(0x7c)) {
  ans <- raw()
  xx = x[sapply(x, FUN = function(x_) length(x_) != 0)]
  if(length(xx) > 0) {
    for(i in seq_along(xx)) {
      ans <- c(ans, collapse, xx[[i]])
    }
    return(ans[-1])
  }
  return(ans)
}
