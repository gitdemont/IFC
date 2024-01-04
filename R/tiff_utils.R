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
  bytes = length(unlist(val_raw, recursive = FALSE, use.names = FALSE))
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
    ifd = c(ifd, sapply(1:4, FUN=function(i) { ifelse(i <= bytes, unlist(val_raw, recursive = FALSE, use.names = FALSE)[i], as.raw(0x00)) }))
    add = raw()
  }
  structure(list(list(raw = as.raw(unlist(ifd, recursive = FALSE, use.names = FALSE)),
                      val = unlist(add, recursive = FALSE, use.names = FALSE),
                      byt = bytes)),
            names = tag)
}

#' @title Image Field Directory Writer
#' @description Writes Image Field Directory (IFD)
#' @param ifd an ifd extracted by cpp_fastTAGS.
#' @param r_con a connection opened for reading.
#' @param w_con a connection opened for writing.
#' @param pos current position within 'w_con'. Default is 0.
#' @param extra extra entries to add to 'ifd'. Default is NULL.
#' @param last whether ifd is last one or not.
#' @param endianness the desired endian-ness ("big" or "little"). Default is .Platform$endian.\cr
#' Endianness describes the bytes order of data stored within the files. This parameter may not be modified.
#' @return the position within 'w_con' after 'IFD' and 'extra' content have been written\cr
#' @keywords internal
writeIFD <- function(ifd, r_con, w_con, pos = 0, extra = NULL, endianness = .Platform$endian, last = FALSE, ...) {
  swap = endianness != .Platform$endian
  
  # extract image byt
  if(any("273" == names(ifd)) &&
     any("279" == names(ifd)) &&
     (ifd[["279"]]$val >= 4)) ifd[["273"]]$byt = ifd[["279"]]$val
  
  # add extra content to ifd
  ifd = c(ifd, extra)
  
  # reorder ifd
  ifd = ifd[order(as.integer(names(ifd)))]
  
  # convert modified number of entries
  n_entries = length(ifd)
  ent = cpp_uint32_to_raw(n_entries)
  if(swap) ent = rev(ent)
  
  # compute offsets of additional content
  pos = 4 + pos + 2 + n_entries * 12
  for(i_tag in seq_along(ifd)) {
    if(ifd[[i_tag]]$byt <= 4) next
    # modify ifd val/offset of minimal content 
    tmp = cpp_uint32_to_raw(pos %% 4294967296)
    if(swap) tmp = rev(tmp)
    ifd[[i_tag]]$raw[9:12] <- tmp
    pos <- pos + ifd[[i_tag]]$byt
  }
  
  # special for features offset
  if(any("33080" == names(ifd))) {
    tmp = cpp_uint32_to_raw(1)
    if(swap) tmp = rev(tmp)
    ifd[["33080"]]$raw[5:8] <- tmp
    tmp = cpp_uint32_to_raw(4)
    if(swap) tmp = rev(tmp)
    ifd[["33080"]]$raw[3:4] <- tmp[1:2]
  }
  
  # convert next offset
  if(last) {
    off = as.raw(c(0x00,0x00,0x00,0x00))
  } else {
    off = cpp_uint32_to_raw(pos %% 4294967296)
    if(swap) off = rev(off)
  }
  
  # write ifd
  writeBin(c(ent[1:2],
             unlist(lapply(seq_along(ifd), FUN = function(i_tag) ifd[[i_tag]]$raw), recursive = FALSE, use.names = FALSE),
             off), con = w_con, endian = endianness)
  
  # write additional content
  for(i_tag in seq_along(ifd)) {
    if(ifd[[i_tag]]$byt <= 4) next
    if(typeof(ifd[[i_tag]]$val) == "raw") {
      writeBin(object = ifd[[i_tag]]$val, con = w_con, endian = endianness)
    } else {
      seek(r_con, ifd[[i_tag]]$val)
      writeBin(object = readBin(con = r_con, what = "raw", n = ifd[[i_tag]]$byt), con = w_con, endian = endianness)
    }
  }
  
  # return current pos
  return(pos)
}

#' @title IFD Type Detection
#' @description Detects IFD type
#' @param ifd an Image Field Directory as extracted by cpp_fastTAGS, cpp_getTAGS or subset of getIFD.
#' @return an integer, with label attribute.
#' @keywords internal
IFDtype <- function(ifd) {
  typ = 7L
  i18 = as.numeric(ifd$tags[["33018"]]$val)
  if(length(i18) != 0) {
    l30 = ifd$tags[["33030"]]$byt
    i70 = as.numeric(ifd$tags[["33070"]]$val)
    if(length(ifd$tags[["33029"]]$val) == 0) {
      if(length(i70) == 0) {
        typ = 1L
      } else {
        if(length(l30) == 0) {
          typ = 4L
        } else {
          if(identical(i70, as.numeric(0))) {
            if(l30 == 0) {
              typ = 2L
            } else{
              typ = 3L
            }
          } else {
            if(identical(i70, i18)) {
              if(l30 == 0) {
                typ = 5L
              } else{
                typ = 6L
              }
            } else {
              typ = 8L
            }
          }
        }
      }
    } else {
      typ = 9L
    }
  }
  N = c("rif",    "sub rif",       "merged rif",
        "cif",    "cif of sub rif","cif of merged rif",
        "non xif","sub cif",       "merged cif")
  structure(typ, label = N[typ])
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
  fsize = file.size(fileName)
  if(fsize >= 2^(cpp_getBits() * 8)) stop("file is too big [",fsize,"] (more than 2^",cpp_getBits() * 8,"-1, consider using 64bits)")
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
#' @param collapse a raw vector used to collapse. Default is as.raw(0x7c)
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
