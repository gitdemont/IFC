#' @title Image Field Directory Builder
#' @description build IFD
#' @param val the value of the ifd
#' @param typ desired ifd type
#' @param tag the desired ifd 'tag'
#' @param endianness the desired endian-ness ("big" or "little"). Default is .Platform$endian.
#' @details if 'val' if of type "character", 'tag' is automatically set to 2.\cr
#' if 'val' is of length 0 NULL is returned.
#' @return NULL or a list of 2 members:\cr
#' -min_content: the minimal ifd content,\cr
#' -add_content: the additional ifd content if 'val' converted to raw does not fit in 4 bytes.
#' @keywords internal
buildIFD = function(val, typ, tag, endianness = .Platform$endian) {
  if(length(val) == 0) return(NULL)
  sizes = c(1,1,2,4,4,1,1,2,4,4,4,8)
  multi = c(1,1,1,1,2,1,1,1,1,2,1,1)
  switch(typeof(val),
         "character" = { 
           typ <- 2
           val = strsplit(x = val, split = "")[[1]]
         })
  val_raw = lapply(val, FUN = function(x) {
    switch(typ,
           { x # 1 BYTE
           },
           { charToRaw(x) # 2 ASCII
           },
           { packBits(intToBits(x),type="raw")[1:2] # 3 SHORT 2 bytes, what happen when endianness is swapped ?
           },
           { packBits(intToBits(x),type="raw") # 4 LONG, 4 bytes
           },
           { packBits(intToBits(x),type="raw") # 5 RATIONAL = 2 LONG
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
  ifd = list(packBits(intToBits(tag),type="raw")[1:2], #tag
             packBits(intToBits(typ),type="raw")[1:2], #typ
             packBits(intToBits(count),type="raw")) #count
  if(endianness != .Platform$endian) {
    ifd = lapply(ifd, rev)
    val_raw = lapply(val_raw, rev)
  }
  if(bytes > 4) {
    ifd = c(ifd, packBits(intToBits(0),type="raw")) #val/offsets
    add = val_raw
  } else {
    ifd = c(ifd, sapply(4:1, FUN=function(i) { ifelse(i > bytes, as.raw(0x00), unlist(val_raw)[bytes-i+1]) }))
    add = raw()
  }
  ifd = list(list(min_content = unlist(ifd), add_content = unlist(add)))
  names(ifd) = tag
  return(ifd)
}
