#' @title DAF File Checksum 
#' @description 
#' This function returns CIF checksum computed from DAF.
#' @param fileName path to file.
#' @param endianness The endian-ness ("big" or "little") of the target system for the file. Default is .Platform$endian.
#' @keywords internal
checksumDAF <- function(fileName, endianness = .Platform$endian) {
  # TODO ask AMNIS how checksum is computed
  if(missing(fileName)) stop("'fileName' can't be missing")
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = "daf")
  if(!file.exists(fileName)) stop(paste("can't find",fileName,sep=" "))
  assert(endianness, len = 1, alw= c("big", "little"))
  fileName = normalizePath(fileName, winslash = "/", mustWork = FALSE)
  toskip=cpp_scanFirst(fname = fileName, target = "</Assay>", start = 0, end = 0)
  if(toskip == 0) stop(paste0(fileName, "\ndoes not seem to be well formatted: </Assay> not found")) 
  toskip = toskip + nchar("</Assay>") - 1
  tmp=read_xml(readBin(con = fileName, what = "raw", n = toskip), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  
  ##### extracts obj_count / chan_number
  obj_count = as.integer(xml_attr(xml_find_first(tmp, "//SOD"), attr = "objcount"))
  chan_number = as.integer(xml_attr(xml_find_first(tmp, "//ChannelPresets"), attr = "count"))
  
  ##### check binary
  is_binary=as.logical(na.omit(xml_attr(xml_find_first(tmp, "//Assay"), attr = "binaryfeatures")))
  if(length(is_binary)==0) {is_binary=FALSE}
  
  ##### always use fast = FALSE method
  fast = FALSE
  if(fast) {
    if(is_binary) {
      ##### open daf for binary extraction
      toread=file(description = fileName, open = "rb")
      on.exit(close(toread))
      
      ##### extracts important values
      seek(toread,toskip+3)
      feat_version=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      feat_number=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      obj_number=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      if(obj_count != obj_number) stop("mismatch between expected object count and features values stored")
      
      ##### extracts images values
      seek(toread, toskip + feat_number*(obj_number*8 + 4) + 15)
      SO_number=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness)) # number of SO
      if(SO_number != obj_number) stop("mismatch between expected object count and images numbers stored")
      images=sapply(1:min(5, SO_number), FUN=function(i_image) {
        id = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
        imgIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
        readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used img offsets are uint32
        mskIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
        readBin(toread, "raw", size = 1, n = 76 + 8 * 4 * chan_number, endian = endianness) # not used for checksumDAF
        c(imgIFD, mskIFD)
      })
      return(sum(images))
    } else {
      #####  keep on reading xml
      nodes = xml_find_all(tmp, "//SO")
      SO_number = length(nodes)
      if(SO_number != obj_count) stop("mismatch between expected object count and images numbers stored")
      images=do.call(what = cbind, args = xml_attrs(nodes[1:min(5,SO_number)]))
      return(sum(as.integer(images[2:3, ])))
    }
  } else {
    ##### initialize values
    obj = c(0,1,2,3,4)[1:min(5, obj_count)]
    i_image = 0
    images = c()
    if(is_binary) {
      ##### open daf for binary extraction
      toread=file(description = fileName, open = "rb")
      on.exit(close(toread))
      
      ##### extracts important values
      seek(toread,toskip+3)
      feat_version=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      feat_number=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      obj_number=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      if(obj_count != obj_number) stop("mismatch between expected object count and features values stored")
      
      ##### extracts images values
      seek(toread, toskip + feat_number*(obj_number*8 + 4) + 15)
      SO_number=cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness)) # number of SO
      if(SO_number != obj_number) stop("mismatch between expected object count and images numbers stored")
      
      ##### start extracting id + img/msk offsets
      while((length(obj) != 0) && (i_image != SO_number)) {
        i_image = i_image + 1
        id = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
        imgIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
        readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used img offsets are uint32
        mskIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
        readBin(toread, "raw", size = 1, n = 76 + 8 * 4 * chan_number, endian = endianness) # not used for checksumDAF
        if(id %in% obj) {
          obj = setdiff(obj, id)
          images = sum(images, imgIFD, mskIFD)
        } else {
          warning("raw object are stored in ascending order.")
        }
      }
    } else {
      #####  keep on reading xml
      nodes = xml_find_all(tmp, "//SO")
      SO_number = length(nodes)
      if(SO_number != obj_count) stop("mismatch between expected object count and images numbers stored")
      
      ##### loop over xml SO nodes for id and img/msk offsets
      while((length(obj) != 0) && (i_image != SO_number)) {
        i_image = i_image + 1
        val = as.integer(xml_attrs(nodes[[i_image]])[1:3])
        if(val[1] %in% obj) {
          obj = setdiff(obj, val[1])
          images = sum(images, val[2], val[3])
        } else {
          warning("raw object are stored in ascending order.")
        }
      }
    }
    return(images)
  }
}
