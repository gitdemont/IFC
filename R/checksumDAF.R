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
  images = data.frame()
  
  is_binary=as.logical(na.omit(xml_attr(xml_find_first(tmp, "//Assay"), attr = "binaryfeatures")))
  if(length(is_binary)==0) {is_binary=FALSE}
  
  if(is_binary) {
    ##### open daf fo binary extration
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
    images=lapply(1:min(5, SO_number), FUN=function(i_image) {
      id = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      imgIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used img offsets are uint32
      mskIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used msk offsets are uint32
      spIFD = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      readBin(toread, "raw", size = 1, n = 4, endian = endianness) # not used ?
      w = readBin(toread, "double", size = 8, n = 1, endian = endianness)
      l = readBin(toread, "double", size = 8, n = 1, endian = endianness)
      fs = readBin(toread, "double", size = 8, n = 1, endian = endianness)
      cl = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      ct = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      objCenterX = readBin(toread, "double", size = 8, n = 1, endian = endianness)
      objCenterY = readBin(toread, "double", size = 8, n = 1, endian = endianness)
      n_s = cpp_int32_to_uint32(readBin(toread, "integer", size = 4,  n = 1, endian = endianness))
      bgstd = readBin(toread, "double", size = 8, n = n_s)
      n_m = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      bgmean = readBin(toread, "double", size = 8, n = n_m)
      n_c = cpp_int32_to_uint32(readBin(toread, "integer", size = 4, n = 1, endian = endianness))
      satcount = readBin(toread, "double", size = 8, n = n_c)
      n_p = cpp_int32_to_uint32(readBin(toread, "integer",size = 4,  n = 1, endian = endianness))
      satpercent = readBin(toread, "double", size = 8, n = n_p)
      if(any(c(n_s, n_m, n_c, n_p) != chan_number)) stop("mismatch between expected channel count and channel values stored")
      c("id"=id,"imgIFD"=imgIFD,"mskIFD"=mskIFD,"spIFD"=spIFD,
        "w"=w,"l"=l,"fs"=fs,
        "cl"=cl,"ct"=ct,
        "objCenterX"=objCenterX,"objCenterY"=objCenterY,
        "bgstd"=bgstd,"bgmean"=bgmean,"satcount"=satcount,"satpercent"=satpercent)
    })
    images=as.data.frame(do.call(what = "rbind", args = images), stringsAsFactors = FALSE)
  } else {
    # keep on reading xml
    nodes = xml_find_all(tmp, "//SO")
    SO_number = length(nodes)
    if(SO_number != obj_count) stop("mismatch between expected object count and images numbers stored")
    images=do.call(what = cbind, args = xml_attrs(nodes[1:min(5,SO_number)]))
    img_tmp_tomodify=c("bgstd","bgmean","satcount","satpercent")
    img_tmp_tokeep=grep(paste0(img_tmp_tomodify, collapse="|"), dimnames(images)[[1]], invert = TRUE)
    img_tmp_new=lapply(img_tmp_tomodify, FUN=function(k) do.call("cbind", strsplit(images[k,],"|", useBytes = TRUE, fixed=TRUE)))
    tryCatch({
      images=rbind(images[img_tmp_tokeep,],do.call("rbind",img_tmp_new))
    }, error = function(e) {
      stop("mismatch between expected channel count and channel values stored")
    })
    images=apply(images, 1, as.numeric)
    rm(list=ls()[grep("img_tmp_",ls())])
    images=as.data.frame(images, stringsAsFactors = FALSE)
    names(images)=c("id","imgIFD","mskIFD","spIFD","w","l","fs","cl","ct","objCenterX","objCenterY",
                    paste0("bgstd",(1:chan_number)),
                    paste0("bgmean",(1:chan_number)),
                    paste0("satcount",(1:chan_number)),
                    paste0("satpercent",(1:chan_number)))
  }
  return(as.integer(sum(images[, c("imgIFD","mskIFD")])))
}
