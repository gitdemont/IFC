#' @title IFC_images Raw Conversion
#' @description 
#' Helper to convert images (IFC_images object) to raw vector.
#' @param images an IFC_images object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param pb_title character string, giving the window title on the dialog box. Default is "".
#' @return an raw vector of images binaries.
#' @keywords internal
toBIN_images = function(images, endianness = .Platform$endian,
                        verbose = FALSE, display_progress = TRUE, pb_title = "") {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("exporting images information as binary")
  assert(images, cla = "IFC_images")
  assert(endianness, alw = c("little", "big"))
  assert(display_progress, alw = c(TRUE, FALSE))
  assert(pb_title, len = 1, typ = "character")
  
  L = nrow(images)
  SO_number = packBits(intToBits(L),"raw")
  bgm = grep("^bgmean", names(images))
  bgs = grep("^bgstd", names(images))
  satc = grep("^satcount", names(images))
  satp = grep("^satpercent", names(images))
  n_m = packBits(intToBits(length(bgm)),"raw")
  n_s = packBits(intToBits(length(bgs)),"raw")
  n_c = packBits(intToBits(length(satc)),"raw")
  n_p = packBits(intToBits(length(satp)),"raw")
  
  if(display_progress) {
    if(.Platform$OS.type == "windows") {
      pb = winProgressBar(title = pb_title, label = "information in %", min = 0, max = 100, initial = 1)
      pb_fun = setWinProgressBar
    } else {
      pb = txtProgressBar(title = pb_title, label = "information in %", min = 0, max = 100, initial = 1, style = 3)
      pb_fun = setTxtProgressBar
    }
    on.exit(close(pb))
    if(endianness == .Platform$endian)  {
      imgs = lapply(1:L, FUN=function(i_image) {
        k=i_image/L*100
        pb_fun(pb, value = k, label = sprintf("Converting images values %.0f%%", k))
        c(packBits(intToBits(images[i_image,"id"]),"raw"), # id
          c(packBits(intToBits(images[i_image,"imgIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00))), # imgIFD
          c(packBits(intToBits(images[i_image,"mskIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00))), # mskIFD
          c(packBits(intToBits(images[i_image,"spIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00))), # spIFD
          writeBin(images[i_image,"w"], con = raw(), endian = endianness, size = 8), # w
          writeBin(images[i_image,"l"], con = raw(), endian = endianness, size = 8), # l
          writeBin(images[i_image,"fs"], con = raw(), endian = endianness, size = 8), # fs
          packBits(intToBits(images[i_image,"cl"]),"raw"), # cl
          packBits(intToBits(images[i_image,"ct"]),"raw"), # ct
          writeBin(images[i_image,"objCenterX"], con = raw(), endian = endianness, size = 8), # objCenterX
          writeBin(images[i_image,"objCenterY"], con = raw(), endian = endianness, size = 8), # objCenterY
          n_s, # number of bgstd
          writeBin(unlist(images[i_image, bgs]), con = raw(), endian = endianness, size = 8), # bgstd
          n_m, # number of bgmean
          writeBin(unlist(images[i_image, bgm]), con = raw(), endian = endianness, size = 8), # bgmean
          n_c, # number of satcount
          writeBin(unlist(images[i_image, satc]), con = raw(), endian = endianness, size = 8), # satcount
          n_p, # number of satpercent
          writeBin(unlist(images[i_image, satp]), con = raw(), endian = endianness, size = 8)) # satpercent
      })
    } else {
      n_m = rev(n_m)
      n_s = rev(n_s)
      n_c = rev(n_c)
      n_p = rev(n_p)
      SO_number = rev(SO_number)
      imgs = lapply(1:L, FUN=function(i_image) {
        k=i_image/L*100
        pb_fun(pb, value = k, label = sprintf("Converting images values %.0f%%", k))
        c(rev(packBits(intToBits(images[i_image,"id"]),"raw")), # id
          rev(c(packBits(intToBits(images[i_image,"imgIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00)))), # imgIFD
          rev(c(packBits(intToBits(images[i_image,"mskIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00)))), # mskIFD
          rev(c(packBits(intToBits(images[i_image,"spIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00)))), # spIFD
          writeBin(images[i_image,"w"], con = raw(), endian = endianness, size = 8), # w
          writeBin(images[i_image,"l"], con = raw(), endian = endianness, size = 8), # l
          writeBin(images[i_image,"fs"], con = raw(), endian = endianness, size = 8), # fs
          rev(packBits(intToBits(images[i_image,"cl"]),"raw")), # cl
          rev(packBits(intToBits(images[i_image,"ct"]),"raw")), # ct
          writeBin(images[i_image,"objCenterX"], con = raw(), endian = endianness, size = 8), # objCenterX
          writeBin(images[i_image,"objCenterY"], con = raw(), endian = endianness, size = 8), # objCenterY
          n_s, # number of bgstd
          writeBin(unlist(images[i_image, bgs]), con = raw(), endian = endianness, size = 8), # bgstd
          n_m, # number of bgmean
          writeBin(unlist(images[i_image, bgm]), con = raw(), endian = endianness, size = 8), # bgmean
          n_c, # number of satcount
          writeBin(unlist(images[i_image, satc]), con = raw(), endian = endianness, size = 8), # satcount
          n_p, # number of satpercent
          writeBin(unlist(images[i_image, satp]), con = raw(), endian = endianness, size = 8)) # satpercent
      })
    }
  } else {
    if(endianness == .Platform$endian)  {
      imgs = lapply(1:L, FUN=function(i_image) {
        c(packBits(intToBits(images[i_image,"id"]),"raw"), # id
          c(packBits(intToBits(images[i_image,"imgIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00))), # imgIFD
          c(packBits(intToBits(images[i_image,"mskIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00))), # mskIFD
          c(packBits(intToBits(images[i_image,"spIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00))), # spIFD
          writeBin(images[i_image,"w"], con = raw(), endian = endianness, size = 8), # w
          writeBin(images[i_image,"l"], con = raw(), endian = endianness, size = 8), # l
          writeBin(images[i_image,"fs"], con = raw(), endian = endianness, size = 8), # fs
          packBits(intToBits(images[i_image,"cl"]),"raw"), # cl
          packBits(intToBits(images[i_image,"ct"]),"raw"), # ct
          writeBin(images[i_image,"objCenterX"], con = raw(), endian = endianness, size = 8), # objCenterX
          writeBin(images[i_image,"objCenterY"], con = raw(), endian = endianness, size = 8), # objCenterY
          n_s, # number of bgstd
          writeBin(unlist(images[i_image, bgs]), con = raw(), endian = endianness, size = 8), # bgstd
          n_m, # number of bgmean
          writeBin(unlist(images[i_image, bgm]), con = raw(), endian = endianness, size = 8), # bgmean
          n_c, # number of satcount
          writeBin(unlist(images[i_image, satc]), con = raw(), endian = endianness, size = 8), # satcount
          n_p, # number of satpercent
          writeBin(unlist(images[i_image, satp]), con = raw(), endian = endianness, size = 8)) # satpercent
      })
    } else {
      n_m = rev(n_m)
      n_s = rev(n_s)
      n_c = rev(n_c)
      n_p = rev(n_p)
      SO_number = rev(SO_number)
      imgs = lapply(1:L, FUN=function(i_image) {
        c(rev(packBits(intToBits(images[i_image,"id"]),"raw")), # id
          rev(c(packBits(intToBits(images[i_image,"imgIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00)))), # imgIFD
          rev(c(packBits(intToBits(images[i_image,"mskIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00)))), # mskIFD
          rev(c(packBits(intToBits(images[i_image,"spIFD"]),"raw"), as.raw(c(0x00, 0x00, 0x00, 0x00)))), # spIFD
          writeBin(images[i_image,"w"], con = raw(), endian = endianness, size = 8), # w
          writeBin(images[i_image,"l"], con = raw(), endian = endianness, size = 8), # l
          writeBin(images[i_image,"fs"], con = raw(), endian = endianness, size = 8), # fs
          rev(packBits(intToBits(images[i_image,"cl"]),"raw")), # cl
          rev(packBits(intToBits(images[i_image,"ct"]),"raw")), # ct
          writeBin(images[i_image,"objCenterX"], con = raw(), endian = endianness, size = 8), # objCenterX
          writeBin(images[i_image,"objCenterY"], con = raw(), endian = endianness, size = 8), # objCenterY
          n_s, # number of bgstd
          writeBin(unlist(images[i_image, bgs]), con = raw(), endian = endianness, size = 8), # bgstd
          n_m, # number of bgmean
          writeBin(unlist(images[i_image, bgm]), con = raw(), endian = endianness, size = 8), # bgmean
          n_c, # number of satcount
          writeBin(unlist(images[i_image, satc]), con = raw(), endian = endianness, size = 8), # satcount
          n_p, # number of satpercent
          writeBin(unlist(images[i_image, satp]), con = raw(), endian = endianness, size = 8)) # satpercent
      })
    }
  }
  return(c(SO_number, unlist(imgs)))
}
