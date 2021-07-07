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

#' @title IFC_images Raw Conversion
#' @description 
#' Helper to convert images (`IFC_images` object) to raw vector.
#' @param images an `IFC_images` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @param ... other arguments to be passed.
#' @return an raw vector of images binaries.
#' @keywords internal
toBIN_images = function(images, endianness = .Platform$endian,
                        verbose = FALSE, display_progress = TRUE, title_progress = "", ...) {
  dots = list(...)
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("exporting images information as binary")
  assert(images, cla = "IFC_images")
  assert(endianness, alw = c("little", "big"))
  assert(display_progress, alw = c(TRUE, FALSE))
  assert(title_progress, len = 1, typ = "character")
  
  L = nrow(images)
  SO_number = cpp_uint32_to_raw(L)
  bgm = grep("^bgmean", names(images))
  bgs = grep("^bgstd", names(images))
  satc = grep("^satcount", names(images))
  satp = grep("^satpercent", names(images))
  n_m = cpp_uint32_to_raw(length(bgm))
  n_s = cpp_uint32_to_raw(length(bgs))
  n_c = cpp_uint32_to_raw(length(satc))
  n_p = cpp_uint32_to_raw(length(satp))
  extra = as.raw(c(0x00, 0x00, 0x00, 0x00))
  
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
    if(endianness == .Platform$endian)  {
      imgs = lapply(1:L, FUN=function(i_image) {
        setPB(pb, value = i_image, title = title_progress, label = "converting images values (binary)")
        c(cpp_uint32_to_raw(images[i_image,"id"]),
          c(cpp_uint32_to_raw(images[i_image,"imgIFD"]), extra),
          c(cpp_uint32_to_raw(images[i_image,"mskIFD"]), extra),
          c(cpp_uint32_to_raw(images[i_image,"spIFD"]), extra),
          writeBin(images[i_image,"w"], con = raw(), endian = endianness, size = 8), # w
          writeBin(images[i_image,"l"], con = raw(), endian = endianness, size = 8), # l
          writeBin(images[i_image,"fs"], con = raw(), endian = endianness, size = 8), # fs
          cpp_uint32_to_raw(images[i_image,"cl"]),
          cpp_uint32_to_raw(images[i_image,"ct"]),
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
        setPB(pb, value = i_image, title = title_progress, label = "converting images values (binary)")
        c(rev(cpp_uint32_to_raw(images[i_image,"id"])),
          rev(c(cpp_uint32_to_raw(images[i_image,"imgIFD"]), extra)),
          rev(c(cpp_uint32_to_raw(images[i_image,"mskIFD"]), extra)),
          rev(c(cpp_uint32_to_raw(images[i_image,"spIFD"]), extra)),
          writeBin(images[i_image,"w"], con = raw(), endian = endianness, size = 8), # w
          writeBin(images[i_image,"l"], con = raw(), endian = endianness, size = 8), # l
          writeBin(images[i_image,"fs"], con = raw(), endian = endianness, size = 8), # fs
          rev(cpp_uint32_to_raw(images[i_image,"cl"])),
          rev(cpp_uint32_to_raw(images[i_image,"ct"])),
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
        c(cpp_uint32_to_raw(images[i_image,"id"]),
          c(cpp_uint32_to_raw(images[i_image,"imgIFD"]), extra),
          c(cpp_uint32_to_raw(images[i_image,"mskIFD"]), extra),
          c(cpp_uint32_to_raw(images[i_image,"spIFD"]), extra),
          writeBin(images[i_image,"w"], con = raw(), endian = endianness, size = 8), # w
          writeBin(images[i_image,"l"], con = raw(), endian = endianness, size = 8), # l
          writeBin(images[i_image,"fs"], con = raw(), endian = endianness, size = 8), # fs
          cpp_uint32_to_raw(images[i_image,"cl"]),
          cpp_uint32_to_raw(images[i_image,"ct"]),
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
        c(rev(cpp_uint32_to_raw(images[i_image,"id"])),
          rev(c(cpp_uint32_to_raw(images[i_image,"imgIFD"]), extra)),
          rev(c(cpp_uint32_to_raw(images[i_image,"mskIFD"]), extra)),
          rev(c(cpp_uint32_to_raw(images[i_image,"spIFD"]), extra)),
          writeBin(images[i_image,"w"], con = raw(), endian = endianness, size = 8), # w
          writeBin(images[i_image,"l"], con = raw(), endian = endianness, size = 8), # l
          writeBin(images[i_image,"fs"], con = raw(), endian = endianness, size = 8), # fs
          rev(cpp_uint32_to_raw(images[i_image,"cl"])),
          rev(cpp_uint32_to_raw(images[i_image,"ct"])),
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
