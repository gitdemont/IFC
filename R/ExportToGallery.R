################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2025 Yohann Demont                                             #
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

resize_img = function(img, new_height = 0, new_width = 0, bg = 0) {
  d = dim(img)
  if(length(d) == 2) {
    return(cpp_resize(img, new_height = new_height, new_width = new_width, bg = bg[1], sd = 0, add_noise = FALSE))
  }
  if(length(d) == 3) {
    simplify2array(lapply(seq_len(d[3]), FUN = function(i) {
      cpp_resize(matrix(img[,,i], ncol = d[2]), new_height = new_height, new_width = new_width, bg = bg[i], sd = 0, add_noise = FALSE)
    }))
  }
}
bind_cols = function(images, add_lines = 0, color = "yellow") {
  images = images[sapply(images, length) != 0]
  d = do.call(cbind, lapply(images, dim))
  if(length(d) == 0) return(NULL)
  stopifnot(length(unique(d[2,])) == 1, nrow(d) == 3 && length(unique(d[3,])) == 1)
  v = sum(d[1,])
  if(add_lines) {
    col = col2rgb(color) / 255
    images = c(images, list(aperm(array(col, dim = c(3, d[2] * length(images) + add_lines, add_lines)), perm = c(3,2,1))))
    v = v + add_lines
  }
  aperm(array(unlist(lapply(images, aperm, perm = c(3,2,1)), use.names = FALSE), c(d[3,1], d[2,1], v)), perm = c(3,2,1))
}
bind_rows = function(images, add_lines = 0, color = "yellow") {
  images = images[sapply(images, length) != 0]
  d = do.call(cbind, lapply(images, dim))
  if(length(d) == 0) return(NULL)
  stopifnot(length(unique(d[1,])) == 1, nrow(d) == 3 && length(unique(d[3,])) == 1)
  v = sum(d[2,])
  if(add_lines) {
    col = col2rgb(color) / 255
    images = c(images, list(aperm(array(col, dim = c(3, add_lines, d[1])), perm = c(3,2,1))))
    v = v + add_lines
  }
  aperm(array(unlist(lapply(images, aperm, perm = c(3,1,2)), use.names = FALSE), c(d[3,1], d[1,1], v)), perm = c(2,3,1))
}

#' @title Sheet Creation
#' @description
#' Creates contact sheet of `IFC_img` / `IFC_msk` objects
#' @param imgs `IFC_img` / `IFC_msk` list as extracted by \code{\link{ExtractImages_toMatrix}} or \code{\link{ExtractMasks_toMatrix}} functions with \code{'mode'="rgb"}.
#' @param layout a character vector of [acquired channels + 'composite' images] members to export. Default is missing to export everything.
#' Note that members can be missing to be removed from final montage.
#' Note that members not found will be automatically removed and a warning will be thrown. 
#' @param main main title that will be displayed on top center of the image.
#' If too large it will be clipped. Default is \code{character(0)}.
#' @param grid grid size of the image montage. Default is \code{c(c(length(imgs),1)[c(byrow,!byrow)],c(length(imgs),1)[c(!byrow,byrow)])}.
#' @param byrow whether montage is filled row-wise or not. Default is \code{TRUE}.
#' @param add_channels whether to add channels names. Default is \code{TRUE}.
#' @param add_lines a length=2 finite integer vector, size of separating lines between objects. Default is \code{c(2,0)}, 1st value control horizontal space and 2nd vertical space.
#' If add_lines < 1, no separating lines are added.
#' @param add_ids integer, index of column to mark objects ids number. Default is \code{0}.
#' If add_ids < 1, no ids are added. It will also be the one use for \code{'fv'}.
#' @param order_ids whether to sort images. Default is \code{FALSE}.
#' @param scale a named list to pass to \code{addScaleBar} function. Default is \code{list()}.
#' \code{'scale'=TRUE}, can eventually be used to draw scale with default parameters.
#' @param fv feature values to display. Default is \code{NULL}. It should be a named vector whose named are objects' ids.
#' @param fv_fmt \code{'fmt'} argument to \code{\link{sprintf}} for \code{'fv'} formatting. Default is \code{"\%05.3f"}.
#' @param fv_color color for \code{'fv'}. Default is \code{"darkblue"}.
#' @param id_color color for ids. Default is \code{"yellow"}.
#' @param ln_color background color for separating lines. Default is \code{"grey20"}.
#' @param bn_color background color for main and channels. Default is \code{"grey6"}. Foreground color will be \code{"white"} or \code{"black"} depending on \code{'bn_color'} luminance.
#' @param rz_color color used when resizing is needed fo fill empty space in the montage. Default is \code{NA_character_} to use image color.
#' @return a rgb array
#' @keywords internal
CreateSheet <- function(imgs, layout, main = character(0),
                        grid = c(c(length(imgs),1)[c(!byrow,byrow)],c(length(imgs),1)[c(byrow,!byrow)]), byrow = TRUE,
                        add_channels = TRUE, add_lines = c(2,0), add_ids = 0, order_ids = FALSE,
                        scale = list(),
                        fv = NULL, fv_fmt = "%05.3f",
                        fv_color = "darkblue", id_color = "yellow", ln_color = "grey20", bn_color = "grey6", rz_color = NA_character_) {
  grid = na.omit(grid)
  add_lines = na.omit(as.integer(add_lines)); add_lines=add_lines[add_lines>=0]
  fv = unlist(fv, recursive = TRUE, use.names = TRUE)
  main = as.character(main)
  if(length(main) != 0) assert(main, len = 1, typ = "character")
  assert(add_lines, len = 2, typ = "integer")
  assert(na.omit(add_lines), len = 2)
  h_lines = ifelse(add_lines[1] < 0, 0, add_lines[1])
  v_lines = ifelse(add_lines[2] < 0, 0, add_lines[2])
  assert(order_ids, len = 1, alw = c(TRUE, FALSE))
  assert(add_channels, len = 1, alw = c(TRUE, FALSE))
  add_ids = na.omit(as.integer(add_ids))
  assert(add_ids, len = 1, typ = "integer")
  assert(grid, len = 2)
  assert(fv_color, len = 1)
  assert(id_color, len = 1)
  checkColor(ln_color)
  checkColor(bn_color)
  if(!is.na(id_color)) checkColor(id_color)
  if(!is.na(fv_color)) checkColor(fv_color)
  if(!is.na(rz_color)) fill_color = checkColor(rz_color) / 255
  
  extract_max = min(prod(grid), length(imgs))
  ids = sapply(names(imgs), FUN = function(x) as.numeric(x))
  channel_id = attr(imgs, "channel_id")
  imgs = head(imgs, extract_max)
  ids =  head(ids, extract_max)
  if(order_ids) {
    imgs = imgs[order(ids)]
    ids = ids[order(ids)]
  }
  
  # change layout
  if(missing(layout)) layout = channel_id
  layout = as.character(layout)
  mess = assert(layout, alw = channel_id, fun = "return")
  if(length(mess) != 0) warning(paste0(mess, "\n - and has been automatically removed from 'layout'"), call. = FALSE, immediate. = TRUE)
  layout = na.omit(unlist(lapply(layout, FUN = function(x) which(channel_id %in% x)), recursive = TRUE, use.names = FALSE))
  if(length(layout) == 0) stop("'layout' is of length 0 which is not allowed")
  
  # check/add object_ids
  if(add_ids > 0) if(!(add_ids %in% seq_len(length(layout)))) {
    warning("can't find 'add_ids' in 'layout'")
    add_ids = 0;
  }
  if(add_ids <= 0) {
    id_color = NA
    fv_color = NA
  }
  
  tx_color = ifelse(getLuminance(bn_color) < 128, "white", "black")
  chan_names = names(imgs[[1]][layout])
  
  BN = col2rgb(bn_color) / 255
  B = c(add_channels, length(main) != 0)
  text_height = 2 + dim(texttomatrix("foo"))[1]
  sep_b = array(BN, c(text_height * sum(B), add_lines, 3))
  
  # determine images indices in montage
  if(byrow) {
    pos = do.call(rbind, lapply(seq_len(grid[1]), FUN=function(r) {
      idx = seq_len(grid[2])+(r-1)*grid[2]
      idx[idx>extract_max] <- NA_integer_
      idx
    }))
  } else {
    pos = do.call(cbind, lapply(seq_len(grid[2]), FUN=function(r) {
      idx = seq_len(grid[1])+(r-1)*grid[1]
      idx[idx>extract_max] <- NA_integer_
      idx
    }))
  }
  
  # compute global image width W and grid lines' heights H
  W = suppressWarnings(max(sapply(imgs[na.omit(c(pos))], FUN = function(img) ifelse(length(img), dim(img[[1]])[2], NA_integer_)), na.rm = TRUE))
  H = suppressWarnings(apply(do.call(cbind, lapply(seq_len(nrow(pos)), FUN = function(i) {
    sapply(imgs[i], simplify = TRUE, FUN = function(img) {
      if(length(img)) return(dim(img[[1]])[1])
      NA_integer_
    })
  })), 2, max, na.rm = TRUE))
  
  banner = NULL
  if(any(B)) {
    banner = bind_rows(lapply(chan_names, FUN  = function(n) {
      if(!add_channels) n = ""
      addText2(image = array(BN, c(text_height * sum(B), W, 3)), text = n, color = tx_color, yoff = -2, anchor = c(0,0.5), vjust = "B", hjust = "C")
    }), v_lines, bn_color)
  }
  
  img_sprites = bind_rows(lapply(seq_len(ncol(pos)), FUN = function(j) {
    bind_cols(c(list(banner), lapply(seq_len(nrow(pos)), FUN = function(i) {
      bind_cols(list(bind_rows(lapply(seq_along(layout), FUN = function(k) {
        if(is.na(pos[i,j])) return(array(BN, c(ifelse(is.finite(H[i]), H[i], 0), W, 3)))
        object_id = names(ids[pos[i,j]])
        img = imgs[[pos[i,j]]][[layout[k]]]
        d = dim(img)
        if(d[1] != H[i] || d[2] != W) {
          if(is.na(rz_color)) {
            bg = objectNormalize(
              mat = matrix(attr(img, "BG_MEAN"),1,1),
              input_range = attr(img, "input_range"),
              full_range = attr(img, "full_range"),
              force_range = attr(img, "force_range"),
              gamma = attr(img, "gamma"))
            bg = objectColorize(mat = matrix(bg[1],1,1), color = attr(img, "color"))
          } else {
            bg = fill_color
          }
          img = resize_img(img, new_height = H[i], new_width = W, bg = bg)
        }
        if(add_ids %in% k) {
          if(!is.na(id_color)) {
            img = addText2(image = img, text = object_id, yoff = 1, xoff = 2, color = id_color, anchor = c(1,0), vjust = "T", hjust = "L")
          }
          if(object_id %in% names(fv)) {
            val = try(sprintf(fmt = fv_fmt, fv[object_id]), silent = TRUE)
            img = addText2(image = img, text = val, yoff = 1 + text_height, xoff = 2, color = fv_color, anchor = c(1,0), vjust = "T", hjust = "L")
          }
        }
        img
      }), v_lines, ln_color)), h_lines, ln_color)
    })), 0, ln_color)
  }), 0, ln_color)
  d = dim(img_sprites)
  img_sprites = img_sprites[seq_len(d[1] - h_lines), seq_len(d[2] - v_lines), 1:3]
  
  if(length(main) != 0) img_sprites = addText2(image = img_sprites, main, color = tx_color, anchor = c(1, 0.5), vjust = "T", hjust = "C")
  try({
    if(length(scale) != 0) {
      if(identical(scale, TRUE)) {
        img_sprites = do.call(what = addScaleBar, args = list(image = img_sprites))
      } else {
        scale = scale[setdiff(names(scale), c("image"))]
        if(length(scale) != 0) img_sprites = do.call(what = addScaleBar, args = c(list(image = img_sprites), scale))
      }
    }
  }, silent = TRUE)
  structure(img_sprites, object_id = ids)
}

#' @title Gallery Export
#' @description
#' Exports gallery of `IFC_img` / `IFC_msk` objects
#' @param ... arguments to be passed to \code{\link{objectExtract}} with the exception of \code{'ifd'} and \code{'bypass'}(=\strong{TRUE}).\cr
#' If \code{'param'} is provided \code{'mode'}(=\strong{"rgb"}) and the above parameters will be overwritten.\cr
#' If \code{'offsets'} are not provided extra arguments can also be passed with \code{...} \code{\link{getOffsets}}.\cr
#' \strong{/!\\} If not any of \code{'fileName'}, \code{'info'} and \code{'param'} can be found in \code{'...'} then \code{attr(offsets, "fileName_image")} will be used as \code{'fileName'} input parameter to pass to \code{\link{objectParam}}.
#' @param objects integer vector, IDEAS objects ids numbers to use.
#' This argument is not mandatory, if missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`. 
#' This argument is not mandatory but it may allow to save time for repeated image export on same file.
#' @param image_type image_type of desired offsets. Either \code{"img"} or \code{"msk"}. Default is \code{"img"}.
#' @param grid grid size of the image montage. Default is missing, resulting in a 1-row or 1-column montage depending on \code{'byrow'} argument.
#' @param byrow whether montage is filled row-wise or not. Default is \code{FALSE}.
#' @param layout a character vector of [acquired channels + 'composite' images] members to export. Default is missing to export everything.
#' Note that members can be missing to be removed from final gallery export.
#' Note that members not found will be automatically removed and a warning will be thrown.
#' @param export export format. Either "file", "matrix", "base64". Default is "matrix".
#' @param write_to used when \code{'export'} is \code{"file"} or \code{"base64"} to compute respectively filename or base64 id attribute.
#' Exported file extension and base64 MIME type will be deduced from this pattern. Allowed export are \code{".bmp"}, \code{".jpg"}, \code{".jpeg"}, \code{".png"}, \code{".tif"}, \code{".tiff"}.
#' \code{'.pdf'} can also be used when \code{'export'} is not \code{"base64"}.
#' Note that \code{".bmp"} are faster but not compressed producing bigger data.\cr
#' Placeholders, if found, will be substituted:\cr
#' -\code{\%d}: with full path directory,\cr
#' -\code{\%p}: with first parent directory,\cr
#' -\code{\%e}: with extension (without leading .),\cr
#' -\code{\%s}: with shortname (i.e. basename without extension),\cr
#' -\code{\%o}: with object_id (at most 10, will be collapse with "_", if more than one).\cr
#' -\code{\%c}: with channel_id (will be collapse with "_", if more than one, composite if any will be bracketed).
#' A good trick is to use:\cr
#' -\code{"\%d/\%s_gallery_Obj[\%o]_Ch[\%c].tiff"}, when \code{'export'} is \code{"file"}\cr
#' -\code{"\%s_gallery.bmp"}, when \code{'export'} is \code{"base64"}.\cr
#' Note that if missing and \code{'export'} is not \code{"file"}, \code{'write_to'} will be set to \code{"\%s_gallery.bmp"}.
#' @param base64_id whether to add id attribute to base64 exported object. Default is \code{FALSE}.\cr
#' Only applied when export is \code{"base64"}.
#' @param base64_att attributes to add to base64 exported object. Default is \code{""}.\cr
#' Only applied when export is \code{"base64"}. For example, use \code{"class='draggable'"}.\cr
#' Note that \code{id} (if \code{'base64_id'} is \code{TRUE}) and \code{width} and \code{height} are already used.
#' @param overwrite only apply when \code{'export'} is \code{"file"} whether to overwrite file or not. Default is \code{FALSE}.
#' @param main main title that will be displayed on top center of the image.
#' If too large it will be clipped. Default is \code{character(0)}, to produce no main banner at all.
#' @param add_channels whether to add channels names. Default is \code{TRUE}.
#' @param add_ids integer, index of column to mark objects ids number. Default is \code{1}.
#' If add_ids < 1, no ids are added.
#' @param add_lines integer vector, size of separating lines between objects, 1st value control horizontal space and 2nd vertical space.
#' If only one value is provided it will be used for both. Default is \code{2}.
#' @param bg_color background color for main, channels and separating lines. Default is \code{"grey20"}.
#' @param dpi integer, the resolution of the image in DPI (dots per inch). Default is \code{300}.\cr
#' Please note that whatever this parameter is final resolution will be 96 dpi.\cr
#' However image will be scaled according this parameter and magnification factor will be equal to this parameter divided by 96.
#' @param scale a named list whose members are \code{'size'}, \code{'style'}, \code{'color'}, \code{'xoff'}, \code{'yoff'}. Default is \code{list()} to draw no scale. Otherwise,\cr
#' -\code{'size'} positive integer. Scale's bar size in micro-meter. Default is \code{'7'}.\cr
#' -\code{'style'} a character string. Scale's bar style, either \code{"dash"} or \code{"line"}. Default is \code{"dash"}.\cr
#' -\code{'color'} a character string. color of the scale. Default is \code{"white"}.\cr
#' -\code{'xoff'} integer. x offset in image to draw scale, starting from bottom left corner.\cr
#' -\code{'yoff'} integer. y offset in image to draw scale, starting from bottom left corner.\cr
#' \code{'scale'=TRUE}, can eventually be used to draw scale with default parameters.
#' @param extract_max maximum number of objects to extract. Default is \code{10}. Use \code{+Inf} to extract all.
#' @param sampling whether to sample objects or not. Default is \code{FALSE}.
#' @param display_progress whether to display a progress bar. Default is \code{TRUE}.
#' @details arguments of \code{\link{objectExtract}} will be deduced from \code{\link{ExportToGallery}} input arguments.
#' \strong{TRICK}: for exporting only ONE \code{'objects'}, set \code{'add_channels'}(=\strong{FALSE}), \code{'add_ids'}(>=\strong{1}), \code{'force_width'}(=\strong{FALSE}), \code{'dpi'}(=\strong{96}); this allows generating image with its original size incrusted with its id number.
#' @return Depending on \code{'export'}:\cr
#' -\code{"matrix"}, a rgb array,\cr
#' -\code{"base64"}, a data-uri string,\cr
#' -\code{"file"}, invisibly returns path of exported file.\cr
#' with `object_id` attribute corresponding to the objects exported.
#' @export
ExportToGallery <- function(...,
                            objects,
                            offsets,
                            image_type = "img", 
                            layout,
                            grid,
                            byrow = TRUE,
                            export = c("file", "matrix", "base64")[2],
                            write_to, 
                            base64_id = FALSE, 
                            base64_att = "", 
                            overwrite = FALSE,
                            main = character(0), 
                            add_channels = TRUE, 
                            add_ids = 1, 
                            add_lines = 2,
                            bg_color = "grey20", 
                            dpi = 300, 
                            scale = list(),
                            extract_max = 10, 
                            sampling = FALSE, 
                            display_progress = TRUE) {
  dots = list(...)
  
  # backup last state of device ask newpage and set to FALSE
  old_ask <- devAskNewPage(ask = FALSE)
  on.exit(devAskNewPage(ask = old_ask), add = TRUE)
  # change locale
  locale_back <- setloc(c("LC_ALL" = "en_US.UTF-8"))
  enc_back <- options("encoding" = "UTF-8")
  on.exit(suspendInterrupts({setloc(locale_back); options(enc_back)}), add = TRUE)
  
  # check mandatory param
  assert(image_type, len = 1, alw = c("img", "msk"))
  assert(export, len = 1, alw = c("file", "matrix", "base64"))
  sampling = as.logical(sampling); assert(sampling, len = 1, alw = c(TRUE,FALSE))
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE,FALSE))
  base64_id = as.logical(base64_id); assert(base64_id, len = 1, alw = c(TRUE,FALSE))
  if(!missing(base64_att)) {
    base64_att = na.omit(as.character(base64_att))
    assert(base64_att, len = 1, typ = "character")
  }
  overwrite = as.logical(overwrite); assert(overwrite, len = 1, alw = c(TRUE,FALSE))
  extract_max = na.omit(as.integer(extract_max)); extract_max=extract_max[extract_max>=0]
  assert(extract_max, len = 1, typ = "integer")
  dpi = na.omit(as.integer(dpi)); dpi=dpi[dpi>=0]
  assert(dpi, len = 1, typ = "integer")
  zoom = dpi / 96
  main = as.character(main);
  if(length(main) != 0) assert(main, len = 1, typ = "character")
  add_lines = suppressWarnings(as.integer(add_lines))
  h_lines = ifelse(length(add_lines) < 1 || is.na(add_lines[1]), 0, add_lines[1])
  v_lines = ifelse(length(add_lines) < 2, h_lines, ifelse(is.na(add_lines[2]), 0, add_lines[2]))
  add_ids = na.omit(as.integer(add_ids)); assert(add_ids, len = 1, typ = "integer")
  if(missing(write_to)) stop("'write_to' can't be missing")
  
  type = getFileExt(write_to)
  if(type == "jpg") type = "jpeg"
  if(type == "tif") type = "tiff"
  alw_type =  c("bmp", "jpeg", "png", "tiff")
  if(export != "base64") alw_type = c(alw_type, "pdf")
  write_to_ = gsub(paste0(type,"$"), "bmp", write_to)
  
  # precompute param
  dots=dots[setdiff(names(dots), c("mode","export"))]
  args=list(mode = "rgb",
            export = "matrix",
            write_to = write_to_,
            overwrite = overwrite,
            base64_id = base64_id,
            base64_att = base64_att)
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  
  param = do.call(what = dotsParam, args = c(dots, args))
  fileName = param$fileName_image
  title_progress = basename(fileName)
  file_extension = getFileExt(fileName)
  
  # check objects
  nobj = as.integer(param$objcount)
  N = nchar(sprintf("%1.f",abs(nobj-1)))
  if(missing(objects)) {
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = na.omit(as.integer(objects))
    tokeep = (objects >= 0) & (objects < nobj)
    if(!all(tokeep)) {
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process:\n", paste0(objects[!tokeep], collapse=", "))
      objects = objects[tokeep]
    }
  }
  extract_max = as.integer(min(extract_max, length(objects)))
  if(sampling) {
    SEED = param$random_seed
    if(!is.list(SEED)) SEED = do.call(fetch_seed, list(seed = SEED))
    with_seed({objects=sample(objects,extract_max)}, SEED$seed, SEED$kind, SEED$normal.kind, SEED$sample.kind)
  } else {
    objects=objects[seq_len(extract_max)]
  }
  
  # check export/write_to
  overwritten = FALSE
  if(export != "matrix") {
    if(missing(write_to)) {
      if(export == "file") stop("'write_to' can't be missing")
      write_to = "%s_gallery.bmp"
    }
    assert(write_to, len = 1, typ = "character")
    assert(type, len = 1, alw = alw_type)
    splitf_obj = splitf(param$fileName_image)
    splitp_obj = splitp(write_to)
    if(length(objects) > 10) {
      obj_text = paste0(sprintf(paste0("%0",N,".f"), objects[1:10]), collapse="_")
      obj_text = paste0(obj_text, "_...")
    } else {
      obj_text = paste0(sprintf(paste0("%0",N,".f"), objects), collapse="_")
    }
    chan_text = paste0(sprintf(paste0("%0",2,".f"), as.integer(param$chan_to_keep)), collapse="_")
    if(length(param$composite) != 0) {
      comp_text = paste0("_(",gsub("/",",",param$composite),")", collapse="_")
    } else {
      comp_text = ""
    }
    write_to = formatn(splitp_obj, 
                       splitf_obj, 
                       object = obj_text,
                       channel = paste0(chan_text,comp_text))
    if(export == "file") {
      dir_name = dirname(write_to)
      if(!dir.exists(dir_name)) if(!dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dir_name))
      if(file.exists(write_to)) {
        if(!overwrite) stop(paste0("file ", write_to, " already exists"))
        write_to = normalizePath(write_to, winslash = "/")
        overwritten = TRUE
      }
      message(paste0("file will be exported in :\n", normalizePath(dirname(write_to), winslash = "/")))
    }
  }
  
  # extract images/masks
  dots=dots[setdiff(names(dots), c("param","mode","objects","display_progress"))]
  args = list(param = param, mode = "rgb", objects = objects, display_progress = display_progress)
  if(!missing(offsets)) args = c(args, list(offsets = offsets))
  if(image_type == "img") { fun = ExtractImages_toMatrix } else { fun = ExtractMasks_toMatrix }
  ans = do.call(what = fun, args = c(dots, args))
  L = length(ans)
  if(L == 0) {
    if(export == "file") {
      warning(paste0("ExportToGallery: No objects to export, check the objects you provided.\n",
                     "Can't create 'write_to' =", write_to, " from file.\n", param$fileName_image),
              immediate. = TRUE, call. = FALSE)
      return(invisible(structure(character(), object_id = structure(numeric(), names = character()))))
    } else {
      warning("ExportToGallery: No objects to export, check the objects you provided.\n", immediate. = TRUE, call. = FALSE)
      if(export == "base64") return(invisible(structure(character(), object_id = structure(numeric(), names = character()))))
      return(structure(array(NA_real_, dim = c(0,0,3)), object_id = structure(numeric(), names = character())))
    }
  }
  tryCatch({
    col = "white"
    if(missing(layout)) layout = attr(ans, "channel_id")
    if(add_ids > 0) if(any(param$brightfield$channel)) if(layout[add_ids] %in% as.character(which(param$brightfield$channel))) col = "black"
    if(!identical(scale, TRUE)) {
      scale = scale[setdiff(names(scale), "pix")]
      scale = c(scale, list("pix" = switch(as.character(param$magnification), "20"=1,"40"=0.5,"60"=0.3)))
    }
    args = list(imgs = ans, main = main, layout = layout,
                byrow = byrow, ln_color = bg_color,
                scale = scale, add_lines = c(h_lines, v_lines), add_channels = add_channels, add_ids = add_ids,
                order_ids = sampling, id_color = col)
    if(!missing(grid)) args = c(args, list(grid = grid))
    ret = do.call(CreateSheet, args = args)
    d = dim(ret)
    
    if(export == "file") {
      write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
      if(type == "pdf") {
        args = list(file = write_to, width = zoom*d[2]/96, height = zoom*d[1]/96)
      } else {
        args = list(filename = write_to, width = zoom*d[2], height = zoom*d[1], units = "px", res = 96)
      }
      do.call(what = type, args = args)
      grid.raster(image = ret, interpolate = TRUE)
      dev.off(dev.cur())
      message(paste0("\n######################\n", write_to, "\nhas been successfully ", ifelse(overwritten, "overwritten", "exported"), "\n"))
      return(invisible(structure(write_to, object_id = attr(ret, "object_id"))))
    }
    if(export == "base64") {
      if(base64_id) {
        return(structure(sprintf("<img id=%s %s width='%s' height='%s' src='data:image/%s;base64,%s'>",
                                 write_to,
                                 base64_att,
                                 ncol(ret),
                                 nrow(ret),
                                 type,
                                 cpp_base64_encode(objectWrite(x = ret, type = type, raw()))), object_id = attr(ret, "object_id")))
      } else {
        return(structure(sprintf("<img %s width='%s' height='%s' src='data:image/%s;base64,%s'>",
                                 base64_att,
                                 ncol(ret),
                                 nrow(ret),
                                 type,
                                 cpp_base64_encode(objectWrite(x = ret, type = type, raw()))), object_id = attr(ret, "object_id")))
      }
    }
    return(structure(ret, object_id = attr(ret, "object_id")))
  }, error = function(e) {
    stop(ifelse(export == "file", paste0(normalizePath(write_to, winslash = "/", mustWork = FALSE), " has been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n", e$message), e$message), call. = FALSE)
  })
}
