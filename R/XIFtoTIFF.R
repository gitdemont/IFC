################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2021 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2021                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Jean-Pierre Marolleau, Loïc Garçon,       #
#                     CHU Amiens                                               #
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

#' @title RIF/CIF File Conversion to TIFF
#' @description
#' Converts and subsets RIF or CIF files to TIFF.
#' This function is experimental.
#' @param fileName path of file to subset.
#' It has to be a '.rif' or '.cif' file.
#' @param write_to pattern used to export file.
#' Placeholders, like "\%d/\%s_fromR.\%e", will be substituted:\cr
#' -\%d: with full path directory of 'fileName'\cr
#' -\%p: with first parent directory of 'fileName'\cr
#' -\%e: with extension of 'fileName' (without leading .)\cr
#' -\%s: with shortname from 'fileName' (i.e. basename without extension).\cr
#' Exported file extension has to be .tif or .tiff.
#' @param objects integer vector, IDEAS objects ids numbers to use. If missing, the default, all objects will be used.
#' @param offsets object of class `IFC_offset`. If missing, the default, offsets will be extracted from 'fileName'.\cr
#' This param is not mandatory but it may allow to save time for repeated XIF export on same file.
#' @param fast whether to fast extract 'objects' or not. Default is TRUE.
#' Meaning that 'objects' will be extracted expecting that 'objects' are stored in ascending order.\cr
#' Note that a warning will be sent if an 'object' is found at an unexpected order.
#' In such a case you may need to rerun function with 'fast' = FALSE.
#' If set to FALSE, all raw object_ids will be scanned from 'fileName' to ensure extraction of desired 'objects'.\cr
#' IMPORTANT: whatever this argument is, features are extracted assuming an ascending order of storage in file.
#' @param endianness The endian-ness ("big" or "little") of the target system for the file. Default is .Platform$endian.\cr
#' Endianness describes the bytes order of data stored within the files. This parameter may not be modified.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @param verbosity quantity of information displayed when verbose is TRUE; 1: normal, 2: rich. Default is 1.
#' @param overwrite whether to overwrite file or not. Default is FALSE.\cr
#' Note that if TRUE, it will overwrite exported file if path of 'fileName' and deduced from 'write_to' arguments are different.
#' Otherwise, you will get an error saying that overwriting source file is not allowed.\cr
#' Note also that an original file, i.e. generated by IDEAS(R) or INSPIRE(R), will never be overwritten.
#' Otherwise, you will get an error saying that overwriting original file is not allowed.\cr
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param add_tracking whether to register files' paths and objects' ids in the exported file. Default is TRUE.
#' @param ... other arguments to be passed.
#' @details This function is experimental and under development inputs and outputs may change in the future
#' @return It invisibly returns full path of exported file.
#' @keywords internal
XIFtoTIFF <- function (fileName, write_to, objects, offsets,
                       fast = TRUE,
                       endianness = .Platform$endian,
                       verbose = FALSE, verbosity = 1, 
                       overwrite = FALSE, display_progress = TRUE,
                       add_tracking = TRUE, ...) {
  dots = list(...)
  # change locale
  locale_back <- setloc(c("LC_ALL" = "en_US.UTF-8"))
  enc_back <- options("encoding" = "UTF-8")
  on.exit(suspendInterrupts({setloc(locale_back); options(enc_back)}), add = TRUE)
  
  # check madatory param
  if(missing(fileName)) stop("'fileName' can't be missing")
  tmp = duplicated(fileName)
  if(any(tmp)) {
    warning(paste0("duplicated files have been removed from 'fileName': ","\n-", paste0(fileName[tmp],collapse="\n-")))
    fileName = fileName[!tmp]
  }
  if(length(fileName) != 1) stop("'fileName' should be of length 1")
  if(!file.exists(fileName)) stop(paste("can't find", fileName, sep = " "))
  if(missing(write_to)) stop("'write_to' can't be missing")
  # extract_features = as.logical(extract_features); assert(extract_features, len = 1, alw = c(TRUE, FALSE))
  extract_features = FALSE
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  add_tracking = as.logical(add_tracking); assert(add_tracking, len = 1, alw = c(TRUE, FALSE))
  if(verbose) {
    verbosity = as.integer(verbosity); assert(verbosity, len = 1, alw = c(1,2))
    VER = ifelse(verbose & (verbosity==2), TRUE, FALSE)
  } else {
    VER = FALSE
  }
  fileName = enc2native(normalizePath(fileName, winslash = "/", mustWork = FALSE))
  r_endian = cpp_checkTIFF(fileName)
  swap = r_endian != endianness
  f_Ext = getFileExt(fileName)
  title_progress = basename(fileName)
  assert(write_to, len = 1, typ = "character")
  splitf_obj = splitf(fileName)
  splitp_obj = splitp(write_to)
  write_to = formatn(splitp_obj, splitf_obj)
  e_Ext = getFileExt(write_to); assert(e_Ext, len = 1, alw = c("tif", "tiff"))
  if(any(splitp_obj$channel > 0)) message("'write_to' has %c argument but channel information can't be retrieved with XIFtoTIFF()")
  if(any(splitp_obj$object > 0)) message("'write_to' has %o argument but channel information can't be retrieved with XIFtoTIFF()")
  dir_name = dirname(write_to)
  if(!dir.exists(dir_name)) if(!dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dir_name))
  
  overwritten = FALSE
  if(file.exists(write_to)) {
    write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
    if(!overwrite) stop(paste0("file ", write_to, " already exists"))
    if(tolower(fileName) == tolower(write_to)) stop("you are trying to overwrite source file which is not allowed")
    tryCatch({
      IFD_export = getIFD(fileName = write_to, offsets = "first", force_trunc = FALSE, trunc_bytes = 8, verbose = FALSE, bypass = FALSE)[[1]] 
    }, error = function(e) {
      stop(paste0(write_to, "\ndoes not seem to be well formatted:\n", e$message), call. = FALSE)
    })
    if(length(IFD_export$tags[["33090"]])==0) stop("you are trying to overwrite an original file which is not allowed")
    tmp_file = tempfile()
    overwritten = TRUE
  }
  file_w = ifelse(overwritten, tmp_file, write_to)
  
  # Extract important values
  IFDs = getIFD(fileName = fileName, offsets = "first", force_trunc = FALSE, trunc_bytes = 8, verbose = FALSE, bypass = FALSE)
  
  compute_offsets = TRUE
  if(!missing(offsets)) {
    if(!("IFC_offset" %in% class(offsets))) {
      warning("provided offsets do not match with expected ones, offsets will be recomputed", immediate. = TRUE, call. = FALSE)
    } else {
      if(attr(offsets, "checksum") != checksumXIF(fileName)) {
        warning("provided offsets do not match with expected ones, offsets will be recomputed", immediate. = TRUE, call. = FALSE)
      } else {
        compute_offsets = FALSE
      }
    }
  }
  if(compute_offsets) {
    offsets = suppressMessages(getOffsets(fileName = fileName, fast = fast, display_progress = display_progress)) 
  } else {
    fast = TRUE
  }
  XIF_test = ifelse(length(attr(offsets, "test")) == 0, testXIF(fileName), attr(offsets, "test"))
  XIF_step = as.integer(XIF_test == 1) + 1L
  nobj = as.integer(attr(x = offsets, which = "obj_count"))
  
  if(missing(objects)) {
    message("All objects will be extracted\n")
    objects = as.integer(0:(nobj - 1))
  } else {
    objects = as.integer(objects)
    tokeep = (objects >= 0) & (objects < nobj)
    if((length(tokeep) == 0) && (XIF_test >= 0)) {
      warning(paste0("XIFtoTIFF: No objects to export, check the objects you provided.\n",
                     "Can't create 'write_to' file.\n", write_to),
              immediate. = TRUE, call. = FALSE)
      return(invisible(NULL))
    }
    if(!all(tokeep)) {
      objects = unique(objects[tokeep]) # unique is very important
      warning("Some objects that are not in ", fileName, " have been automatically removed from extraction process")
    }
  }
  objects = objects[order(objects)] # since it appears that objects are stored with increasing number
  
  # Initialize values
  obj_id = -1
  L = length(objects)
  objects_1 = objects + 1
  offsets = subsetOffsets(offsets = offsets, objects = objects, image_type = c("img", "msk"))
  offsets = c(IFDs[[1]]$curr_IFD_offset, offsets)
  off_number = length(offsets)
  fname2 = charToRaw(fileName)
  # unwanted tags
  # TODO should we remove features and db (Inspire, Assay, ASSISTdb, dates, ..., tags 33027,33064,33078) ?
  #   258 corresponds to compression, will be overwritten with 8 or 16 depending on bit depth of uncompressed data
  #   259 corresponds to compression, will be overwritten with 1 for uncompressed data
  #   273 corresponds to strip offset, will be overwritten with uncompressed data
  #   277 corresponds to SamplesPerPixel, will be overwritten with 1
  #   279 corresponds to strip bytes count, will be overwritten with uncompressed data size
  # 33004 corresponds to file date
  # 33005 corresponds to user
  # 33018 corresponds to total object number
  # 33030 stores file path in case of merged
  # 33080 corresponds to offset of Features values, will be overwritten if features are found
  # 33081 appears in merged file, it has same val = 33080, but is of typ = 2 and map NULL 
  # 33082 corresponds to binary Features version, will be overwritten if features are found
  # 33083 corresponds to Features values in merged or subset
  # 33090, 33091, 33092, 33093, 33094 corresponds to tags we add to track objects origin
  unwanted = c(33004, 33005, 33018, 
               33027, 33064, 33078,                 # db
               33080, 33081, 33082, 33083,          # features
               33090, 33091, 33092, 33093, 33094)
  if(XIF_test >= 0) unwanted = c(258,259,273,277,279, unwanted)
  if(!add_tracking) unwanted = c(33030, unwanted)
  
  # open connections for reading and writing
  toread = file(description = fileName, open = "rb")
  on.exit(close(toread), add = TRUE)
  tryCatch(suppressWarnings({
    towrite = file(description = file_w, open = "wb")
  }), error = function(e) {
    stop(paste0(ifelse(overwritten,"temp ","'write_to' "), "file: ", file_w, "\ncan't be created: check name ?"))
  })
  write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
  tryCatch(expr = {
    # go to beginning of file
    seek(toread, 0)
    # read and writes magick number
    writeBin(object = readBin(toread, what = "raw", n = 4, endian = r_endian), con = towrite, endian = r_endian)
    # define writing position
    pos = 8
    tmp = cpp_uint32_to_raw(pos %% 4294967296)
    if(endianness != .Platform$endian) tmp = rev(tmp)
    writeBin(object = tmp, con = towrite, endian = endianness)
    
    if(display_progress) {
      pb1 = newPB(title = title_progress, label = "objects subsetting", min = 0, max = L, initial = 0, style = 3)
      on.exit(endPB(pb1), add = TRUE)
    }
    OBJECT_ID = NULL
    for(i_off in 1:off_number) {
      extra = list()
      # extract IFD
      IFD = cpp_fastTAGS(fname = fileName, offset = offsets[i_off], swap = swap)
      ifd = IFD$tags[!(names(IFD$tags) %in% unwanted)]
      ifd = ifd[!sapply(ifd, simplify = TRUE, USE.NAMES = FALSE, FUN = function(i) i$byt == 0)] # removes NULL tags
      TYPE = IFD$tags[["33002"]]$val
      if(any(TYPE == 2)) OBJECT_ID = IFD$tags[["33003"]]$val
      # if((length(TYPE) != 0) && (TYPE == 2)) OBJECT_ID = OBJECT_ID
      if(length(OBJECT_ID) == 1) {
        if(fast) {
          expected_obj = as.integer(gsub("msk_", "", gsub("img_", "", names(offsets[i_off]), fixed = TRUE), fixed = TRUE))
          if(TYPE == 2) if(OBJECT_ID != expected_obj) {
            warning("Extracted object_id [",OBJECT_ID,"] differs from expected one [",expected_obj,"]", call. = FALSE, immediate. = TRUE)
          }
        }
      }
      if(display_progress) setPB(pb1, value = (obj_id+1)/XIF_step, title = title_progress, label = "objects subsetting")
      
      # extract features values in 1st IFD
      if(i_off == 1) {
        # we need to reintroduce 258 bit depth, 259: compression, 273: strip offset, 277: SamplesPerPixel, 279: strip byte count tags
        features = list()
        if(XIF_test >= 0) {
          extra = c(extra,
                    buildIFD(val = 8, typ = 3, tag = 258, endianness = r_endian),
                    buildIFD(val = 1, typ = 3, tag = 259, endianness = r_endian),
                    buildIFD(val = 8, typ = 4, tag = 273, endianness = r_endian),
                    buildIFD(val = 1, typ = 3, tag = 277, endianness = r_endian),
                    buildIFD(val = 1, typ = 4, tag = 279, endianness = r_endian)
          )
        }
        extra = c(extra, 
                  # change now time
                  buildIFD(val = format(Sys.time(), "%d-%m-%Y %H:%M:%S %p"), typ = 2, tag = 33004, endianness = r_endian),
                  # change user
                  buildIFD(val = "IFC package", typ = 2, tag = 33005, endianness = r_endian),
                  # modify total objects count
                  buildIFD(val = L, typ = 4, tag = 33018, endianness = r_endian),
                  # add pkg version tag /!\ mandatory to prevent overwriting original file
                  buildIFD(val = paste0(unlist(packageVersion("IFC")), collapse = "."), typ = 2, tag = 33090, endianness = r_endian))
        if(add_tracking) {
          extra = c(extra,
                    # add fileName tag /!\ allow to track where exported objects are coming from
                    buildIFD(val = collapse_raw(c(list(suppressWarnings(getFullTag(IFD = structure(list(IFD), class = "IFC_ifd_list", "fileName_image" = fileName), which = 1, tag = "33091", raw = TRUE))),
                                                  list(fname2)),
                                                collapse=as.raw(0x3e)),
                             typ = 2, tag = 33091, endianness = r_endian),
                    # add checksum tag /!\ allow to track where exported objects are coming from
                    buildIFD(val = collapse_raw(c(list(suppressWarnings(getFullTag(IFD = structure(list(IFD), class = "IFC_ifd_list", "fileName_image" = fileName), which = 1, tag = "33092", raw = TRUE))),
                                                  list(charToRaw(num_to_string(checksumXIF(fileName))))),
                                                collapse=as.raw(0x3e)),
                             typ = 2, tag = 33092, endianness = r_endian))
        }
      } else {
        if(XIF_test >= 0) {
          # we need to reintroduce 258 bit depth, 259: no compression, 273: strip offset, 277: SamplesPerPixel, 279: strip byte count tags
          extra = c(extra,
                    buildIFD(val = 16, typ = 3, tag = 258, endianness = r_endian),
                    buildIFD(val = 1, typ = 3, tag = 259, endianness = r_endian),
                    buildIFD(val = 1, typ = 4, tag = 273, endianness = r_endian),
                    buildIFD(val = 1, typ = 3, tag = 277, endianness = r_endian)
          )
          # we decompress images to raw
          extra[["273"]]$val = cpp_rawdecomp(fname=fileName,
                                             offset=IFD$tags[["273"]]$val,
                                             nbytes=IFD$tags[["279"]]$val,
                                             imgWidth=IFD$tags[["256"]]$val,
                                             imgHeight=IFD$tags[["257"]]$val,
                                             compression=IFD$tags[["259"]]$val,
                                             swap = endianness!=r_endian,
                                             verbose=FALSE)
          extra = c(extra, buildIFD(val = length(extra[["273"]]$val), typ = 4, tag = 279, endianness = r_endian))
          extra[["273"]]$byt = length(extra[["273"]]$val)
          # correct FillOrder
          tmp = cpp_uint32_to_raw(1)
          if(endianness!=r_endian) tmp = rev(tmp)
          if(length(ifd[["266"]])!=0) ifd[["266"]]$raw[9:12] <- tmp
        }
        if(add_tracking) {
          extra = c(extra, 
                    # register current object id in new tag to be able to track it
                    buildIFD(val = collapse_raw(c(list(suppressWarnings(getFullTag(IFD = structure(list(IFD), class = "IFC_ifd_list", "fileName_image" = fileName), which = 1, tag = "33093", raw = TRUE))),
                                                  list(charToRaw(num_to_string(OBJECT_ID)))),
                                                collapse=as.raw(0x3e)),
                             typ = 2, tag = 33093, endianness = r_endian),
                    
                    # register current fil_ori in new tag to be able to track it
                    buildIFD(val = collapse_raw(c(list(suppressWarnings(getFullTag(IFD = structure(list(IFD), class = "IFC_ifd_list", "fileName_image" = fileName), which = 1, tag = "33094", raw = TRUE))),
                                                  list(fname2)),
                                                collapse=as.raw(0x3e)),
                             typ = 2, tag = 33094, endianness = r_endian)
          )
        }
        # TODO if we remove mask tags it may be better to use offsets names as tracking object id
        
        # modify object id
        tmp = cpp_uint32_to_raw(floor(obj_id/(XIF_step)))
        if(endianness!=r_endian) tmp = rev(tmp)
        if(length(ifd[["33003"]])!=0) ifd[["33003"]]$raw[9:12] <- tmp
        
        # TODO ask amnis what to do with 33024
      }
      pos = writeIFD(ifd, r_con = toread, w_con = towrite, pos = pos, extra = extra, endianness = r_endian, last = (i_off == off_number))
      obj_id = obj_id + 1
    }
  }, error = function(e) {
    stop(paste0("Can't create 'write_to' file.\n", write_to,
                ifelse(overwritten,"\nFile was not modified.\n","\n"),
                "See pre-file @\n", normalizePath(file_w, winslash = "/", mustWork = FALSE), "\n",
                e$message), call. = FALSE)
  }, finally = close(towrite))
  if(overwritten) {
    mess = paste0("\n######################\n", write_to, "\nhas been successfully overwritten\n")
    if(!suppressWarnings(file.rename(to = write_to, from = file_w))) { # try file renaming which is faster
      if(!file.copy(to = write_to, from = file_w, overwrite = TRUE)) { # try file.copy if renaming is not possible
        stop(paste0("Can't copy temp file@\n", normalizePath(file_w, winslash = "/"), "\n",
                    "Can't create 'write_to' file.\n", write_to,
                    "\nFile was not modified.\n"), call. = FALSE)
      } else {
        file.remove(file_w, showWarnings = FALSE)
      }
    }
  } else {
    mess = paste0("\n######################\n", write_to, "\nhas been successfully exported\n")
  }
  message(mess)
  return(invisible(write_to))
}
