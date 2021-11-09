################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2021 Yohann Demont                                             #
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

#' @title FCS File Parser
#' @description
#' Parse data from Flow Cytometry Standard (FCS) compliant files.
#' @param fileName path to file.
#' @param options list of options used to parse FCS file. It should contain:\cr
#' - header, a list whose members define the "at" offset from header$start$at and the "n" number of bytes to extract:\cr
#' -- start: where start reading FCS dataset.                  Default is list(at = 0,  n = 6),\cr
#' -- text_beg: where to retrieve file text segment beginning. Default is list(at = 10, n = 8),\cr
#' -- text_end: where to retrieve file text segment end.       Default is list(at = 18, n = 8),\cr
#' -- data_beg: where to retrieve file text segment beginning. Default is list(at = 26, n = 8),\cr
#' -- data_end: where to retrieve file text segment end.       Default is list(at = 34, n = 8),\cr
#' - apply_scale, whether to apply data scaling. It only applies when fcs file is stored as DATATYPE "I". Default is TRUE.\cr
#' - first_only, whether to extract only first. Default is FALSE
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other arguments to be passed.
#' @details 'options' may be tweaked according to file type, instrument and software used to generate it.\cr
#' Default 'options' should allow to read most files.\cr
#' 'apply_scale' and 'first_only' can also be passed to 'options' thanks to ...
#' @source Data File Standard for Flow Cytometry, version FCS 3.1 from Spidlen J. et al. available at \url{https://onlinelibrary.wiley.com/doi/10.1002/cyto.a.20825}.
#' @return a list whose elements are lists for each dataset stored within the file.\cr
#' each sub-list contains:\cr
#' - header, list of header information corresponding to 'options'\cr
#' - delimiter, unique character used to separate keyword - values\cr
#' - text, list of keywords values,\cr
#' - data, data.frame of values.
#' @export
readFCS <- function(fileName, options = list(header = list(start = list(at = 0, n = 6),
                                                           text_beg = list(at = 10, n = 8),
                                                           text_end = list(at = 18, n = 8),
                                                           data_beg = list(at = 26, n = 8),
                                                           data_end = list(at = 34, n = 8)),
                                             apply_scale = TRUE,
                                             first_only = TRUE),
                    display_progress = TRUE, ...) {
  dots = list(...)
  
  if(missing(fileName)) stop("'fileName' can't be missing")
  assert(display_progress, len = 1, alw = c(TRUE,FALSE))
  assert(fileName, len = 1)
  if(!file.exists(fileName)) stop(paste("can't find",fileName,sep=" "))
  title_progress = basename(fileName)
  
  # ensure options names are valid
  if(!identical(sort(names(options)), c("apply_scale", "first_only", "header"))) stop("'options' should be a named list containing \"header\", \"apply_scale\", and \"first_only\" entries")
  if(!identical(sort(names(options$header)), c("data_beg", "data_end", "start", "text_beg", "text_end"))) stop("'options$header' should be a named list containing \"start\", \"text_beg\", \"text_end\", \"data_beg\", and \"data_end\" entries")
  if(!(all(sapply(options$header, FUN = function(x) identical(sort(names(x)), c("at","n")))))) stop("each 'options$header' members should be a named list containing \"at\" and \"n\" entries")
  
  # ensure start at is valid
  at = suppressWarnings(as.integer(options$header[["start"]]$at[1]))
  at = na.omit(at[at >=0 & at < file.size(fileName)])
  assert(at, len=1)
  
  # create connection binary reading
  toread = file(description = fileName, open = "rb")
  on.exit(close(toread), add = TRUE)
  
  # we will read offsets from options
  header = sapply(names(options$header), simplify = FALSE, FUN = function(x) {
    seek(toread, options$header[[x]]$at[1] + ifelse(x == "start", 0, at))
    raw = trimws(x = rawToChar(readBin(toread, what = "raw", n = options$header[[x]]$n)))
    if(x != "start") raw = suppressWarnings(na.omit(as.integer(raw) + at))
    raw
  })
  
  # FIXME, should we validate each header entry ?
  # e.g. header$start == FCSx.x
  # and so on ...
  
  # check if we can find options arguments in dots
  if("first_only" %in% names(dots)) options$first_only <- dots$first_only
  if("apply_scale" %in% names(dots)) options$apply_scale <- dots$apply_scale
  
  # now we can extract text segment,
  # the primary text segment has to be in within bytes 58 - 99,999,999
  off1 = header$text_beg
  off2 = header$text_end
  seek(toread, off1)
  # first byte of text segment has to be the delimiter
  delimiter = rawToChar(readBin(toread, what = "raw", n = 1))
  text = rawToChar(readBin(toread, what = "raw", n = off2 - off1))
  # when same character as delimiter is used within keyword-value pair it has to be escaped (repeated twice)
  # we generate a 20 random characters delim_esc that does not contain delimiter
  # we also ensure that this delim is not found elsewhere in the file
  found = 1
  while(found) {
    delim_esc = random_name(n = 20)
    delim_esc = strsplit(x = delim_esc, split = delimiter, fixed = TRUE)[[1]]
    delim_esc = delim_esc[delim_esc!=""]
    delim_esc = paste0(delim_esc, collapse="")
    found = cpp_scanFirst(fileName, delim_esc)
  }
  # we 1st look at double delimiter instance and substitute it with delim_esc 
  text = gsub(pattern = paste0(delimiter,delimiter), replacement = delim_esc, x = text, fixed = TRUE)
  # then text is split with delimiter
  text = strsplit(x = text, split = delimiter, fixed = TRUE)[[1]]
  # it can happen that splitting results in whitespace(s) only keyword so we remove it
  while(trimws(text[1]) == "") { text = text[-1] }
  # then escaped double delimiter is replaced with only one delimiter
  text = gsub(pattern = delim_esc, replacement = delimiter, x = text, fixed = TRUE)
  # finally keyword-value pairs are converted to named list
  id_val = seq(from = 2, to = length(text), by = 2)
  id_key = id_val-1
  text = structure(as.list(text[id_val]), names = text[id_key])
  
  # now we can extract additional text segment
  # we will use text to extract supplemental text segment offsets
  extra_off1 = suppressWarnings(na.omit(as.integer(text[["$BEGINSTEXT"]]) + at))
  extra_off2 = suppressWarnings(na.omit(as.integer(text[["$ENDSTEXT"]])   + at))
  if((length(extra_off1) != 0) &&
     (length(extra_off2) != 0) &&
     (extra_off1 != off1) && 
     (extra_off2 != off2) &&
     (extra_off2 > extra_off1)) {
    # we apply same process as for previously (see text segment)
    seek(toread, extra_off1)
    extra_text = rawToChar(readBin(toread, what = "raw", n = extra_off2 - extra_off1))
    extra_text = gsub(pattern = paste0(delimiter,delimiter), replacement = delim_esc, x = extra_text, fixed = TRUE)
    extra_text = strsplit(x = extra_text, split = delimiter, fixed = TRUE)[[1]]
    while(trimws(extra_text[1]) == "") { extra_text = extra_text[-1] }
    extra_text = gsub(pattern = delim_esc, replacement = delimiter, x = extra_text, fixed = TRUE)
    id_val = seq(from = 2, to = length(extra_text), by = 2)
    id_key = id_val-1
    extra_text = structure(as.list(extra_text[id_val]), names = extra_text[id_key])
    tmp = names(extra_text) %in% names(text)
    if(any(tmp)) warning("supplemental text segment contains keyword(s) already found in text", call. = FALSE, immediate. = TRUE)
    text = c(text, extra_text[!tmp])
  }
  if(!any("$FIL" == names(text))) text[["$FIL"]] <- fileName
  text[["@IFC_file"]] <- basename(text[["$FIL"]]) # internal filename if found, otherwise fileName
  text[["@IFC_fileName"]] <- fileName
  text[["@IFC_dataset"]] <- length(options$header$start$at)
  text[["@IFC_version"]] <- paste0(unlist(packageVersion("IFC")), collapse = ".")
  text[["@IFC_FCSversion"]] <- header$start
  
  # now we can extract data segment
  # we will use text to extract data segment offsets
  off1 = suppressWarnings(na.omit(as.integer(text[["$BEGINDATA"]]) + at))
  off2 = suppressWarnings(na.omit(as.integer(text[["$ENDDATA"]])   + at))
  # if not found in text despite being mandatory, we will use header
  if(length(off1) == 0) off1 = header$data_beg
  if(length(off2) == 0) off2 = header$data_end

  data_bytes = off2 - off1 + 1
  # prepare default returned value for data
  data = data.frame()
  if(off2 > off1) {
    seek(toread, off1)
    # retrieve info to extract data
    type = text[["$DATATYPE"]]
    if(!(type %in% c("A","I","F","D"))) stop("non-compatible data type:", type)
    mode = text[["$MODE"]]
    if(mode != "L") stop("data stored in mode[",mode,"] are not supported") # mode "C" and "U" have been deprecated in FCS spe
    n_obj = as.integer(text[["$TOT"]])
    n_par = as.integer(text[["$PAR"]])
    features_names = grep("^\\$P\\d.*N$", names(text), value = TRUE)
    
    # hereafter we create several bit_* variables
    # bit_v : PnB, bits depth of the value
    # bit_r : PnR, bits range of the value
    # bit_n : number of bytes to read 
    # bit_d : bits depth to read ( = 8 * bit_n )
    # bit_o : bytes order to read
    # bit_m : bits mask, for instance if bit_v is 10 bits but the value is read from 16 bits then 6 bits are not used
    
    bit_v = unlist(lapply(text[paste0("$P",1:n_par,"B")], FUN = function(x) suppressWarnings(as.integer(x))))
    if(n_par != length(features_names)) stop("mismatch between found vs expected number of parameters")
    if(n_par != length(bit_v)) stop("mismatch between found vs expected number of parameters")
    
    # type "A" is deprecated in newer version of FCS specifications
    if(type == "A") {
      if((data_bytes + off1) > file.size(fileName)) stop("data length points to outside of file")
      if(length(unique(bit_v)) == 1) { # no need for conversion when we have to extract same number of bytes
        data = gsub(paste0("(.{",bit_v,"})"), "\\1 ", readBin(toread, what = "character", n = data_bytes))
      } else {
        raw = readBin(toread, what = "raw", n = data_bytes)
        if(display_progress) {
          pb = newPB(session = dots$session, min = 0, max = n_par, initial = 0, style = 3)
          tryCatch({
            data = sapply(1:n_par, FUN = function(i_par) {
              setPB(pb, value = i_par, title = title_progress, label = "data-type[A]: extracting values")
              bits = as.integer(text[[paste0("$P",i_par,"B")]]) # each PnB determines number of bytes to extract
              off = (i_par - 1) * n_par
              if((off + n_obj) > data_bytes) stop("buffer overrun")
              sapply(1:n_obj, FUN = function(i_obj) {
                as.numeric(readBin(con = raw[i_obj + off], what = "character", n = bits))
              })
            })
          }, finally = endPB(pb))
        } else {
          data = sapply(1:n_par, FUN = function(i_par) {
            bits = as.integer(text[[paste0("$P",i_par,"B")]]) # each PnB determines number of bytes to extract
            off = (i_par - 1) * n_par
            if((off + n_obj) > data_bytes) stop("buffer overrun")
            sapply(1:n_obj, FUN = function(i_obj) {
              as.numeric(readBin(con = raw[i_obj + off], what = "character", n = bits))
            })
          })
        }
      }
    } else {
      # some files register wrong dataend offset resulting in an off-by-one byte
      # the following should allow to correct it
      if((data_bytes %% 8) %% 2) data_bytes = data_bytes - 1
      if((data_bytes %% 8) %% 2) data_bytes = data_bytes - 1
      if((data_bytes %% 8) %% 2) stop("number of data bytes does not respect fcs specification")
      
      # extract order
      # it is not clear from FCS spe if type == "I" can be 8, 16, 32 or 64,
      # so we choose to use byte order as the most trustable parameter to define number of bytes to extract (bit_n and bit_d)
      # however it is clearly mentioned that "F" has to be 32 and "D" 64 so we throw an error it if byteorder does not match
      b_ord = text[["$BYTEORD"]]
      bit_o = as.integer(strsplit(b_ord, split=",", fixed=TRUE)[[1]])
      bit_n = length(bit_o) 
      b_ord = paste0(bit_o, collapse = ",")
      
      endian = "unk"
      endian_l = paste0(1:bit_n, collapse = ",")
      endian_b = paste0(bit_n:1, collapse = ",")
      bit_d = 8 * bit_n
      args = list(what = "numeric", size = bit_n)
      
      # here we correct data_length, if needed
      if(n_obj * n_par * bit_n != data_bytes) warning("number of data bytes have been corrected", call. = FALSE, immediate. = TRUE)
      data_bytes = n_obj * n_par * bit_n
      if((data_bytes + off1) > file.size(fileName)) stop("data length points to outside of file")
      
      if(type == "I") { 
        # with readBin unsigned integer can only be extracted for 8bits and 16bits. So, 
        # for 32bits and 64bits we have to extract signed integers and convert them afterwards
        args = list(what = "integer", size = bit_n, signed = bit_n > 2)
        # force bits depth, shall never be > bit_d
        tmp = bit_v > bit_d
        if(any(tmp)) {
          warning(paste0("some 'PnB' keyword(s) has been forced to ", bit_d, ":\n",
                                    paste0(paste0("\t- ", names(bit_v)[tmp]), collapse = "\n")), call. = FALSE, immediate. = TRUE)
          for(i in names(bit_v)[tmp]) text[[i]] <- as.character(bit_d)
        }
        bit_v[tmp] <- bit_d
        
        # it is not clear how to perform tightbit packing
        # bit_p = xxxxx
        
        # compute bit mask
        bit_r = sapply(text[paste0("$P",1:n_par,"R")], FUN = function(x) suppressWarnings(ceiling(log2(as.integer(x)))))
        bit_m = unlist(lapply(1:n_par, FUN = function(i_par) packBits(as.raw(sapply(1:bit_d, FUN = function(i) i <= min(bit_r[i_par],bit_v[i_par]))))))
      } else {
        # fcs spe mentions:
        # No bit mask shall be applied when reading type "F" or "D" data
        if(bit_d != ifelse(type == "F", 32L, 64L)) stop("mismatch between bytes order and bits depth")
        # force bits depth, shall always be max allowed depth
        tmp = bit_v != bit_d
        if(any(tmp)) {
          warning(paste0("some 'PnB' keyword(s) has been forced to ", bit_d, ":\n",
                         paste0(paste0("\t- ", names(bit_v)[tmp]), collapse = "\n")), call. = FALSE, immediate. = TRUE)
          for(i in names(bit_v)[tmp]) text[[i]] <- as.character(bit_d)
        }
        bit_v <- bit_d 
        bit_r <- bit_d
        bit_m = NULL
      }
      if(endian_l == b_ord) endian = "little"
      if(endian_b == b_ord) endian = "big"
      
      if((endian != "unk") && all(bit_r == bit_d)) {
        # for type == "I" we are forced to apply bit masking if a PnR is not equal to bit_d
        data = do.call(what = readBin, args = c(list(con = toread, endian = endian, n = data_bytes / args$size), args))
      } else { 
        # create bit mask and order vector for all parameters
        bit_o = c(outer(bit_o, seq(from = 0, to = args$size * (n_par - 1), by = args$size), `+`))
        if(length(bit_m) != 0) # no bit_m for type F nor D
          if(length(bit_o) != length(bit_m)) stop("mismatch beetween bit order and bit mask")
        # applying bit ordering and masking is time-consuming
        if(display_progress) {
          lab = sprintf("data-type[%s]: extracting values", type)
          pb = newPB(session = dots$session, min = 0, max = n_obj, initial = 0, style = 3)
          tryCatch({
            if(length(bit_m) == 0) {
              data = sapply(1:n_obj, FUN = function(i_obj) {
                setPB(pb, value = i_obj, title = title_progress, label = lab)
                do.call(what = readBin, args = c(list(endian = "little", n = n_par, 
                                                      con = readBin(toread, what = "raw", n = n_par * args$size)[bit_o]), args))
              })
            } else {
              data = sapply(1:n_obj, FUN = function(i_obj) {
                setPB(pb, value = i_obj, title = title_progress, label = lab)
                do.call(what = readBin, args = c(list(endian = "little", n = n_par, 
                                                      con = readBin(toread, what = "raw", n = n_par * args$size)[bit_o] & bit_m), args))
              })
            }
          }, finally = endPB(pb))
        } else {
          if(length(bit_m) == 0) {
            data = sapply(1:n_obj, FUN = function(i_obj) {
              do.call(what = readBin, args = c(list(endian = "little", n = n_par, 
                                                    con = readBin(toread, what = "raw", n = n_par * args$size)[bit_o]), args))
            })
          } else {
            data = sapply(1:n_obj, FUN = function(i_obj) {
              do.call(what = readBin, args = c(list(endian = "little", n = n_par, 
                                                    con = readBin(toread, what = "raw", n = n_par * args$size)[bit_o] & bit_m), args))
            })
          }
        } 
      } 
    }
    if(type == "I") { # for 32bits and 64bits integers we need to convert signed to unsigned
      if(bit_n == 4) data = as.integer(sapply(c(data), cpp_int32_to_uint32))
      if(bit_n == 8) data = as.integer(sapply(c(data), cpp_int64_to_uint64))
    }
    
    # convert data to data.frame
    data = matrix(data, ncol = n_par, nrow = n_obj, byrow = TRUE)
    feat_names = NULL
    if(n_par > 0) feat_names = sapply(1:n_par, FUN = function(i) {
      N = text[[paste0("$P",i,"N")]]
      S = text[[paste0("$P",i,"S")]]
      if(length(S) != 0) return(paste(N , paste0("< ",S," >")))
      return(N)
    })
    data = structure(data.frame(data, check.names = FALSE), names = feat_names)
    
    # scale data only for type I, ISAC spe mentions:
    # When linear scale is used, $PnE/0,0/ shall be entered if the floating point data type is used i.e. "F" or "D"
    # meaning that no scaling shall be used for type "F" and "D". Besides type "A" is deprecated
    if(options$apply_scale && (type == "I") && (n_par > 0)) {
      sapply(1:n_par, FUN = function(i) {
        PnE = paste0("$P",i,"E")
        PnR = paste0("$P",i,"R")
        trans = text[[PnE]]
        ran = na.omit(as.numeric(text[[PnR]]))
        # FIXME should we warn user when scaling values are not valid ?
        if((length(trans) == 0) || (length(ran) == 0)) return(NULL) # no scaling info
        trans = na.omit(as.numeric(strsplit(trans, split = ",", fixed = TRUE)[[1]]))
        if((length(trans) != 2) || any(trans == 0)) return(NULL) # invalid PnE info
        data[,i] <<- 10^(trans[1] * data[,i] / ran) * trans[2]
        return(NULL)
      })
    }
  }
  
  # TODO retrieve analysis segment ?
  # # we will use text to extract analysis segment offsets
  # off1 = suppressWarnings(as.integer(text[["$BEGINANALYSIS"]]))
  # off2 = suppressWarnings(as.integer(text[["$ENDANALYSIS"]]))
  # # if not found in text despite being mandatory, we will use header
  # if(length(off1) == 0) off1 = suppressWarnings(as.integer(header$data_beg))
  # if(length(off2) == 0) off2 = suppressWarnings(as.integer(header$data_end))
  # anal=raw()
  
  ans = list(list(header=header,
                  delimiter=delimiter,
                  # anal=raw(),
                  text=text, 
                  data=data))
  
  # extract additional FCS set if any
  more = integer()
  if(!options$first_only) more = suppressWarnings(na.omit(as.integer(text[["$NEXTDATA"]]) + at))
  if((length(more) != 0) && !(more %in% options$header$start$at)) {
    options$header$start$at <- c(more, options$header$start$at)
    ans = c(ans, readFCS(fileName = fileName, options = options,
                         display_progress = display_progress, ...))
  }
  return(structure(ans, class = "IFC_fcs", fileName = fileName))
}

#' @title FCS Object Data Sets Merging
#' @description
#' Merges FCS data object with various data sets.
#' @param fcs `IFC_fcs` object as extracted by readFCS().
#' @param ... other arguments to be passed.
#' @details in data can contain extra columns named 'import_file' and 'import_subfile' intended to allow file/dataset identification
#' @return a list of list containing:\cr
#' - header, list of header information corresponding to 'options'\cr
#' - delimiter, unique character used to separate keyword - values\cr
#' - text, list of keywords values,\cr
#' - data, data.frame of values.
#' @keywords internal
FCS_merge_dataset <- function(fcs, ...) {
  dots = list(...)
  display_progress = dots$display_progress
  if(length(display_progress) == 0) display_progress = TRUE
  assert(display_progress, len=1, alw = c(TRUE, FALSE))
  
  L = length(fcs)
  if(L > 1) {
    if(display_progress) {
      pb = newPB(session = dots$session, label = "FCS", min = 0, max = L)
      on.exit(endPB(pb))
    }
    features = Reduce(function(x, y) {
      Nx = names(x)
      Ny = names(y)
      com <- Nx[Nx %in% Ny]
      Nxx <- Nx[!Nx %in% Ny]
      Nyy <- Ny[!Ny %in% Nx]
      xx = x[, Nxx, drop = FALSE]
      xx = cbind.data.frame(xx, matrix(NA, nrow = nrow(xx), ncol = length(Nyy)))
      names(xx) = c(Nxx, Nyy)
      yy = y[, Nyy, drop = FALSE]
      yy = cbind.data.frame(yy, matrix(NA, nrow = nrow(yy), ncol = length(Nxx)))
      names(yy) = c(Nyy, Nxx)
      aa = rbind.data.frame(xx[, c(Nxx, Nyy), drop = FALSE], yy[, c(Nxx, Nyy), drop = FALSE], make.row.names = FALSE)
      if(length(com) != 0) aa = structure(cbind.data.frame(aa, rbind.data.frame(x[, com, drop = FALSE],
                                                                                y[, com, drop = FALSE],
                                                                                make.row.names = FALSE)),
                                          names = c(Nxx, Nyy, com))
      aa
    },
    lapply(1:L, FUN = function(i) {
      if(display_progress) setPB(pb, value = i, title = "FCS", label = "Merging Data Sets")
      dat = fcs[[i]]$data
      if(!any("import_file" == names(dat))) dat = cbind.data.frame(dat, "import_file"=rep(fcs[[i]]$text[["@IFC_file"]], nrow(dat)))
      if(!any("import_subfile" == names(dat))) dat = cbind.data.frame(dat, "import_subfile"=rep(fcs[[i]]$text[["@IFC_dataset"]], nrow(dat)))
      dat
    }))
  } else {
    features = fcs[[1]]$data
    if(!any("import_file" == names(features))) features = cbind.data.frame(features, "import_file"=rep(fcs[[1]]$text[["@IFC_file"]], nrow(features)))
    if(!any("import_subfile" == names(features))) features = cbind.data.frame(features, "import_subfile"=rep(1, nrow(features)))
  }
  
  ans = list(list(header=fcs[[1]]$header,
                  delimiter=fcs[[1]]$delimiter,
                  text=fcs[[1]]$text, 
                  data = features))
  class(ans) <- "IFC_fcs"
  attr(ans, "fileName") <- attr(fcs, "fileName")
  bar <- unique(features[, "import_file"])
  if(length(bar) > 1) attr(ans, "Merged_fcs") <- bar
  ans
}

#' @title FCS Object Samples Merging
#' @description
#' Merges FCS data object with various samples.
#' @param fcs `IFC_fcs` object as extracted by readFCS().
#' @param ... other arguments to be passed.
#' @details in data can contain extra columns named 'import_file' and 'import_subfile' intended to allow file/dataset identification
#' @return a list of list containing:\cr
#' - header, list of header information corresponding to 'options'\cr
#' - delimiter, unique character used to separate keyword - values\cr
#' - text, list of keywords values,\cr
#' - data, data.frame of values.
#' @keywords internal
FCS_merge_sample <- function(fcs, ...) {
  dots = list(...)
  display_progress = dots$display_progress
  if(length(display_progress) == 0) display_progress = TRUE
  assert(display_progress, len=1, alw = c(TRUE, FALSE))
  
  L = length(fcs)
  if(L > 1) {
    if(display_progress) {
      pb = newPB(session = dots$session, label = "FCS", min = 0, max = L)
      on.exit(endPB(pb))
    }
    features = Reduce(function(x, y) {
      Nx = names(x)
      Ny = names(y)
      # FIXME should we add a check_names argument to ensure that each fcs sample have exactly the same names ?
      if(sum(nchar(Nx)) > sum(nchar(Ny))) {N = Nx} else {N = Ny}
      names(x) = N
      names(y) = N
      rbind.data.frame(x, y, make.row.names = FALSE)
    },
    lapply(1:L, FUN = function(i) {
      if(display_progress) setPB(pb, value = i, title = "FCS", label = "Merging Samples")
      dat = fcs[[i]]$data
      if(!any("import_file" == names(dat))) dat = cbind.data.frame(dat, "import_file"=rep(fcs[[i]]$text[["@IFC_file"]], nrow(dat)))
      if(!any("import_subfile" == names(dat))) dat = cbind.data.frame(dat, "import_subfile"=rep(fcs[[i]]$text[["@IFC_dataset"]], nrow(dat)))
      dat
    }))
  } else {
    features = fcs[[1]]$data
    if(!any("import_file" == names(features))) features = cbind.data.frame(features, "import_file"=rep(fcs[[1]]$text[["@IFC_file"]], nrow(features)))
    if(!any("import_subfile" == names(features))) features = cbind.data.frame(features, "import_subfile"=rep(1, nrow(features)))
  }
  
  ans = list(list(header=fcs[[1]]$header,
                  delimiter=fcs[[1]]$delimiter,
                  text=fcs[[1]]$text, 
                  data = features))
  class(ans) <- "IFC_fcs"
  attr(ans, "fileName") <- attr(fcs, "fileName")
  bar <- unique(features[, "import_file"])
  if(length(bar) > 1) attr(ans, "Merged_fcs") <- bar
  ans
}

#' @title Spillover Converter
#' @description
#' Converts spillover matrix to spillover keyword and reversely
#' @param spillover either a spillover matrix or a spillover keyword
#' @return if 'spillover' is a matrix, it returns a string. If 'spillover' is a string, it returns a matrix. In all cases if spillover is of length 0, it will return NULL.
#' @keywords internal
convert_spillover <- function(spillover) {
  if(length(spillover) == 0) return(NULL)
  if(is.matrix(spillover)) {
    if(length(rownames(spillover)) == 0) stop("'spillover' should have rownames")
    feat_n = parseFCSname(rownames(spillover))
    return(paste(ncol(spillover), paste0(feat_n$PnN, collapse=","), paste0(as.vector(spillover), collapse=","), sep=","))
  } else {
    foo = strsplit(spillover, split=",", fixed=TRUE)[[1]]
    feat_l = as.integer(foo[1])
    feat_n = foo[2:(feat_l+1)]
    vals = foo[-(1:(feat_l+1))]
    if(length(vals) != feat_l^2) stop("'spillover' keyword does not fulfill fcs specifications")
    return(matrix(as.numeric(vals), ncol=feat_l, nrow=feat_l, dimnames=list(NULL, feat_n)))
  }
}

#' @title FCS Object Converter
#' @description
#' Converts FCS data object to `IFC_data` object.
#' @param fcs `IFC_fcs` object as extracted by readFCS().
#' @param ... other arguments to be passed.
#' @details in data can contain extra columns named 'import_file' and 'import_subfile' intended to allow file/dataset identification
#' @return A named list of class `IFC_data`, whose members are:\cr
#' -description, a list of descriptive information,\cr
#' -Merged_fcs, character vector of path of files used to create fcs, if input was a merged,\cr
#' -fileName, path of fileName input,\cr
#' -fileName_image, path of .cif image fileName is referring to,\cr
#' -features, a data.frame of features,\cr
#' -features_def, a describing how features are defined,\cr
#' -graphs, a list of graphical elements found,\cr
#' -pops, a list describing populations found,\cr
#' -regions, a list describing how regions are defined,\cr
#' -images, a data.frame describing information about images,\cr
#' -offsets, an integer vector of images and masks IFDs offsets,\cr
#' -stats, a data.frame describing populations count and percentage to parent and total population,\cr
#' -checksum, a checksum integer.
#' @keywords internal
FCS_to_data <- function(fcs, ...) {
  # create structure
  dots = list(...)
  display_progress = dots$display_progress
  if(length(display_progress) == 0) display_progress = TRUE
  assert(display_progress, len=1, alw = c(TRUE, FALSE))
  min_data = list("description" = list("Assay" = data.frame("date" = NULL, "IDEAS_version" = NULL, "binaryfeatures" = NULL),
                                       "ID" = data.frame("file" = NULL, "creation" = NULL, "objcount" = NULL, "checksum" = NULL),
                                       "Images" = data.frame("name" = NULL, "color" = NULL, "physicalChannel" = NULL, "xmin" = NULL,
                                                             "xmax" = NULL, "xmid" = NULL, "ymid"= NULL, "scalemin"= NULL, "scalemax"= NULL,
                                                             "tokens"= NULL, "baseimage"= NULL, "function"= NULL),
                                       "masks" = data.frame()),
                  "Merged_fcs" = character(),
                  "fileName" = character(),
                  "fileName_image" = character(),
                  "features" = structure(.Data = list(), class = c("data.frame", "IFC_features")),
                  "features_def" = structure(.Data =  list(), class = c("IFC_features_def")),
                  "graphs" = structure(.Data =  list(), class = c("IFC_graphs")),
                  "pops" = structure(.Data =  list(), class = c("IFC_pops")),
                  "regions" = structure(.Data = list(), class = c("IFC_regions")),
                  "images" = structure(.Data = list(), class = c("data.frame", "IFC_images")),
                  "offsets" = structure(.Data = integer(), class = c("IFC_offsets")),
                  "stats" = data.frame(),
                  "checksum" = integer())
  class(min_data) = c("IFC_data")
  
  # define features categories which requires no param
  No_Param = c("Time", "Object Number", "Raw Centroid X", "Raw Centroid Y",  "Flow Speed", "Camera Line Number", "Camera Timer", "Objects per mL", "Objects per sec")
  
  features = FCS_merge_dataset(fcs, ...)[[1]]$data
  
  identif = names(features) %in% c("import_file", "import_subfile")
  idx = features[, identif, drop = FALSE]
  if(!"import_file" %in% names(idx)) idx$import_file = fcs[[1]][["@IFC_file"]]
  if(!"import_subfile" %in% names(idx)) idx$import_subfile = 1
  obj_count = as.integer(nrow(features))
  
  multiple = prod(length(unique(idx[, 1])), length(unique(idx[, 2]))) > 1
  # if several files creates pops to identify who is who
  if(multiple) {
    idx$count = 1:obj_count
    all_obj = rep(FALSE, obj_count)
    pops = by(idx, idx[, c("import_file", "import_subfile")], FUN = function(x) {
      if(length(unique(idx[, "import_subfile"])) == 1) {
        name = unique(x$import_file)
      } else {
        name = paste(unique(x$import_file), "dataset", unique(x$import_subfile), sep = "_")
      }
      obj = all_obj
      obj[x$count] <- TRUE
      buildPopulation(name = name, type = "T", color = "White", lightModeColor = "Black", obj = obj)
    })
    pops = pops[sapply(pops, FUN = function(p) length(p) != 0)]
  }
  
  features = subset(x = features, select = !identif)
  features_def = lapply(names(features), FUN = function(i_feat) {
    # TODO check it it correctly imports linear values
    if(i_feat %in% No_Param) return(buildFeature(name = gsub("LOG$", "LIN", i_feat, ignore.case = TRUE), val = features[, i_feat], def = i_feat))
    return(buildFeature(name = gsub("LOG$", "LIN", i_feat, ignore.case = TRUE), val = features[, i_feat]))
  })
  
  # fill min object
  instrument = sapply(fcs, FUN = function(x) {
    tmp = x$text$`$CYT`
    if((length(tmp) == 0) || (tmp == "")) return("unk")
    return(tmp)
  })
  FCS_version = sapply(fcs, FUN = function(x) {
    tmp = x$text[["@IFC_FCSversion"]]
    if((length(tmp) == 0) || (tmp == "")) return("unk")
    return(tmp)
  })
  spillover = lapply(fcs, FUN = function(x) {
    tmp = x$text[c("$SPILLOVER","SPILL","spillover")]
    tmp = tmp[sapply(tmp, length) != 0]
    if(length(tmp) == 0) return(NULL)
    return(tmp)
  })
  if(is.list(spillover) && length(spillover) == 1) spillover = spillover[[1]]
  if(is.list(spillover) && length(spillover) == 1) spillover = spillover[[1]]
  if(!is.list(spillover) && length(spillover) == 1 && spillover == "") spillover = NULL
  if(length(spillover) != 0) {
    features_names = parseFCSname(names(features))
    spillover = convert_spillover(spillover)
    rownames(spillover) <- names(features)[apply(sapply(colnames(spillover), FUN = function(x) {
      x==features_names$PnN
    }), 2, FUN = function(i) which(i)[1])]
  }
  # checksum = sapply(fcs, FUN = function(x) {
  #   tmp = x$description[[1]]$`$ENDDATA`
  #   if((length(tmp) == 0) || (tmp == "")) return("unk")
  #   return(tmp)
  # })
  
  min_data$fileName = normalizePath(attr(fcs, "fileName"), winslash = "/", mustWork = FALSE)
  bar <- unique(idx[, "import_file"])
  if(length(bar) > 1) min_data$Merged_fcs <- bar
  min_data$description$Assay = data.frame(date = file.info(min_data$fileName)$mtime, FCS_version = paste0(FCS_version, collapse = ", "), stringsAsFactors = FALSE)
  min_data$description$ID = data.frame(file = min_data$fileName,
                                       creation = format(file.info(min_data$fileName)$ctime, format = "%d-%b-%y %H:%M:%S"),
                                       objcount = obj_count,
                                       stringsAsFactors = FALSE)
  min_data$description$FCS = min_data$description$ID
  min_data$checksum = integer()
  min_data$features = structure(data.frame("Object Number" = 0:(obj_count-1), check.names = FALSE), class = c("data.frame", "IFC_features"))
  min_data$features_def = structure(list(buildFeature(name = "Object Number", val = 0:(obj_count-1), def = "Object Number")[1:4]), names = "Object Number", class = c("list", "IFC_features_def"))
  # foo = grep("^\\$P|^\\@P|^\\$D|^@D|^flowCore", names(fcs@description), value = TRUE, invert = TRUE)
  # min_data$info = fcs@description[foo]
  min_data$description$FCS = c(min_data$description$ID, list(instrument = paste0(instrument, collapse = ", "), spillover = spillover))
  min_data = suppressWarnings(IFC::data_add_features(obj = min_data, features = features_def, session = dots$session))
  min_data = IFC::data_add_pops(obj = min_data,
                                pops = list(buildPopulation(name = "All", type = "B",
                                                            color = "White", lightModeColor = "Black",
                                                            obj = rep(TRUE, obj_count))),
                                session = dots$session)
  # min_data$features_comp = min_data$features[, grep("^FS.*$|^SS.*$|LOG|^Object Number$|TIME", names(min_data$features), value = TRUE, invert = TRUE, ignore.case = TRUE)]
  if(multiple) {
    min_data = IFC::data_add_pops(obj = min_data, pops = pops, session = dots$session)
  }
  K = class(min_data$pops)
  min_data$pops = lapply(min_data$pops, FUN = function(p) {
    attr(p, "reserved") <- TRUE
    return(p)
  })
  class(min_data$pops) <- K
  pops = min_data$pops
  stats = data.frame(stringsAsFactors = FALSE,
                     check.rows = FALSE,
                     check.names = FALSE,
                     t(sapply(names(pops),
                              FUN = function(p) {
                                count = sum(pops[[p]]$obj)
                                base = pops[[p]]$base
                                type = pops[[p]]$type
                                if (base == "") base = "All"
                                parent = sum(pops[[base]]$obj)
                                c(type = type, parent = base, count = count,
                                  perc_parent = count/parent * 100,
                                  perc_tot = count/obj_count * 100)
                              })))
  stats[, 3] = as.numeric(stats[, 3])
  stats[, 4] = as.numeric(stats[, 4])
  stats[, 5] = as.numeric(stats[, 5])
  min_data$stats = stats
  return(min_data)
}

#' @title FCS File Reader
#' @description
#' Extracts data from Flow Cytometry Standard (FCS) Files.
#' @param fileName path(s) of file(s). If multiple files are provided they will be merged and 
#' populations will be created to identify each single file within returned `IFC_data` object.
#' @source Data File Standard for Flow Cytometry, version FCS 3.1 from Spidlen J. et al. available at \url{https://onlinelibrary.wiley.com/doi/10.1002/cyto.a.20825}.
#' @param ... other arguments to be passed to readFCS function.
#' @return A named list of class `IFC_data`, whose members are:\cr
#' -description, a list of descriptive information,\cr
#' -Merged_fcs, character vector of path of files used to create fcs, if input was a merged,\cr
#' -fileName, path of fileName input,\cr
#' -fileName_image, path of .cif image fileName is referring to,\cr
#' -features, a data.frame of features,\cr
#' -features_def, a describing how features are defined,\cr
#' -graphs, a list of graphical elements found,\cr
#' -pops, a list describing populations found,\cr
#' -regions, a list describing how regions are defined,\cr
#' -images, a data.frame describing information about images,\cr
#' -offsets, an integer vector of images and masks IFDs offsets,\cr
#' -stats, a data.frame describing populations count and percentage to parent and total population,\cr
#' -checksum, a checksum integer.
#' @export
ExtractFromFCS <- function(fileName, ...) {
  # create structure
  dots = list(...)
  display_progress = dots$display_progress
  if(length(display_progress) == 0) display_progress = TRUE
  assert(display_progress, len=1, alw = c(TRUE, FALSE))
  fileName = normalizePath(path = fileName, winslash = "/", mustWork = TRUE)
  
  # read the fcs file and extract features and description
  L = length(fileName)
  if(display_progress) {
    pb = newPB(session = dots$session, label = "reading files", min = 0, max = L)
    on.exit(endPB(pb))
  }
  fcs = lapply(1:L, FUN = function(i_file) {
    if(display_progress) setPB(pb, value = i_file, title = "Extracting FCS", label = "reading files")
    FCS_merge_dataset(readFCS(fileName = fileName[[i_file]], ...))[[1]]
  })
  attr(fcs, "fileName") <- fileName[1]
  FCS_to_data(fcs, ...)
}

#' @title FCS File Writer
#' @description
#' Writes an `IFC_data` object to a Flow Cytometry Standard (FCS) file.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param write_to pattern used to export file.
#' Placeholders, like "\%d/\%s_fromR.\%e", will be substituted:\cr
#' -\%d: with full path directory of 'obj$fileName'\cr
#' -\%p: with first parent directory of 'obj$fileName'\cr
#' -\%e: with extension of 'obj$fileName' (without leading .)\cr
#' -\%s: with shortname from 'obj$fileName' (i.e. basename without extension).\cr
#' Exported file extension will be deduced from this pattern. Note that it has to be a .fcs.
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' Note that if TRUE, it will overwrite exported file if path of 'fileName' and deduced from 'write_to' arguments are different.
#' Otherwise, you will get an error saying that overwriting source file is not allowed.\cr
#' Note also that an original file will never be overwritten.
#' @param delimiter an ASCII character to separate the FCS keyword-value pairs. Default is : "/".
#' @param cytometer string, if provided it will be used to fill $CYT keyword.\cr
#' However, when missing $CYT will be filled with obj$description$FCS$instrument if found, or "Image Stream" otherwise.
#' @param ... other arguments to be passed. keyword-value pairs can be passed here.
#' @return invisibly returns full path of exported file.
#' @export
ExportToFCS <- function(obj, write_to, overwrite = FALSE, delimiter="/", cytometer = "Image Stream", ...) {
  dots = list(...)
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  now = format(Sys.time(), format = "%d-%b-%y %H:%M:%S")
  
  # check mandatory param
  assert(obj, cla = "IFC_data")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  if(missing(write_to)) stop("'write_to' can't be missing")
  assert(write_to, len = 1, typ = "character")
  assert(overwrite, len = 1, alw = c(TRUE, FALSE))
  assert(cytometer, len = 1, typ = "character")
  # assert(display_progress, c(TRUE, FALSE))
  raw_delimiter = charToRaw(delimiter)
  if(length(raw_delimiter) != 1) stop("'delimiter' should be of size 1")
  if(readBin(raw_delimiter, what = "int", signed = FALSE, size = 1) > 127) stop("'delimiter' should be an ASCII character")
  delimiter_esc = paste0(delimiter, delimiter)
  
  # tests if file can be written
  fileName = normalizePath(obj$fileName, winslash = "/", mustWork = FALSE)
  title_progress = basename(fileName)
  splitf_obj = splitf(fileName)
  splitp_obj = splitp(write_to)
  write_to = formatn(splitp_obj, splitf_obj)
  file_extension = getFileExt(write_to)
  assert(file_extension, len = 1, alw = "fcs")
  if(any(splitp_obj$channel > 0)) message("'write_to' has %c argument but channel information can't be retrieved with data_to_DAF()")
  if(any(splitp_obj$object > 0)) message("'write_to' has %o argument but channel information can't be retrieved with data_to_DAF()")
  
  overwritten = FALSE
  if(file.exists(write_to)) {
    write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
    if(!overwrite) stop(paste0("file ",write_to," already exists"))
    if(tolower(fileName) == tolower(write_to)) stop("you are trying to overwrite source file which is not allowed")
    tryCatch({
      fcs = readFCS(fileName = write_to)
    }, error = function(e) {
      stop(paste0(write_to, "\ndoes not seem to be well formatted")) 
    })
    if(length(fcs[[1]]$text[["@IFC_version"]]) == 0) stop("you are trying to overwrite an original file which is not allowed")
    tmp_file = normalizePath(tempfile(), winslash = "/", mustWork = FALSE)
    overwritten = TRUE
  }
  
  dir_name = dirname(write_to)
  if(!dir.exists(dir_name)) if(!dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dir_name))
  file_w = ifelse(overwritten, tmp_file, write_to)
  tryCatch(suppressWarnings({
    towrite = file(description = file_w, open = "wb")
  }), error = function(e) {
    stop(paste0(ifelse(overwritten,"temp ","'write_to' "), "file: ", file_w, "\ncan't be created: check name ?"))
  })
  close(towrite)
  write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
  
  # defines some variables
  pkg_ver = paste0(unlist(packageVersion("IFC")), collapse = ".")
  # is_fcs = length(obj$description$FCS)!=0
  
  # init header
  header = list(version  = "FCS3.0",
                space    =  "    ",
                text_beg = "      58",
                text_end = character(),
                data_beg = character(),
                data_end = character(),
                anal_beg = "       0",
                anal_end = "       0")
  
  # we modify features to add populations
  all_pops = do.call(what = cbind, args = lapply(obj$pops, FUN = function(p) p$obj))
  colnames(all_pops) = names(obj$pops)
  features = cbind(obj$features[, setdiff(names(obj$features), colnames(all_pops))], all_pops)
  # need to replace non finite values by something; IDEAS is using 0 so we use 0 also
  # TODO maybe replace -Inf by features min and +Inf by features max ?
  features = as.data.frame(apply(features, 2, FUN = function(x) {x[!is.finite(x)] <- 0; x}), stringsAsFactors = TRUE)
  
  # determines length of text_segment2
  feat_names = parseFCSname(names(features))
  N = feat_names$PnN
  S = feat_names$PnS
  L = length(N)
  text2_length = 0
  text_segment2 = lapply(1:L, FUN = function(i) {
    # each time we write /$PnN/<feature_name>/$PnB/32/$PnE/0, 0/$PnR/<feature_max_value> and /$PnS/<feature_alt-name> if PnS is not empty
    # since we write type "F" this has no importance
    bar = max(features[, i], na.rm = TRUE)
    # if(ceiling(bar) == bar) bar = bar + 1
    # FIXME shall we use 262144, as it is commonly used ?
    foo = c("", paste0("$P",i,"N"), N[i], paste0("$P",i,"B"), "32", paste0("$P",i,"E"), "0, 0", paste0("$P",i,"R"), ceiling(bar))
    if(S[i] != "") foo = c(foo, paste0("$P",i,"S"), S[i])
    foo = gsub(pattern = delimiter, x = foo, replacement = delimiter_esc, fixed=TRUE)
    foo = charToRaw(paste(foo, collapse = delimiter))
    text2_length <<- text2_length + length(foo)
    return(foo)
  })
  cyt = obj$description$FCS$instrument
  if((length(cyt) == 0 ) || (cyt == "")) cyt = "Image Stream"
  if(!missing(cytometer)) cyt = cytometer
  if(length(obj$fileName_image) == 0) obj$fileName_image = ""
  
  # init text segment with mandatory + custom parameters #* = mandatory
  text_segment1 = list("$BEGINSTEXT" = "0",                                                      #*
                       "$ENDSTEXT" = "0",                                                        #*
                       "$BEGINANALYSIS" = "0",                                                   #*
                       "$ENDANALYSIS" = "0",                                                     #*
                       "$BYTEORD" = c("4,3,2,1", "1,2,3,4")[(.Platform$endian == "little") + 1], #*
                       "$DATATYPE" = "F",                                                        #*
                       "$MODE" = "L",                                                            #*
                       "$NEXTDATA" = "0",                                                        #*
                       "$TOT" = obj$description$ID$objcount,                                     #*
                       "$PAR" = L,                                                               #*
                       #* BEGINDATA and ENDDATA are also mandatory and will be added afterwards
                       #* PnB, PnE, PnN, and PnR are also mandatory and are part of text_segment2
                       "$CYT" = cyt,
                       "@IFC_fileName" = obj$fileName,
                       "@IFC_fileName_image" = obj$fileName_image,
                       "@IFC_version" = pkg_ver,
                       "@IFC_date" = now
  )
  
  # adds extra keywords from ...
  # removes keywords that are already in text_segment1 or that are named session or that have no name
  extra_keywords = setdiff(names(dots), c(names(text_segment1), "session", ""))
  # removes keywords that are already in text_segment2 ($PnN, $PnS, $PnB, $PnE, $PnR)
  extra_keywords = grep("^\\$P\\d.*[N|S|B|E|R|]$", extra_keywords, value = TRUE, invert = TRUE)
  # gets keywords in dots
  extra_keywords = dots[extra_keywords]
  # removes keywords whose values are NULL
  extra_keywords = extra_keywords[sapply(extra_keywords, length) != 0]
  # adds to text_segment1
  text_segment1 = c(text_segment1, extra_keywords)
  
  # determines length of data
  data_length = 4 * L * nrow(features)
  
  # determines length of mandatory param
  N = names(text_segment1)
  text1_length = sum(c(nchar("$BEGINDATA"), 2, # 2 for additional delimiters
                       nchar("$ENDDATA"),   2, # 2 for additional delimiters, there is already one at the beg of text2
                       sapply(1:length(text_segment1), FUN = function(i) {
                         length(charToRaw(paste(c("",
                                                  gsub(delimiter, delimiter_esc, N[i], fixed=TRUE), # FIXME should we escape keyword ?
                                                  gsub(delimiter, delimiter_esc, text_segment1[[i]], fixed=TRUE)), # delimiter if present in value should be escaped
                                                collapse = delimiter)))
                       }), text2_length,    1, # 1 for additional delimiters, to terminate text2
                                           57  # 57 for end of header
  ))
  
  # compute missing offsets 
  # ENDSTEXT
  # determining text_end is tricky since it depends on it own length
  # so we use a function to optimize it
  f = function(x, text_length) {
    data_beg = x + 1
    if(data_beg >= 1e9) data_beg = 0
    data_end = x + data_length - 1
    if(data_end >= 1e9) data_end = 0
    ans = text_length + nchar(data_beg) + nchar(data_end) 
    if(ans != x) ans = f(x = ans, text_length = text_length)
    return(ans)
  }
  text_end = f(x = text1_length, text_length = text1_length)
  if(text_end >= 1e9) stop("primary text segment is too big")
  header$text_end = sprintf("%8i", text_end)
  
  # BEGINDATA
  data_beg = text_end + 1 # +1 because data start just after text segment end
  if(data_beg >= 1e9) {
    header$data_beg = sprintf("%8i", 0)
  } else {
    header$data_beg = sprintf("%8i", data_beg)
  }
  text_segment1 = c(text_segment1, list("$BEGINDATA" = as.character(data_beg)))
  
  # ENDDATA
  data_end = data_beg + data_length - 1 #-1 because last data byte is at minus one from total length
  if(data_end >= 1e9) {
    header$data_end = sprintf("%8i", 0)
  } else {
    header$data_end = sprintf("%8i", data_end)
  }
  text_segment1 = c(text_segment1, list("$ENDDATA" = as.character(data_end)))
  
  towrite = file(description = file_w, open = "wb")
  tryCatch({
    # writes header
    writeBin(object = charToRaw(paste0(header, collapse="")), con = towrite)
    
    # writes text segment1
    N = names(text_segment1)
    lapply(1:length(text_segment1), FUN = function(i) {
      writeBin(object = charToRaw(paste(c("", 
                                          gsub(delimiter, delimiter_esc, N[i], fixed=TRUE), # FIXME should we escape keyword ?
                                          gsub(delimiter, delimiter_esc, text_segment1[i], fixed=TRUE)), # delimiter if present in value should be escaped
                                        collapse = delimiter)), con = towrite)
    })
    
    # writes text segment2
    lapply(1:length(text_segment2), FUN = function(i) {
      writeBin(object = text_segment2[[i]], con = towrite)
    })
    writeBin(object = charToRaw(delimiter), con = towrite) # we add final delimiter after the last keyword-value
    
    # export features values in data segment
    apply(features, 1, FUN = function(x) writeBin(x, con = towrite, size = 4))
    
  }, error = function(e) {
    stop(paste0("Can't create 'write_to' file.\n", write_to,
                ifelse(overwritten,"\nFile was not modified.\n","\n"),
                "See pre-file @\n", normalizePath(file_w, winslash = "/"), "\n",
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
