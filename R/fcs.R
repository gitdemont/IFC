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
    vals = foo[-(seq_len(feat_l+1))]
    if(length(vals) != feat_l^2) stop("'spillover' keyword does not fulfill fcs specifications")
    return(matrix(as.numeric(vals), ncol=feat_l, nrow=feat_l, dimnames=list(NULL, feat_n)))
  }
}

#' @title FCS Keyword Checker
#' @description
#' Helper to check that FCS keyword-value pairs are compliant with specifications
#' @param text a named list of keywords values.
#' @param delimiter delimiter used to separate keyword-value pairs. /!\ NOTE that files with 0x00 'delimiter' can _NOT_ be parsed.
#' @param version version to check keywords compliance against. Default is 3.0.
#' @param encoding name of the encoding for raw to character conversion. Default is "UTF-8".
#' @param fun function to execute when mandatory parameters are not met. Default is "warning". Allowed are "stop","warning","message","return".
#' @param ... other arguments to be passed.
#' @keywords internal
FCS_check_keywords <- function(text, delimiter, version = 3.0, encoding = "UTF-8", fun = "warning", ...) {
  # prepare return message
  msg = c()
  
  # check inputs
  assert(text, typ = "list")
  assert(fun, len = 1, alw = c("stop","warning","message","return"))
  version = suppressWarnings(as.numeric(version)); version = na.omit(version); assert(version, len=1)
  encoding = na.omit(encoding); assert(encoding, len=1, typ="character")
  assert(delimiter, len=1, typ="character")
  raw_delimiter = attr(delimiter, "raw")
  if((length(raw_delimiter) != 1) || (typeof(raw_delimiter) != "raw")) stop("can't find valid \"raw\" attribute of 'delimiter'")
  
  # modify encoding FCS keywords values should be UTF-8 and keywords themselves ASCII (which are part of UTF-8)
  old_enc <- options("encoding")
  on.exit(options(old_enc))
  options("encoding" = encoding)
  
  # set keywords to upper case
  names(text) = toupper(names(text))
  N = names(text)
  
  # check required keywords
  key_mandatory = c("$DATATYPE","$PAR","$MODE","$BYTEORD","$NEXTDATA",
                    "$BEGINSTEXT", "$ENDSTEXT",
                    "$BEGINANALYSIS", "$ENDANALYSIS",
                    "$BEGINDATA", "$ENDDATA", 
                    "$CYT", "$TOT")
  if(version == 2.0) key_mandatory = key_mandatory[1:5]
  if(version <  3.1) key_mandatory = setdiff(key_mandatory, c("$CYT"))
  if(version >= 3.2) key_mandatory = setdiff(key_mandatory, c("$MODE", "$BEGINSTEXT", "$ENDSTEXT", "$BEGINANALYSIS", "$ENDANALYSIS"))
  tmp = key_mandatory %in% N
  if(!all(tmp)) msg = c(msg, paste0("`REQUIRED not found`:\n\t- ", paste0(key_mandatory[!tmp], collapse = "\n\t- ")))
  
  # check type/mode/byteord
  type = text[["$DATATYPE"]]
  byteord = text[["$BYTEORD"]]
  if(!(type %in% c("A","I","F","D"))) msg = c(msg, paste0("`non-compatible $DATATYPE[", type,"] (allowed are \"A\",\"I\",\"F\",\"D\")`"))
  if(version >= 3.1) if(type %in% c("A"))  msg = c(msg, paste0("`deprecated $DATATYPE[", type,"]`"))
  mode = text[["$MODE"]]
  if(version <= 3.1) {
    if(!(mode %in% c("L","C","U"))) msg = c(msg, paste0("`non-compatible $MODE[",mode,"] (allowed are \"L\",\"C\",\"U\")`"))
    if(mode %in% c("C","U"))  msg = c(msg, paste0("`deprecated $MODE[", mode,"]`"))
  } else {
    if("$MODE" %in% N) {
      msg = c(msg, paste0("`$MODE is a deprecated keyword`"))
      if(!(mode %in% c("L"))) msg = c(msg, paste0("`non-compatible $MODE[",mode,"] (allowed is \"L\")`"))
    }
  }
  if(version >= 3.1) if(!(byteord %in% c("1,2,3,4","4,3,2,1"))) msg = c(msg, paste0("`non-compatible $BYTEORD[",byteord,"] (allowed are \"1,2,3,4\",\"4,3,2,1\")`"))
  
  # check number of objects/parameters
  n_obj = na.omit(suppressWarnings(as.integer(text[["$TOT"]])))
  n_par = na.omit(suppressWarnings(as.integer(text[["$PAR"]])))
  if((length(n_obj) == 1) && (n_obj == 0)) msg = c(msg, "`$TOT is 0`")
  if((length(n_par) == 1) && (n_par == 0)) msg = c(msg, "`$PAR is 0`")
  if(version > 2.0) {
    foo = grep("^\\$P\\d+N$", N, value = TRUE, ignore.case = TRUE)
    if(length(foo) != n_par) msg = c(msg, paste0("`$PnN mismatch between found[",length(foo),"] vs expected[",n_par,"]`"))
  }
  
  # check uniqueness
  tmp = duplicated(N)
  if(any(tmp)) msg = c(msg, paste0("`non unique`:\n\t- ", paste0(N[tmp], collapse = "\n\t- ")))
  
  # check uniqueness of PnN
  PnN = text[paste0("$P",seq_len(n_par),"N")]
  foo = PnN[sapply(PnN, length) != 0]
  tmp = duplicated(foo)
  if(any(tmp)) {
    bar = sapply(unique(foo[tmp]), FUN = function(x) paste0(names(which(x == foo)), collapse = ","))
    if(any(tmp)) msg = c(msg, paste0("`non unique $PnN`:\n\t- ", paste0(paste0(bar, "[",unique(foo[tmp]), "]"), collapse = "\n\t- ")))
  }
  
  # check non empty
  tmp = sapply(text, length) == 0
  if(any(tmp)) msg = c(msg, paste0("`empty value",ifelse(sum(tmp)==1,"","s")," found`:\n\t- ", paste0(N[tmp], collapse = "\n\t- ")))
  
  # check keywords names use only 0x20 - 0xFE characters
  if(version >= 3.0) {
    # starting ver 3.0 spe says: The TEXT part should not contain return (ASCII 13), line feed (ASCII 10) or other unpritable characters (unless they are value or delimiter)
    tmp = sapply(N, FUN = function(x) {v = charToRaw(x); all(v >= 0x20 & v <= 0xFE) })
    if(!all(tmp)) msg = c(msg, paste0("`bad BYTE in name allowed ASCII are [0x20-0xFE (32-126)]`:\n\t- ", paste0(N[!tmp], collapse = "\n\t- ")))
  } else {
    # for ver 2.0 spe mentions : The TEXT part should not contain 'carriage return' or 'line feed' characters  (unless they are value or delimiter)
    tmp = sapply(N, FUN = function(x) {v = charToRaw(x); all(v != 0x0A & v != 0x0D) })
    if(!all(tmp)) msg = c(msg, paste0("`bad BYTE in name 0x0A and 0x0D are not allowed`:\n\t- ", paste0(N[!tmp], collapse = "\n\t- ")))
  }
  
  # check numeric keywords are not padded with characters other than 0
  bar = TRUE
  old_loc <- setloc(c("LC_NUMERIC" = "C"))
  tryCatch({
    foo = sapply(text, FUN = function(x) (length(x) !=0) && (!is.na(suppressWarnings(as.numeric(x)))))
    bar = sapply(text[foo], FUN = function(x) {
      xx = tolower(x)
      if(xx %in% c("true","false")) return(TRUE) # to handle TRUE/FALSE to num conversion
      xx = strsplit(x = xx, split = "", fixed = TRUE)[[1]]
      if(any(xx == "x")) return(TRUE) # to handle raw to num conversion
      all(xx %in% c("0","1","2","3", "4", "5", "6", "7", "8", "9", "+", "-",".","e"))
    })
  }, finally = suppressWarnings(setloc(old_loc)))
  if(!all(bar)) msg = c(msg, paste0("`padded numeric (only padding with 0 is allowed)`:\n\t- ", paste0(paste0(N[foo][!bar], "[",text[foo][!bar], "]"), collapse = "\n\t- ")))
  
  # check keyword-value pairs do not start with delimiter
  # delimiter = rawToChar(delimiter)
  foo = sapply(N, FUN = function(x) {
    if(length(x) == 0) return(FALSE)
    substr(x, 1,1) == delimiter
  })
  bar = sapply(text, FUN = function(x) {
    if(length(x) == 0) return(FALSE)
    substr(x, 1,1) == delimiter
  })
  tmp = foo | bar
  if(any(tmp)) msg = c(msg, paste0("`start with delimiter`:\n\t- ", paste0(paste0(N[tmp], "[",text[tmp], "]"), collapse = "\n\t- ")))
  
  # check range and amplification
  foo = sapply(seq_len(n_par), FUN = function(i) {
    msg = c()
    PnB = paste0("$P",i,"B")
    PnE = paste0("$P",i,"E")
    PnG = paste0("$P",i,"G")
    PnN = paste0("$P",i,"N")
    PnR = paste0("$P",i,"R")
    PnT = paste0("$P",i,"DATATYPE")
    TYPE = type
    if((version >= 3.2) && (PnT %in% N)) {
      TYPE = text[[PnT]]
      if(length(TYPE) == 0) TYPE = ""
      if(!(TYPE %in% c("I","F","D"))) msg = c(msg, paste0("invalid ", PnT , "[", TYPE,"] (allowed are \"I\",\"F\",\"D\")"))
    }
    if(length(text[[PnB]]) == 0) {
      msg = c(msg, paste0(PnB, " not found (REQUIRED)"))
    } else {
      bit = na.omit(suppressWarnings(as.integer(text[[PnB]])))
      if(length(bit) == 0) {
        if(!(identical(text[[PnB]],"*") && identical(type, "A"))) msg = c(msg, paste0("invalid ", PnB, "[",text[[PnB]],"]"))
      } else {
        if(version >= 3.2 && !identical(type, "A") && as.logical(bit %% 8)) msg = c(msg, paste0("invalid ", PnB, "[", text[[PnB]],"] not divisible by 8 are deprecated in FCS >= 3.2"))
        if(TYPE == "F" && bit != 32) if(version >= 3.1) msg = c(msg, paste0("invalid ", PnB, " should be \"32\" with $DATATYPE[F]"))
        if(TYPE == "D" && bit != 64) if(version >= 3.1) msg = c(msg, paste0("invalid ", PnB, " should be \"64\" with $DATATYPE[D]"))
      }
    }
    if(length(text[[PnR]]) == 0) {
      msg = c(msg, paste0(PnR, " not found (REQUIRED)"))
    } else {
      ran = na.omit(suppressWarnings(as.numeric(text[[PnR]])))
      if(length(ran) == 0) msg = c(msg, paste0(PnR,"[",text[[PnR]],"] is not valid"))
    }
    if(length(text[[PnN]]) == 0) {
      if(version > 2.0) msg = c(msg, paste0(PnN, " not found (REQUIRED)")) # $PnN is not listed as REQUIRED in FCS2.0
    } else {
      if(grepl(",", text[[PnN]], fixed = TRUE)) c(msg, paste0(PnN,"[",text[[PnN]],"] should not contain ','")) 
    }
    if(length(text[[PnE]]) == 0) {
      if(version > 2.0) msg = c(msg, paste0(PnE, " not found (REQUIRED)")) # $PnE is not listed as REQUIRED in FCS2.0
    } else {
      trans = na.omit(suppressWarnings(as.numeric(strsplit(text[[PnE]], split = ",", fixed = TRUE)[[1]])))
      if(length(trans) != 2) {
        msg = c(msg, paste0(PnE,"[",text[[PnE]],"] should be \"f1,f2\" with numeric f1 and f2"))
      } else {
        if((version >= 3.1) && (TYPE == "D" || TYPE == "F")) { # starting from FCS3.1 $PnE has to be "0,0" for DATATYPE D and F
          if((trans[1] != 0) || (trans[2] != 0)) msg = c(msg, paste0(PnE,"[",text[[PnE]],"] should be \"f1,f2\" where f1=0 AND f2=0 with $DATATYPE[",TYPE,"]"))
        } else { # otherwise $PnE always (? for 1.0) had to be f1,f2 with  f1=f2=0 OR (f1!=0 AND f2!=0)
          if(!(all(trans == 0) || all(trans != 0))) msg = c(msg, paste0(PnE,"[",text[[PnE]],"] should be \"f1,f2\" where (f1=0 AND f2=0, for linear) OR (f1!=0 AND f2!=0, for log) with $DATATYPE[",TYPE,"]"))
        }
      }
    }
    gain = na.omit(suppressWarnings(as.numeric(text[[PnG]])))
    if(length(gain) != 0) {
      if((version >= 3.2) && !identical(TYPE, "I")) msg = c(msg, paste0(PnG,"[",text[[PnG]],"] shall be used for integer type only not $DATATYPE[",TYPE,"]"))
      if(length(text[[PnE]]) != 0) {
        trans = na.omit(suppressWarnings(as.numeric(strsplit(text[[PnE]], split = ",", fixed = TRUE)[[1]])))
        if((length(trans) == 2) && (gain != 1) && ((trans[1] != 0) || (trans[2] != 0))) msg = c(msg, paste0(PnG,"[",text[[PnG]],"] should not be used with logarithmic amplification ($PnE != \"0,0\")",PnE,"[",text[[PnE]],"]"))
      }
    }
    return(paste0(msg, collapse="|"))
  })
  bar = unlist(recursive = FALSE, use.names = FALSE, lapply(foo, FUN = function(x) x != ""))
  if(any(bar)) {
    if(version >= 3.2) {
      msg = c(msg, paste0("`bad PnB/PnDATATYPE/PnE/PnG/PnN/PnR keywords`:\n\t- ", paste0(foo[bar], collapse="\n\t- "))) 
    } else {
      msg = c(msg, paste0("`bad PnB/PnE/PnG/PnN/PnR keywords`:\n\t- ", paste0(foo[bar], collapse="\n\t- ")))
    }
  }
  
  # check spillover
  if("$SPILLOVER" %in% N) {
    tryCatch({
      sp = convert_spillover(text[["$SPILLOVER"]])
      if(!all(colnames(sp) %in% unlist(recursive = FALSE, use.names = FALSE, PnN))) stop("'spillover' is defined with non PnN names [", paste0(colnames(sp),collapse = ","),"]")
    }, error = function(e) {
      msg <<- c(msg, paste0("$SPILLOVER ", e$message))
    })
  }
  if(length(msg) != 0) msg = c(sprintf("non FCS%.1f compliant keywords", version), msg)
  msg = gsub("\r","",msg,fixed=TRUE)
  if(fun == "return") return(msg)
  if(length(msg) != 0) {
    args = list(paste0(msg, collapse = "\n"))
    if(fun == "warning") args = c(args, list(call. = FALSE, immediate. = TRUE))
    do.call(what = fun, args = args)
  }
}

#' @title FCS Header Parser
#' @description
#' Helper to parse header segment from Flow Cytometry Standard (FCS) compliant files.
#' @param fileName path to file.
#' @param header, a list whose members define the "at" offset from header$start$at and the "n" number of bytes to extract:\cr
#' - start: where to start reading current FCS dataset.       Default is list(at = 0,  n = 6),\cr
#' - space: where to retrieve space.                          Default is list(at = 6,  n = 4),\cr
#' - text_beg: where to retrieve file text segment beginning. Default is list(at = 10, n = 8),\cr
#' - text_end: where to retrieve file text segment end.       Default is list(at = 18, n = 8),\cr
#' - data_beg: where to retrieve file data segment beginning. Default is list(at = 26, n = 8),\cr
#' - data_end: where to retrieve file data segment end.       Default is list(at = 34, n = 8).
#' @param encoding name of the encoding for raw to character conversion. Default is "UTF-8".
#' @param ... other arguments to be passed.
#' @keywords internal
readFCSheader <- function(fileName, header, encoding = "UTF-8", ...) {
  # prepare fileName
  if(missing(fileName)) stop("'fileName' can't be missing")
  assert(fileName, len = 1)
  fileName = enc2native(normalizePath(fileName, winslash = "/", mustWork = FALSE))
  if(!file.exists(fileName)) stop(paste("can't find", fileName, sep=" "))
  fsize = file.size(fileName)
  encoding = na.omit(encoding); assert(encoding, len=1, typ="character")
  
  opt_default = eval(formals(readFCS)$options)
  if(missing(header)) {
    header = opt_default$header
  } else {
    for(i in c("start", "space", "text_beg", "text_end", "data_beg", "data_end")) {
      if(!(i %in% names(header))) header[[i]] <- opt_default$header[[i]]
    }
  }
  
  # ensure start$at is valid
  at = suppressWarnings(as.numeric(header[["start"]]$at[1]))
  at = na.omit(at[at >=0])
  if((length(at) == 0) || (at > file.size(fileName))) stop("HEADER segment: start$at[",at,"] points to outside of the file")
  assert(at, len=1)
  
  # create connection binary reading
  toread = file(description = fileName, open = "rb")
  on.exit(close(toread))
  # modify encoding FCS keywords values should be UTF-8 and keywords themselves ASCII (which are part of UTF-8)
  old_enc <- options("encoding")
  options("encoding" = encoding) # FIXME should we read $UNICODE before parsing for FCS < 3.1
  tryCatch({
    # we will read offsets from options
    # FIXME, should we validate each header entry ?
    # e.g. header$start == FCSx.x
    # and so on ...
    header = sapply(names(header), simplify = FALSE, FUN = function(x) {
      goto = header[[x]]$at[1] + ifelse(x == "start", 0, at)
      if(goto > fsize) stop("HEADER segment: ",x,"[",goto,"] points to outside of the file")
      seek(toread, header[[x]]$at[1] + ifelse(x == "start", 0, at))
      if(x %in% c("start", "space")) {
        raw = rawToChar(readBin(toread, what = "raw", n = header[[x]]$n))
      } else {
        raw = trimws(x = rawToChar(readBin(toread, what = "raw", n = header[[x]]$n)))
        raw = suppressWarnings(na.omit(as.integer(raw) + at))
      }
      raw
    })
  }, finally = options(old_enc))
  return(structure(header, "offset" = at))
}

#' @title FCS Delimiter Reader
#' @description
#' Helper to extract delimiter from Flow Cytometry Standard (FCS) compliant files.
#' @param fileName path to file.
#' @param at offset of delimiter. Default is 58.
#' @param version version to check keywords compliance against. Default is 3.0.
#' @param encoding name of the encoding for raw to character conversion. Default is "UTF-8".
#' @param ... other arguments to be passed.
#' @keywords internal
readFCSdelimiter <- function(fileName, at = 58, version = 3.0, encoding = "UTF-8", ...) {
  # prepare fileName
  if(missing(fileName)) stop("'fileName' can't be missing")
  assert(fileName, len = 1)
  fileName = enc2native(normalizePath(fileName, winslash = "/", mustWork = FALSE))
  if(!file.exists(fileName)) stop(paste("can't find", fileName, sep=" "))
  encoding = na.omit(encoding); assert(encoding, len=1, typ="character")
  version = suppressWarnings(as.numeric(version)); version = na.omit(version); assert(version, len=1)
  
  # ensure start$at is valid
  at = suppressWarnings(as.numeric(at))
  at = na.omit(at[at >=0])
  if((length(at) == 0) || (at > file.size(fileName))) stop("DELIMITER segment: at[",at,"] points to outside of the file")
  assert(at, len=1)
  
  # create connection binary reading
  toread = file(description = fileName, open = "rb")
  on.exit(close(toread))
  # modify encoding FCS keywords values should be UTF-8 and keywords themselves ASCII (which are part of UTF-8)
  old_enc <- options("encoding")
  options("encoding" = encoding) # FIXME should we read $UNICODE before parsing for FCS < 3.1
  delimiter = character()
  raw_delimiter = raw()
  tryCatch({
    seek(toread, at)
    raw_delimiter = readBin(con = toread, what = "raw", n = 1)
    delimiter = rawToChar(raw_delimiter)
  }, finally = options(old_enc))
  
  if(length(raw_delimiter) == 1) {
    # for FCS2.0 delimiter is a BYTE
    # for FCS3.0 delimiter is an ASCII
    # starting from 3.1 delimiter is an ASCII that is not 0x00 nor 0xFF
    if(version == 3.0) if(raw_delimiter > 0x7F) warning("DELIMITER segment: should be any ASCII [0x00-0x7F (0-127)]",
                                                        call. = FALSE, immediate. = TRUE)
    if(version >= 3.1) if((raw_delimiter == 0x00) || (raw_delimiter >= 0x7E)) warning("DELIMITER segment: should be a [0x01-0x7E (1-126)] ASCII character",
                                                                                      call. = FALSE, immediate. = TRUE)
  }
  return(structure(delimiter, "raw" = raw_delimiter))
}

#' @title FCS Text Parser
#' @description
#' Helper to parse text segment from Flow Cytometry Standard (FCS) compliant files.
#' @param fileName path to file.
#' @param delimiter delimiter used to separate keyword-value pairs. /!\ NOTE that files with 0x00 'delimiter' can _NOT_ be parsed.
#' @param start offset of text start. Default is 0.
#' @param end offset of text end. Default is 0.
#' @param encoding name of the encoding for raw to character conversion. Default is "UTF-8".
#' @param empty whether to allow empty values when parsing text segment. Default is FALSE.
#' @param trim remove whitespace in keywords names. Default is "none". Allowed are "both", "left", "right" and "none".
#' @param ... other arguments to be passed.
#' @keywords internal
readFCStext <- function(fileName, delimiter, start = 0, end = 0, encoding = "UTF-8", empty = FALSE, trim = "none", ...) {
  if(missing(fileName)) stop("'fileName' can't be missing")
  assert(fileName, len=1)
  fileName = enc2native(normalizePath(fileName, winslash = "/", mustWork = FALSE))
  if(!file.exists(fileName)) stop(paste("can't find", fileName, sep=" "))
  assert(delimiter, len=1, typ="character")
  raw_delimiter = attr(delimiter, "raw")
  if((length(raw_delimiter) != 1) || (typeof(raw_delimiter) != "raw")) stop("can't find valid \"raw\" attribute of 'delimiter'")
  start = na.omit(as.numeric(start))
  end = na.omit(as.numeric(end))
  if((length(start) != 1) || (length(end) != 1) || (end <= start)) stop("bad TEXT segment offsets")
  encoding = na.omit(encoding); assert(encoding, len=1, typ="character")
  assert(empty, len=1, alw = c(TRUE, FALSE))
  assert(trim, len=1, alw = c("none","left","right","both"))
  if(raw_delimiter == raw(1)) stop("delimiter 0x00 is not supported")
  
  # modify encoding FCS keywords values should be UTF-8 and keywords themselves ASCII (which are part of UTF-8)
  old_enc <- options("encoding")
  on.exit(options(old_enc))
  options("encoding" = encoding)
  
  # create connection binary reading
  toread = file(description = fileName, open = "rb")
  on.exit(close(toread), add = TRUE)
  last = raw()
  nend = end - 2
  while(!(raw_delimiter %in% last) && (abs(nend - end) <= 2)) {
    nend = nend + 1
    seek(toread, nend)
    last = readBin(toread, what = "raw", n = 1)
  }
  if(nend != end) {
    if(abs(nend - end) == 1) {
      warning("TEXT (or sup. TEXT) segment: offset is off by ", nend - end, " byte",
              call. = FALSE, immediate. = TRUE)
    } else {
      warning("TEXT (or sup. TEXT) segment: can't find final delimiter",
              call. = FALSE, immediate. = TRUE)
    }
  }
  
  seek(toread, start)
  text = rawToChar(readBin(toread, what = "raw", n = nend - start))
  # when same character as delimiter is used within keyword-value pair it has to be escaped (repeated twice)
  # according to FCS spe 2.0:
  # -if the separator appears in a keyword or in a keyword value, it must be "quoted" by being repeated
  # -since null (zero length) keywords or keyword values are not permitted, two consecutive separators can never occur between a value and a keyword
  # we generate a 20 random characters delim_esc that does not contain delimiter
  # we also ensure that this delim is not found elsewhere in the TEXT segment
  found = 1
  while(found) {
    # back compatible with old R version, no need for accuracy since it is just for finding a non existing string that allow parsing
    delim_esc = gen_altnames("foo", random_seed = list(seed=found,"Mersenne-Twister", "Inversion", "Rounding"))
    delim_esc = strsplit(x = delim_esc, split = delimiter, fixed = TRUE)[[1]]
    delim_esc = delim_esc[delim_esc!=""]
    delim_esc = paste0(delim_esc, collapse="")
    found = cpp_scanFirst(fileName, charToRaw(delim_esc), start = start, end = end)
  }
  # we 1st look at double delimiter instance and substitute it with delim_esc
  text = gsub(pattern = paste0(delimiter,delimiter), replacement = delim_esc, x = text, fixed = TRUE)
  # then text is split with delimiter
  text = strsplit(x = text, split = delimiter, fixed = TRUE)[[1]]
  # then escaped double delimiter is replaced with only one delimiter
  text = gsub(pattern = delim_esc, replacement = delimiter, x = text, fixed = TRUE)
  # remove 1st empty value (this happen when 1st keyword starts with delimiter)
  is_1st_empty = FALSE
  while((length(text) >= 1) && (text[1] == "")) {
    if(!is_1st_empty) warning("TEXT (or sup. TEXT) segment: 1st keyword starts with delimiter",
                              call. = FALSE, immediate. = TRUE)
    is_1st_empty = TRUE
    text = text[-1]
  }
  # finally keyword-value pairs are converted to named list
  id_val = seq(from = 2, to = length(text), by = 2)
  id_key = id_val-1
  text = structure(as.list(text[id_val]), names = text[id_key])
  # try to fix empty values
  # it happens, See Bras A.E. and van der Velden V.H.J. at \doi{10.1002/cyto.a.24187}, that fcs writers produce non compliant files
  # with empty values so we try to detect them and fix them with the reasoning that delimiter never occur in keys (=keywords names)
  found = grepl(delimiter, names(text), fixed = TRUE)
  if(any(found)) {
    if(empty) {
      text_ok = text[!found]
      text_nok = text[found]
      text_nok = strsplit(names(text_nok), delimiter, fixed = TRUE)
      text_nok = unlist(recursive = FALSE, use.names = TRUE,
                        lapply(seq_along(text[found]), FUN = function(ii) {
                          L = length(text_nok[[ii]])
                          structure(lapply(seq_along(text_nok[[ii]]), FUN = function(kk) {
                            if(kk == L) return(text[found][[ii]])
                            character()
                          }), names = text_nok[[ii]])
                        }))
      text = c(text_ok, text_nok)
    } else {
      warning("TEXT (or sup. TEXT) segment: found delimiter in keywords, please consider using `options$text_empty` = TRUE",
              call. = FALSE, immediate. = TRUE)
    }
  }
  msg = paste0("TEXT (or sup. TEXT) segment: found standard keywords padded with whitespace(s), please consider using `options$text_trim` != \"",trim,"\"\n\t-")
  names(text) = switch(trim,
                       "both" = { 
                         trimws(names(text), which = "both")
                       },
                       "left" = { 
                         found = grep("^\\$[[:alnum:]]+[[:space:]]+$", names(text), ignore.case = TRUE, value = TRUE)
                         if(length(found) != 0) warning(paste0(msg, paste0(found, collapse = "\n\t-")), call. = FALSE, immediate. = TRUE)
                         trimws(names(text), which = "left")
                       },
                       "right" = { 
                         found = grep("^[[:space:]]+\\$[[:alnum:]]$", names(text), ignore.case = TRUE, value = TRUE)
                         if(length(found) != 0) warning(paste0(msg, paste0(found, collapse = "\n\t-")), call. = FALSE, immediate. = TRUE)
                         trimws(names(text), which = "right")
                       },
                       {
                         found = c(grep("^\\$[[:alnum:]]+[[:space:]]+$", names(text), ignore.case = TRUE, value = TRUE),
                                   grep("^[[:space:]]+\\$[[:alnum:]]+$", names(text), ignore.case = TRUE, value = TRUE))
                         if(length(found) != 0) warning(paste0(msg, paste0(found, collapse = "\n\t-")), call. = FALSE, immediate. = TRUE)
                         names(text)
                       })
  return(text)
}

#' @title FCS Data Parser
#' @description
#' Helper to parse data segment from Flow Cytometry Standard (FCS) compliant files.
#' @param fileName path to file.
#' @param text a named list of keywords values.
#' @param version version of FCS specification.
#' @param start offset of text start. Default is 0.
#' @param end offset of text end. Default is 0.
#' @param scale whether to apply data scaling. It only applies when fcs file is stored as DATATYPE "I". Default is TRUE.\cr
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other arguments to be passed.
#' @keywords internal
readFCSdata <- function(fileName, text, version, start = 0, end = 0, scale = TRUE, display_progress = TRUE, ...) {
  if(missing(fileName)) stop("'fileName' can't be missing")
  assert(fileName, len=1)
  fileName = enc2native(normalizePath(fileName, winslash = "/", mustWork = FALSE))
  if(!file.exists(fileName)) stop(paste("can't find", fileName, sep=" "))
  assert(text, typ = "list")
  start = na.omit(as.numeric(start))
  end = na.omit(as.numeric(end))
  if((length(start) != 1) || (length(end) != 1) || (end <= start)) stop("DATA segment: bad offsets")
  assert(scale, len=1, alw = c(TRUE,FALSE))
  title_progress = basename(fileName)
  assert(display_progress, len = 1, alw = c(TRUE,FALSE))
  
  # force local
  locale_back <- setloc(c("LC_ALL" = "en_US.UTF-8"))
  enc_back <- options("encoding" = "UTF-8")
  on.exit(suspendInterrupts({setloc(locale_back); options(enc_back)}), add = TRUE)
  setloc(c("LC_NUMERIC" = "C"))
  
  # retrieve info to extract data
  type = text[["$DATATYPE"]]
  if((length(type) != 1) || !(type %in% c("A","I","F","D"))) stop("DATA segment: non-compatible $DATATYPE[",type,"]")
  n_obj = na.omit(suppressWarnings(as.integer(text[["$TOT"]])))
  n_par = na.omit(suppressWarnings(as.integer(text[["$PAR"]])))
  if(length(n_obj) == 0) stop("DATA segment: $TOT can't be found")
  if(length(n_par) == 0) stop("DATA segment: $PAR can't be found")
  features_names = grep("^\\$P\\d+N$", names(text), value = TRUE, ignore.case = TRUE)
  
  # hereafter we create several bit_* variables
  # bit_t : PnDATATYPE, type of the value FCS >= 3.2
  # bit_v : PnB, bits depth of the value
  # bit_r : PnR, bits range of the value
  # bit_n : number of bytes to read 
  # bit_d : bits depth to read ( = 8 * bit_n )
  # bit_o : bytes order to read
  # bit_m : bits mask, for instance if bit_v is 10 bits but the value is read from 16 bits then 6 bits are not used
  # FIXME it is not clear how to perform bit packing 
  bit_t = sapply(paste0("$P",seq_len(n_par),"DATATYPE"), USE.NAMES = FALSE, simplify = TRUE, FUN = function(x) {
    ans = text[[x]]
    if(length(ans) == 0) ans = type
    return(ans)
  })
  bad_b = c()
  bit_v = unlist(recursive = FALSE, use.names = FALSE, lapply(seq_len(n_par), FUN = function(i) {
    foo = suppressWarnings(as.integer(text[[paste0("$P",i,"B")]]))
    ans = foo
    if(identical(bit_t[i], "F")) ans = 32L
    if(identical(bit_t[i], "D")) ans = 64L
    if(!identical(ans, foo)) bad_b <<- c(bad_b,  i)
    if(length(ans) == 0) ans = NA_integer_
    return(ans)
  }))
  if(length(bad_b) != 0) {
    warning("DATA segment: $PnB keyword",ifelse(length(bad_b) == 1, " has", "s have")," been forced\n",
            paste0(paste0("\t- ", "$P",bad_b,"B"), ": ", bit_v[bad_b], collapse = "\n"),
            call. = FALSE, immediate. = TRUE)
  }
  
  msg = c()
  if(n_par != length(features_names)) msg = c(msg, paste0("DATA segment: mismatch between found[",length(features_names),"] vs expected[",n_par,"] number of $PnN"))
  if(n_par != length(na.omit(bit_v))) msg = c(msg, paste0("DATA segment: mismatch between found[",length(bit_v),"] vs expected[",n_par,"] number of $PnB"))
  if(length(msg) != 0) stop(paste0(msg, collapse = "\n"))
  data_bytes = end - start + 1
  
  if((version <= 3.1) || (type == "A")) {
    if((version != 0) && (length(unique(bit_t)) != 1)) stop("DATA segment: file with various $PnDATATYPE is not allowed")
    # create connection binary reading
    toread = file(description = fileName, open = "rb")
    on.exit(close(toread))
    # go to data start position
    seek(toread, start)
    # type "A" is deprecated in newer version of FCS specifications
    if(type == "A") {
      bit_t = rep("A", n_par)
      if((data_bytes + start) > file.size(fileName)) stop("DATA segment: points to outside of the file")
      if(length(unique(bit_v)) == 1) {
        if(text[["$P1B"]] == "*") {
          data = setdiff(strsplit(x = readBin(toread, what = "character", n = data_bytes),
                                  split = paste0(sapply(as.raw(c("0x20","0x09","0x2C","0x0D","0x0A")), rawToChar), collapse = "|"))[[1]],
                         "")
        } else {
          if(is.na(bit_v[1])) stop("DATA segment: bad $PnB definition for $DATATYPE[A]")
          data = gsub(paste0("(.{",bit_v,"})"), "\\1 ", readBin(toread, what = "character", n = data_bytes))
        }
      } else {
        raw = readBin(toread, what = "raw", n = data_bytes)
        if(display_progress) {
          pb = newPB(min = 0, max = n_par, initial = 0, style = 3)
          on.exit(endPB(pb), add = TRUE)
          data = sapply(seq_len(n_par), FUN = function(i_par) {
            setPB(pb, value = i_par, title = title_progress, label = "$DATATYPE[A]: extracting values")
            bits = suppressWarnings(as.integer(text[[paste0("$P",i_par,"B")]])) # each PnB determines number of bytes to extract
            # FIXME it is not clear how to deal with a mix of PnB == * and PnB == integer ?
            if(is.na(bits)) stop("DATA segment: bad $P",i_par,"B definition for $DATATYPE[A] ", text[[paste0("$P",i_par,"B")]])
            off = (i_par - 1) * n_par
            if((off + n_obj) > data_bytes) stop("DATA segment: buffer overrun")
            sapply(seq_len(n_obj), FUN = function(i_obj) {
              as.numeric(readBin(con = raw[i_obj + off], what = "character", n = bits))
            })
          })
        } else {
          data = sapply(seq_len(n_par), FUN = function(i_par) {
            bits = suppressWarnings(as.integer(text[[paste0("$P",i_par,"B")]])) # each PnB determines number of bytes to extract
            off = (i_par - 1) * n_par
            if((off + n_obj) > data_bytes) stop("DATA segment: buffer overrun")
            sapply(seq_len(n_obj), FUN = function(i_obj) {
              as.numeric(readBin(con = raw[i_obj + off], what = "character", n = bits))
            })
          })
        }
      }
    } else {
      # some files register wrong dataend offset resulting in an off-by-one byte
      # the following should allow to correct it
      if((data_bytes %% 8) %% 2) data_bytes = data_bytes + 1 # go +1
      if((data_bytes %% 8) %% 2) data_bytes = data_bytes - 1 # go  0
      if((data_bytes %% 8) %% 2) data_bytes = data_bytes - 1 # go -1
      off_by = (end - start + 1) - data_bytes
      if(off_by != 0) warning("DATA segment: offset is off by ",off_by," byte",
                              call. = FALSE, immediate. = TRUE)
      if((data_bytes %% 8) %% 2) stop("DATA segment: number of bytes does not respect fcs specifications")
      
      # extract order, endianness
      b_ord = text[["$BYTEORD"]]
      if(version >= 3.1) if(!any(identical(b_ord, "1,2,3,4"), identical(b_ord, "4,3,2,1"))) stop("DATA segment: bad endianness[",b_ord,"], FCS >= 3.1 spe. only allow `1,2,3,4` or `4,3,2,1`")
      bit_o = as.integer(strsplit(b_ord, split=",", fixed=TRUE)[[1]])
      b_ord = paste0(bit_o, collapse = ",")
      
      # determines endianness of the file
      endian = "unk"
      endian_l = paste0(seq_along(bit_o), collapse = ",")
      endian_b = paste0(rev(seq_along(bit_o)), collapse = ",")
      if(endian_l == b_ord) endian = "little"
      if(endian_b == b_ord) endian = "big"
      
      # register bit_v values
      bit_v_back = bit_v
      
      # try bit_v correction
      alw_b = 2^(3:6) # 8,16,32,64, sizes handled by R
      bit_corrected = FALSE
      if(any(!(bit_v %in% alw_b))) {
        bit_corrected = TRUE
        bit_v = sapply(bit_v, FUN = function(x) alw_b[x <= alw_b][1])
      }
      bit_n = unname(bit_v %/% 8)
      
      # try data_length correction and check accuracy of bit_v correction, if any
      if(sum(n_obj * bit_n) != data_bytes) { # data_length is not OK
        if(bit_corrected) {                  # bit_v was already corrected
          bit_v = bit_v_back
          bit_n = unname(bit_v %/% 8)
          bit_corrected = FALSE
          if(sum(n_obj*bit_n)!=data_bytes) { # bit_v is reverted but data_length still does not match
            stop("DATA segment: can't determine bit depth and DATA length")
          }
        } else {                             # bit_v was not corrected, we apply data_length correction
          data_bytes = sum(n_obj * bit_n)
          warning("DATA segment: number of bytes has been corrected",
                  call. = FALSE, immediate. = TRUE)
        }
      } else {                               # data_length is OK
        if(bit_corrected)             {      # check if bit_v was corrected to show a warning
          warning("DATA segment: bits depth have been corrected to next allowed value",
                  call. = FALSE, immediate. = TRUE)
        }
      }
      bit_d = unique(bit_v)
      
      args = sapply(seq_len(n_par), simplify = FALSE, USE.NAMES = TRUE, FUN = function(i) {
        if(identical(bit_t[i], "I")) return(list(what = "integer", size = bit_n[i], signed = bit_n[i] > 2))
        return(list(what = "numeric",  size = bit_n[i], signed = TRUE))
      })
      bit_r = sapply(seq_len(n_par), FUN = function(i) {
        x = text[[paste0("$P",i,"R")]]
        if(identical(bit_t[i], "I")) return(suppressWarnings(ceiling(log2(as.numeric(x))))) # as.integer("4294967296") results in NA so we use as.numeric
        return(c(32L, 64L)[identical(bit_t[i], "D") + 1L])
      })
      bit_m = lapply(seq_len(n_par), FUN = function(i_par) packBits(as.raw(sapply(seq_len(bit_v[i_par]), FUN = function(i) i <= min(bit_r[i_par],bit_v[i_par])))))
      
      # check data_length before starting reading
      if((data_bytes + start) > file.size(fileName)) stop("DATA segment: points to outside of the file")
      
      # fcs specifications mention that type "I" use unsigned integers only
      # but readBin can only extracts 8bits and 16bits unsigned integer. So, 
      # for 32bits and 64bits we have to extract signed integers and convert them afterwards
      # with cpp_v_intxx_to_uintxx functions
      if((length(unique(bit_t)) == 1) &&   # all same type
         (endian != "unk") &&              # endianness is either "little" or "big" and can be passed to readBin without reordering
         (length(bit_d) == 1) &&           # every channels have same bits depth
         (bit_d %in% c(8,16,32,64))) {     # bits depth is a size handled by R
        if(all(bit_r == bit_d)) {          # there is no masking to apply
          data = do.call(args = c(list(con = toread, 
                                       endian = endian, 
                                       n = data_bytes / args[[1]]$size),
                                  args[[1]]),
                         what = readBin)
        } else { # we are forced to apply bit masking if a PnR is not equal to bit_d
          # extract order for the whole data
          ord_ = cpp_get_bytes_order(n_obj, bit_n, bit_o, .Platform$endian != "little") # FIXME
          # extract mask for the whole data
          msk_ = rep(unlist(recursive = FALSE, use.names = FALSE, bit_m), n_obj)
          # extract the whole data in "raw"
          data = readBin(con = toread, what = "raw", n = data_bytes)
          # process data according to args, order and mask and make the conversion
          data = do.call(args = c(list(con = data[ord_] & msk_,
                                       endian = "little", # FIXME
                                       n = data_bytes / args[[1]]$size),
                                  args[[1]]),
                         what = readBin)
        }
        # convert to unsigned integers if needed
        if(type == "I") {
          if(args[[1]]$size == 4) data = cpp_v_int32_to_uint32(data)
          if(args[[1]]$size == 8) data = cpp_v_int64_to_uint64(data)
        }
        data = matrix(data, ncol = n_par, nrow = n_obj, byrow = TRUE)
      } else {
        # extract order for the whole data
        ord_ = cpp_get_bytes_order(n_obj, bit_n, bit_o, .Platform$endian != "little") # FIXME
        # extract mask for the whole data
        msk_ = rep(unlist(recursive = FALSE, use.names = FALSE, bit_m), n_obj)
        # define grouping parameter for the whole data
        spl_ = rep(unlist(recursive = FALSE, use.names = FALSE, mapply(FUN = rep, seq_along(bit_n), bit_n, SIMPLIFY = FALSE)), times = n_obj)
        
        if(display_progress) {
          lab = sprintf("$DATATYPE[%s]: extracting values", type)
          pb = newPB(min = 0, max = n_par, initial = 0, style = 3)
          on.exit(endPB(pb), add = TRUE)
        }
        
        # extract the whole data in "raw"
        data = readBin(con = toread, what = "raw", n = data_bytes)
        # process data according to args, order and mask and make the conversion
        data = unlist(recursive = FALSE, use.names = FALSE, by(list(v = data[ord_] & msk_, split = spl_), spl_, FUN = function(x) {
          i_par = x$split[1]
          if(display_progress) setPB(pb, value = i_par, title = title_progress, label = lab)
          if(args[[i_par]]$size %in% c(3,5,6,7)) { # for sizes not handled by R 
            M = c(1,2,4,8)[args[[i_par]]$size <= c(1,2,4,8)][1]
            m = M - args[[i_par]]$size
            bar = split(x$v,ceiling(seq_along(x$v)/args[[i_par]]$size))
            bar = lapply(bar, FUN = function(i) c(i, rep(as.raw(0x00), m)))
            foo = readBin(con = unlist(recursive = FALSE, use.names = FALSE, bar),
                          endian = "little", # FIXME ?
                          what = args[[i_par]]$what,
                          n = n_obj,
                          size = M,
                          signed = TRUE)
          } else { # size is handled by R
            foo = readBin(con = x$v,
                          endian = "little", # FIXME ?
                          what = args[[i_par]]$what,
                          n = n_obj,
                          size = args[[i_par]]$size,
                          signed = args[[i_par]]$signed)
          }
          # convert to unsigned integers if needed
          if(args[[i_par]]$what == "integer") {
            if(args[[i_par]]$size == 4) return(cpp_v_int32_to_uint32(foo))
            if(args[[i_par]]$size == 8) return(cpp_v_int64_to_uint64(foo))
          }
          foo
        }))
        data = matrix(data, ncol = n_par, nrow = n_obj, byrow = FALSE)
      }
    }
  } else {
    data = try(silent = TRUE, expr = {
      b_ord = text[["$BYTEORD"]]
      endian = "unk"
      if(b_ord == "1,2,3,4") endian = "little" 
      if(b_ord == "4,3,2,1") endian = "big" 
      if(endian == "unk") stop("DATA segment: bad endianness[",b_ord,"], FCS >= 3.1 spe. only allow `1,2,3,4` or `4,3,2,1`")
      
      # some files register wrong dataend offset resulting in an off-by-one byte
      # the following should allow to correct it
      if((data_bytes %% 8) %% 2) data_bytes = data_bytes + 1 # go +1
      if((data_bytes %% 8) %% 2) data_bytes = data_bytes - 1 # go  0
      if((data_bytes %% 8) %% 2) data_bytes = data_bytes - 1 # go -1
      off_by = (end - start + 1) - data_bytes
      if(off_by != 0) warning("DATA segment: offset is off by ",off_by," byte", call. = FALSE, immediate. = TRUE)
      if((data_bytes %% 8) %% 2) stop("DATA segment: number of bytes does not respect fcs specifications")
      
      # backup bit_v values
      bit_v_back = bit_v
      
      # try bit_v correction
      alw_b = as.integer(8*(1:8)) # divisible by 8
      if(any(!(bit_v %in% alw_b))) {
        bit_v = sapply(bit_v, FUN = function(x) alw_b[x <= alw_b][1])
      }
      bit_n = unname(bit_v %/% 8)
      
      types = c(0L,1L,2L,3L)[match(bit_t, c("A", "F", "D", "I"))] # A is not allowed
      sizes = c(1L,2L,3L,4L,5L,6L,7L,8L)[match(bit_n, c(1, 2, 3, 4, 5, 6, 7, 8))]
      sizes[bit_t %in% "A"] <- bit_v_back[bit_t %in% "A"]
      if(!identical(bit_v_back[types %in% 3L], bit_v[types %in% 3L])) stop("DATA segment: bad $PnB, FCS >= 3.2 spe. Values not divisble by 8 are deprecated")
      masks = sapply(seq_len(n_par), FUN = function(i) {
        if(identical(types[i], 3L)) return(suppressWarnings(as.integer(min(8*bit_n[i],ceiling(log2(as.numeric(text[[paste0("$P",i,"R")]])))))))
        return(c(0L,32L,64L)[types[i] + 1L])
      })
      if(anyNA(c(types, sizes, masks))) stop("DATA segment: instructions for bytes lead to NA")
      args = list(fname = fileName, offset = start, events = n_obj,
                  b_typ = types,
                  b_siz = sizes,
                  b_msk = masks,
                  swap = endian != .Platform$endian)
      if(all((bit_n * 8) == masks)) args$b_msk = integer(0)
      matrix(do.call(what = cpp_readFCS, args = args), ncol = n_par, nrow = n_obj, byrow = TRUE)
    })
    
    if(inherits(data, "try-error")) {
      warning(attr(data, "condition")$message, "\n", "retrying to parse FCS3.2", call. = FALSE, immediate. = TRUE)
      data = readFCSdata(fileName = fileName, text = text, version = 0, start = start, end = end, scale = scale, display_progress = display_progress)
    }
  }
  
  # convert data to data.frame
  feat_names = unlist(recursive = FALSE, use.names = FALSE, lapply(seq_len(n_par), FUN = function(i) {
    N = text[[paste0("$P",i,"N")]]
    S = text[[paste0("$P",i,"S")]]
    if(length(S) != 0) return(paste(N , paste0("< ",S," >")))
    return(N)
  }))
  data = structure(data.frame(data, check.names = FALSE), names = feat_names)
  
  # scale data only for type I, ISAC spe mentions:
  # When linear scale is used, $PnE/0,0/ shall be entered if the floating point data type is used i.e. "F" or "D"
  # meaning that no scaling shall be used for type "F" and "D". Besides type "A" is deprecated
  if(scale) {
    for(i in seq_len(n_par)) { # log amplification scaling
      if(!identical(bit_t[i], "I")) next                  # type of value is not integer, SKIP
      PnE = paste0("$P",i,"E")
      PnR = paste0("$P",i,"R")
      trans = text[[PnE]]
      ran = na.omit(suppressWarnings(as.numeric(text[[PnR]])))
      if((length(trans) == 0) || (length(ran) == 0)) next # no scaling info, SKIP
      trans = na.omit(suppressWarnings(as.numeric(strsplit(trans, split = ",", fixed = TRUE)[[1]])))
      # FIXME should we detect and force range to max range ? in case of mismatch with PnR and actual data range ?
      if(length(trans) != 2 || trans[1] == 0) next        # invalid PnE info, SKIP
      if(trans[2] == 0) trans[2] <- 1                     # invalid PnE info, but we apply correction
      data[,i] <- trans[2] * 10^(trans[1] * data[,i] / ran)
    }
    for(i in seq_len(n_par)) {   # gain scaling
      if((version >= 3.2) && !identical(bit_t[i], "I")) next# gain scaling only apply to "I" in FCS3.2
      PnE = paste0("$P",i,"E")
      PnG = paste0("$P",i,"G")
      gain = na.omit(suppressWarnings(as.numeric(text[[PnG]])))
      if(length(gain) == 0) next # no scaling info, SKIP
      trans = "0,0"
      if(length(text[[PnE]]) != 0) trans = text[[PnE]]
      trans = na.omit(suppressWarnings(as.numeric(strsplit(trans, split = ",", fixed = TRUE)[[1]])))
      # FIXME should we apply PnG when PnE is not valid i.e. PnE = f1 (missing f2)
      if((length(trans) == 2) && (trans[1] == 0) && (trans[2] == 0)) data[,i] <- data[,i] / gain
    }
  }
  return(data)
}

#' @title FCS Dataset Parser
#' @description
#' Helper to parse dataset from Flow Cytometry Standard (FCS) compliant files.
#' @param fileName path to file.
#' @param options list of options used to parse FCS file. It should contain (otherwise, it will be filled with the default values listed below):\cr
#' - header, a list whose members define the "at" offset from header$start$at and the "n" number of bytes to extract:\cr
#' -- start: where to start reading current FCS dataset.       Default is list(at = 0,  n = 6),\cr
#' -- space: where to retrieve space.                          Default is list(at = 6,  n = 4),\cr
#' -- text_beg: where to retrieve file text segment beginning. Default is list(at = 10, n = 8),\cr
#' -- text_end: where to retrieve file text segment end.       Default is list(at = 18, n = 8),\cr
#' -- data_beg: where to retrieve file data segment beginning. Default is list(at = 26, n = 8),\cr
#' -- data_end: where to retrieve file data segment end.       Default is list(at = 34, n = 8),\cr
#' - apply_scale, whether to apply data scaling. It only applies when fcs file is stored as DATATYPE "I". Default is TRUE.\cr
#' - dataset, (coerced to) an ordered vector of unique indices of desired dataset(s) to extract. Default is 1 to extract only the first dataset, whereas NULL allows to extract all available datasets.\cr
#' - force_header, whether to force the use of header to determine the position of data segment. Default is FALSE, for using positions found in "$BEGINDATA" and "$ENDDATA" keywords.\cr
#' - text_only, whether to only extract text segment. Default is FALSE.\cr
#' - text_check, whether to check if text segment is compliant with FCS specifications. Default is FALSE.\cr
#' - text_empty, whether to allow empty values when parsing text segment. Default is FALSE.\cr
#' - text_trim, remove whitespace in keywords names. Default is "none". Allowed are "both", "left", "right" and "none".
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other arguments to be passed.
#' @details 'options' may be tweaked according to file type, instrument and software used to generate it.\cr
#' Default 'options' should allow to read most files.\cr
#' 'options' members with the exception of 'header' may be passed thanks to '...'.
#' @return a list containing:\cr
#' - options, list of 'options' used\cr
#' - header, list of header information corresponding to 'options'\cr
#' - delimiter, unique character used to separate keyword-value pairs\cr
#' - text, list of keywords values,\cr
#' - data, data.frame of values.
#' @keywords internal
readFCSdataset <- function(fileName, options, display_progress = TRUE, ...) {
  dots = list(...)
  # prepare fileName
  if(missing(fileName)) stop("'fileName' can't be missing")
  assert(fileName, len = 1)
  fileName = enc2native(normalizePath(fileName, winslash = "/", mustWork = FALSE))
  if(!file.exists(fileName)) stop(paste("can't find", fileName, sep=" "))
  
  # fill options with default values when not provided
  opt_default = eval(formals(readFCS)$options)
  if(missing(options)) {
    options = opt_default
  } else {
    for(i in c("header", "apply_scale", "dataset", "force_header", "text_only", "text_check", "text_empty")) {
      if(!(i %in% names(options))) options[[i]] <- opt_default[[i]]
    }
    for(i in c("start", "space", "text_beg", "text_end", "data_beg", "data_end")) {
      if(!(i %in% names(options$header))) options$header[[i]] <- opt_default$header[[i]]
    }
  }
  # check if we can find options arguments in dots
  if("text_only" %in% names(dots)) options$text_only <- dots$text_only
  if("text_check" %in% names(dots)) options$text_check <- dots$text_check
  if("text_empty" %in% names(dots)) options$text_empty <- dots$text_empty
  if("text_trim" %in% names(dots)) options$text_trim <- dots$text_trim
  if("dataset" %in% names(dots)) options$dataset <- dots$dataset
  if("apply_scale" %in% names(dots)) options$apply_scale <- dots$apply_scale
  if("force_header" %in% names(dots)) options$force_header <- dots$force_header
  assert(options[["text_only"]], len = 1, alw = c(TRUE, FALSE))
  assert(options[["text_check"]], len = 1, alw = c(TRUE, FALSE))
  assert(options[["text_empty"]], len = 1, alw = c(TRUE, FALSE))
  assert(options[["text_trim"]], len = 1, alw = c("none", "both", "left", "right"))
  options[["dataset"]] = sort(unique(unname(options[["dataset"]])), na.last = FALSE)
  if(length(options[["dataset"]]) == 0) options[["dataset"]] = integer()
  assert(options[["apply_scale"]], len = 1, alw = c(TRUE, FALSE))
  assert(options[["force_header"]], len = 1, alw = c(TRUE, FALSE))
  
  # extract HEADER
  header = readFCSheader(fileName = fileName, header = options$header, encoding = "UTF-8")
  at = attr(header, "offset")
  version = suppressWarnings(na.omit(as.numeric(substr(header$start,4,6))))
  if(length(version) == 0) stop("can't determine FCS version: ", header$start)
  
  # extract DELIMITER
  delimiter = readFCSdelimiter(fileName = fileName, at = header$text_beg, version = version, encoding = "UTF-8")
  
  # extract TEXT segment,
  # the primary TEXT segment has to be in within bytes 58 - 99,999,999
  text = readFCStext(fileName = fileName, delimiter = delimiter, 
                     start = 1 + header$text_beg, end = header$text_end,
                     encoding = "UTF-8", empty = options$text_empty, trim = options$text_trim)
  NTEXT = names(text)
  names(text) = toupper(names(text))
  text_bck = text
  
  # extract supplemental TEXT segment
  # we will use text to extract supplemental TEXT segment offsets
  extra_off1 = suppressWarnings(na.omit(as.numeric(text[["$BEGINSTEXT"]]) + at))
  extra_off2 = suppressWarnings(na.omit(as.numeric(text[["$ENDSTEXT"]])   + at))
  if((length(extra_off1) != 0) &&
     (length(extra_off2) != 0) &&
     (extra_off1 != header$text_beg) && 
     (extra_off2 != header$text_end) &&
     (extra_off2 > extra_off1)) {
    tryCatch({
      extra_text = readFCStext(fileName = fileName, delimiter = delimiter, 
                               start = extra_off1, end = extra_off2,
                               encoding = "UTF-8", empty = options$text_empty, trim = options$text_trim)
      NEXTRA = names(extra_text)
      names(extra_text) = toupper(names(extra_text))
      tmp = names(extra_text) %in% names(text)
      if(any(tmp)) warning("TEXT (or sup. TEXT) segment: supplemental TEXT segment contains keyword",ifelse(sum(tmp)==0,"","s")," already found in TEXT segment that will be discarded:\n\t- ",
                           paste0(NEXTRA[tmp], collapse = "\n\t- "),
                           call. = FALSE, immediate. = TRUE)
      text = c(text, extra_text[!tmp])
      NTEXT = c(NTEXT, NEXTRA[!tmp])
    }, error = function(e) {
      text = text_bck
      warning("TEXT (or sup. TEXT) segment: supplemental TEXT segment is not readable:\n",
              e$message, call. = FALSE, immediate. = TRUE)
    })
  }
  
  # get cur fileName and dataset
  if(!any("$FIL" == names(text))) {
    text[["$FIL"]] <- fileName
    NTEXT = c(NTEXT, "$FIL")
  }
  data_fil <- text[["$FIL"]]
  data_set = num_to_string(length(options$header$start$at))
  
  # check that keywords fulfill FCS spe
  if((options$text_check) &&
     (
       (length(options$dataset) == 0) ||
       (data_set %in% options$dataset)
     )) FCS_check_keywords(text = text, delimiter = delimiter, version = version, encoding = "UTF-8")
  
  # retrieve DATA segment offsets
  # $BEGINDATA and $ENDDATA are not REQUIRED keywords in old (<= 2.0) FCS spe
  # so we use HEADER segment for data offsets and silently set 'force_header' = TRUE,
  # otherwise, for >= 3.0, preferred location is within $BEGINDATA and $ENDDATA from TEXT segment
  off1 = numeric()
  off2 = numeric()
  has_been_forced = FALSE
  to_string = function(x) {if(length(x)==0) { return(NULL) } else { return(num_to_string(x))}}
  if(version >= 3.0) {
    off1 = suppressWarnings(na.omit(as.numeric(text[["$BEGINDATA"]]) + at))
    off2 = suppressWarnings(na.omit(as.numeric(text[["$ENDDATA"]])   + at))
    # if not found in text despite being mandatory, we will force_header
    if(any(c(off1, off2, length(off1), length(off2)) == 0)) {
      if(!options$force_header) warning("can't determine DATA offsets from TEXT keywords $BEGINDATA[",
                                        to_string(off1),"] or $ENDDATA[",to_string(off2),
                                        "], 'force_header' has been forced to TRUE",
                                        call. = FALSE, immediate. = TRUE)
      has_been_forced = TRUE
      options$force_header = TRUE
    }
    # we check consistency between HEADER and TEXT
    if((!has_been_forced) &&                                                     # are we forced to read data offsets from header only
       (!options$text_only) &&                                                   # does user want data to be extracted ?
       ((length(options$dataset) == 0) || (data_set %in% options$dataset)) &&    # should the current dataset be extracted ?
       ((!any(0 %in% header$data_beg)) && (!any(0 %in% header$data_end))) &&     # data offsets can be found in header and is not 0
       ((!any(off1 %in% header$data_beg)) || (!any(off2 %in% header$data_end)))) # data offsets from keywords can be found and differs from header
    {
      warning("discrepancies between HEADER[",
              to_string(header$data_beg),"-",to_string(header$data_end),
              "] and TEXT[",
              to_string(off1),"-",to_string(off2),
              "] segments for determining DATA offsets\n",
              "/!\\ you should manually validate results with 'force_header' set to TRUE or FALSE",
              call. = FALSE, immediate. = TRUE)
    }
  } else {
    options$force_header = TRUE
  }
  if(options$force_header) {
    if(any(c(header$data_beg, header$data_end, length(header$data_beg), length(header$data_end)) == 0)) {
      warning("can't 'force_header' because HEADER indicates 0 for DATA offset(s)",
              call. = FALSE, immediate. = TRUE)
    } else {
      off1 = header$data_beg
      off2 = header$data_end
    }
  }
  if(any(c(off1, off2, length(off1), length(off2)) == 0) || (off2 <= off1)) {
    warning("can't determine data offsets, 'text_only' has been forced to TRUE",
            call. = FALSE, immediate. = TRUE)
    options$text_only = TRUE
  }
  
  # check $MODE, we can only extract DATA segment in "L" mode
  mode = text[["$MODE"]]
  has_been_forced = FALSE
  if(version <= 3.1) { # $MODE is REQUIRED in FCS <= 3.1
    if((length(mode) != 1) || (mode != "L")) has_been_forced = TRUE
  } else {             # $MODE is DEPRECATED starting FCS >= 3.2, but if here it should be "L"
    if((length(mode) == 1) && (mode != "L")) has_been_forced = TRUE
  }
  if(has_been_forced) {
    options$text_only = TRUE
    warning("DATA stored in $MODE[",mode,"] are not supported, 'text_only' has been forced to TRUE", # mode "C" and "U" have been deprecated in FCS spe
            call. = FALSE, immediate. = TRUE)
  }
  
  # extract DATA segment
  data = data.frame()                                 # prepare default returned value for data
  if(!options$text_only &&                            # user only wants text segment
     (                                                
       (length(options$dataset) == 0) ||
       (data_set %in% options$dataset)  # no need to extract data if user doesn't need it
     )) data = readFCSdata(fileName = fileName, text = text, version = version,
                           start = off1, end = off2, 
                           scale = options$apply_scale, display_progress = display_progress)
  
  # TODO retrieve analysis segment ?
  # # we will use text to extract analysis segment offsets
  # off1 = suppressWarnings(as.integer(text[["$BEGINANALYSIS"]]))
  # off2 = suppressWarnings(as.integer(text[["$ENDANALYSIS"]]))
  # # if not found in text despite being mandatory, we will use header
  # if(length(off1) == 0) off1 = suppressWarnings(as.integer(header$data_beg))
  # if(length(off2) == 0) off2 = suppressWarnings(as.integer(header$data_end))
  # anal=raw()
  
  # recover original TEXT names (i.e. not forced to upper case)
  for(k in c("file", "fileName", "fileName_image", "date", "dataset", "version", "FCSversion")) {
    kk = paste0("@IFC_", k)
    foo = toupper(kk) == names(text)
    if(any(foo)) NTEXT[foo] <- kk
  }
  names(text) = NTEXT
  text[["@IFC_file"]] <- basename(data_fil) # internal filename if found, otherwise fileName
  text[["@IFC_fileName"]] <- fileName
  text[["@IFC_dataset"]] <- data_set
  text[["@IFC_version"]] <- paste0(unlist(recursive = FALSE, use.names = FALSE, packageVersion("IFC")), collapse = ".")
  text[["@IFC_FCSversion"]] <- header$start
  return(list(options=options,
              header=header,
              delimiter=delimiter,
              # anal=raw(),
              text=text, 
              data=data))
}

#' @title FCS File Parser
#' @description
#' Parse data from Flow Cytometry Standard (FCS) compliant files.
#' @param fileName path to file.
#' @param options list of options used to parse FCS file. It should contain (otherwise, it will be filled with the default values listed below):\cr
#' - header, a list whose members define the "at" offset from header$start$at and the "n" number of bytes to extract:\cr
#' -- start: where to start reading current FCS dataset.       Default is list(at = 0,  n = 6),\cr
#' -- space: where to retrieve space.                          Default is list(at = 6,  n = 4),\cr
#' -- text_beg: where to retrieve file text segment beginning. Default is list(at = 10, n = 8),\cr
#' -- text_end: where to retrieve file text segment end.       Default is list(at = 18, n = 8),\cr
#' -- data_beg: where to retrieve file data segment beginning. Default is list(at = 26, n = 8),\cr
#' -- data_end: where to retrieve file data segment end.       Default is list(at = 34, n = 8),\cr
#' - apply_scale, whether to apply data scaling. It only applies when fcs file is stored as DATATYPE "I". Default is TRUE.\cr
#' - dataset, (coerced to) an ordered vector of unique indices of desired dataset(s) to extract. Default is 1 to extract only the first dataset, whereas NULL allows to extract all available datasets.\cr
#' - force_header, whether to force the use of header to determine the position of data segment. Default is FALSE, for using positions found in "$BEGINDATA" and "$ENDDATA" keywords.\cr
#' - text_only, whether to only extract text segment. Default is FALSE.\cr
#' - text_check, whether to check text segment is compliant with FCS specifications. Default is FALSE.\cr
#' - text_empty, whether to allow empty values when parsing text segment. Default is FALSE.\cr
#' - text_trim, remove whitespace in keywords names. Default is "none". Allowed are "both", "left", "right" and "none".
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other arguments to be passed.
#' @details 'options' may be tweaked according to file type, instrument and software used to generate it.\cr
#' Default 'options' should allow to read most files.\cr
#' 'options' members with the exception of 'header' may be passed thanks to '...'.\cr
#' Experimental (as of v0.2.1.300), readFCS could handle FCS 3.2 files. However, it is important to note that R has no native support for 64bits unsigned integers which are defined in the FCS 3.2 standard.
#' So, those integers are extracted as double (8 bytes) and precision loss will happen for > 2^53 integers on 64bits platforms.
#' @source Data File Standard for Flow Cytometry, version FCS 3.2 from Spidlen J. et al. available at \doi{10.1002/cyto.a.24225}.
#' @return a list whose elements are lists for each dataset stored within the file.\cr
#' each sub-list contains:\cr
#' - header, list of header information corresponding to 'options'\cr
#' - delimiter, unique character used to separate keyword-value pairs\cr
#' - text, list of keywords values,\cr
#' - data, data.frame of values.
#' @export
readFCS <- function(fileName,
                    options = list(header = list(start = list(at = 0, n = 6),
                                                 space = list(at = 6, n = 4),
                                                 text_beg = list(at = 10, n = 8),
                                                 text_end = list(at = 18, n = 8),
                                                 data_beg = list(at = 26, n = 8),
                                                 data_end = list(at = 34, n = 8)),
                                   apply_scale = TRUE,
                                   dataset = 1,
                                   force_header = FALSE,
                                   text_only = FALSE,
                                   text_check = FALSE,
                                   text_empty = FALSE,
                                   text_trim = "none"),
                    display_progress = TRUE, ...) {
  # extract dataset(s)
  ans = list()
  more = 0L
  while(length(more) != 0) {
    dat = readFCSdataset(fileName = fileName, options = options, display_progress = display_progress, ...)
    options <- dat$options
    ans = c(ans, list(dat[-1]))
    text = ans[[length(ans)]]$text
    names(text) = toupper(names(text))
    fileName = text[["@IFC_FILENAME"]]
    more = suppressWarnings(na.omit(as.numeric(text[["$NEXTDATA"]]) + options$header$start$at[1]))
    more = setdiff(more, options$header$start$at)
    if((length(options$dataset) != 0) && all(options$dataset %in% 1L)) more = numeric()
    if((length(more) != 0) && (more >= file.size(fileName))) {
      more = numeric()
      warning("can't extract all datasets: $NEXTDATA points to outside of the file",
              call. = FALSE, immediate. = TRUE)
    }
    options$header$start$at <- c(more, options$header$start$at)
  }
  if(length(options$dataset) == 0) options$dataset = seq_along(ans)
  tmp = options$dataset %in% seq_along(ans)
  if(!all(tmp)) warning("dataset",
                        ifelse(sum(!tmp)==1,"","s"),
                        " [",
                        paste0(options$dataset[!tmp], collapse = ",")
                        ,"] can't be found in\n",
                        fileName,
                        call. = FALSE, immediate. = TRUE)
  return(structure(ans[sapply(seq_along(ans), FUN = function(ii) ans[[ii]]$text[["@IFC_dataset"]] %in% options$dataset)],
                   class = "IFC_fcs", fileName = fileName))
}

#' @title FCS Object Data Sets Merging
#' @description
#' Merges FCS data object with various data sets.
#' @param fcs `IFC_fcs` object as extracted by readFCS().
#' @param ... other arguments to be passed.
#' @details in data can contain extra columns named 'import_file' and 'import_subfile' intended to allow file/dataset identification
#' @return a list of list containing:\cr
#' - header, list of header information corresponding to 'options'\cr
#' - delimiter, unique character used to separate keyword-value pairs\cr
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
      pb = newPB(label = "FCS", min = 0, max = L)
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
      if(length(com) != 0) {
        if(length(c(Nxx, Nyy)) == 0) {
          aa = structure(rbind.data.frame(x[, com, drop = FALSE],
                                          y[, com, drop = FALSE],
                                          make.row.names = FALSE),
                         names = com)
        } else {
          aa = structure(cbind.data.frame(aa, rbind.data.frame(x[, com, drop = FALSE],
                                                               y[, com, drop = FALSE],
                                                               make.row.names = FALSE)),
                         names = c(Nxx, Nyy, com))
        }
      }
      aa
    },
    lapply(seq_len(L), FUN = function(i) {
      if(display_progress) setPB(pb, value = i, title = "FCS", label = "Merging Data Sets")
      dat = fcs[[i]]$data
      if(!any("import_file" == names(dat))) {
        dat[,"import_file"]=rep(fcs[[i]]$text[["@IFC_file"]], nrow(dat))
        # dat = cbind.data.frame(dat, "import_file"=rep(fcs[[i]]$text[["@IFC_file"]], nrow(dat)))
      }
      if(!any("import_subfile" == names(dat))) {
        dat[,"import_subfile"]=rep(fcs[[i]]$text[["@IFC_dataset"]], nrow(dat))
        # dat = cbind.data.frame(dat, "import_subfile"=rep(fcs[[i]]$text[["@IFC_dataset"]], nrow(dat)))
      }
      dat
    }))
  } else {
    features = fcs[[1]]$data
    if(!any("import_file" == names(features))) {
      features[,"import_file"]=rep(fcs[[1]]$text[["@IFC_file"]], nrow(features))
      # features = cbind.data.frame(features, "import_file"=rep(fcs[[1]]$text[["@IFC_file"]], nrow(features)))
    }
    if(!any("import_subfile" == names(features))) {
      features[,"import_subfile"]=rep(1, nrow(features))
      # features = cbind.data.frame(features, "import_subfile"=rep(1, nrow(features)))
    }
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
#' - delimiter, unique character used to separate keyword-value pairs\cr
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
      pb = newPB(label = "FCS", min = 0, max = L)
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
    lapply(seq_len(L), FUN = function(i) {
      if(display_progress) setPB(pb, value = i, title = "FCS", label = "Merging Samples")
      dat = fcs[[i]]$data
      if(!any("import_file" == names(dat))) {
        dat[,"import_file"]=rep(fcs[[i]]$text[["@IFC_file"]], nrow(dat))
        # dat = cbind.data.frame(dat, "import_file"=rep(fcs[[i]]$text[["@IFC_file"]], nrow(dat)))
      }
      if(!any("import_subfile" == names(dat))) {
        dat[,"import_subfile"]=rep(fcs[[i]]$text[["@IFC_dataset"]], nrow(dat))
        # dat = cbind.data.frame(dat, "import_subfile"=rep(fcs[[i]]$text[["@IFC_dataset"]], nrow(dat)))
      }
      dat
    }))
  } else {
    features = fcs[[1]]$data
    if(!any("import_file" == names(features))) {
      features[,"import_file"]=rep(fcs[[1]]$text[["@IFC_file"]], nrow(features))
      # features = cbind.data.frame(features, "import_file"=rep(fcs[[1]]$text[["@IFC_file"]], nrow(features)))
    }
    if(!any("import_subfile" == names(features))) {
      features[,"import_subfile"]=rep(1, nrow(features))
      # features = cbind.data.frame(features, "import_subfile"=rep(1, nrow(features)))
    }
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

#' @title FCS Object Converter
#' @description
#' Converts FCS data object to `IFC_data` object.
#' @param fcs `IFC_fcs` object as extracted by readFCS().
#' @param ... other arguments to be passed.
#' @details in data can contain extra columns named 'import_file' and 'import_subfile' intended to allow file/dataset identification
#' @return A named list of class `IFC_data`, whose members are:\cr
#' -description, a list of descriptive information,\cr
#' -Merged_fcs, character vector of path of files used to create fcs, if input was a merged,\cr
#' -Keywords, a named-list of keywords values, only keywords from 1st 'fcs' segment will be retrieved\cr
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
                                       "masks" = data.frame(matrix(character(),nrow = 0, ncol = 3, dimnames = list(list(),list("type","name","def"))))),
                  "Merged_fcs" = character(),
                  "Keywords" = list(),
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
    idx$count = seq_len(obj_count)
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
    names(x$text) = toupper(names(x$text))
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
    names(x$text) = toupper(names(x$text))
    tmp = x$text[c("$SPILLOVER","SPILL","SPILLOVER")]
    tmp = tmp[sapply(tmp, length) != 0]
    if(length(tmp) == 0) return(NULL)
    return(tmp)
  })
  spillover = spillover[sapply(spillover, length) != 0]
  if(is.list(spillover) && length(spillover) == 1) spillover = spillover[[1]]
  if(is.list(spillover) && length(spillover) == 1) spillover = spillover[[1]]
  if(!is.list(spillover) && length(spillover) == 1 && spillover == "") spillover = NULL
  if(length(spillover) != 0) {
    features_names = parseFCSname(names(features))
    spillover = try(convert_spillover(spillover), silent = TRUE)
    if(inherits(spillover, "try-error")) {
      spillover = NULL
    } else {
      rownames(spillover) <- names(features)[apply(sapply(colnames(spillover), FUN = function(x) {
        x==features_names$PnN
      }), 2, FUN = function(i) which(i)[1])]
    }
  }
  # checksum = sapply(fcs, FUN = function(x) {
  #   tmp = x$description[[1]]$`$ENDDATA`
  #   if((length(tmp) == 0) || (tmp == "")) return("unk")
  #   return(tmp)
  # })
  
  min_data$fileName = normalizePath(attr(fcs, "fileName"), winslash = "/", mustWork = FALSE)
  bar <- unique(idx[, "import_file"])
  if(length(bar) > 1) min_data$Merged_fcs <- bar
  min_data$Keywords <- fcs[[1]]$text
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
  min_data = suppressWarnings(IFC::data_add_features(obj = min_data, features = features_def))
  min_data = IFC::data_add_pops(obj = min_data,
                                pops = list(buildPopulation(name = "All", type = "B",
                                                            color = "White", lightModeColor = "Black",
                                                            obj = rep(TRUE, obj_count))),
                                display_progress = display_progress)
  # min_data$features_comp = min_data$features[, grep("^FS.*$|^SS.*$|LOG|^Object Number$|TIME", names(min_data$features), value = TRUE, invert = TRUE, ignore.case = TRUE)]
  if(multiple) {
    min_data = IFC::data_add_pops(obj = min_data, pops = pops, display_progress = display_progress)
  }
  K = class(min_data$pops)
  min_data$pops = lapply(min_data$pops, FUN = function(p) {
    attr(p, "reserved") <- TRUE
    return(p)
  })
  class(min_data$pops) <- K
  min_data$stats = get_pops_stats(min_data$pops, obj_count)
  return(min_data)
}

#' @title FCS File Reader
#' @description
#' Extracts data from Flow Cytometry Standard (FCS) Files.
#' @param fileName path(s) of file(s). If multiple files are provided they will be merged and 
#' populations will be created to identify each single file within returned `IFC_data` object.
#' @source Data File Standard for Flow Cytometry, version FCS 3.2 from Spidlen J. et al. available at \doi{10.1002/cyto.a.24225}.
#' @param ... other arguments to be passed to readFCS function, with the exception of 'options$text_only'.\cr
#' Experimental (as of v0.2.1.300), ExtractFromFCS could handle FCS 3.2 files. However, it is important to note that R has no native support for 64bits unsigned integers which are defined in the FCS 3.2 standard.
#' So, large integers are extracted as double (8 bytes) and precision loss will happen for > 2^53 integers on 64bits platforms.
#' @return A named list of class `IFC_data`, whose members are:\cr
#' -description, a list of descriptive information,\cr
#' -Merged_fcs, character vector of path of files used to create fcs, if input was a merged,\cr
#' -Keywords, a named-list of keywords values, only keywords from 1st 'fcs' segment will be retrieved\cr
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
  dots = dots[!(names(dots) %in% c("fcs","text_only"))]
  display_progress = dots$display_progress
  if(length(display_progress) == 0) display_progress = TRUE
  assert(display_progress, len=1, alw = c(TRUE, FALSE))
  fileName = enc2native(normalizePath(path = fileName, winslash = "/", mustWork = TRUE))
  
  # read the fcs file and extract features and description
  L = length(fileName)
  if(display_progress) {
    pb = newPB(label = "reading files", min = 0, max = L)
    on.exit(endPB(pb))
  }
  fcs = lapply(seq_len(L), FUN = function(i_file) {
    if(display_progress) setPB(pb, value = i_file, title = "Extracting FCS", label = "reading files")
    do.call(what = FCS_merge_dataset, args = c(dots, list(fcs = quote(do.call(what = readFCS,  args=c(dots, list(fileName = fileName[[i_file]])))))))[[1]]
  })
  attr(fcs, "fileName") <- fileName[1]
  do.call(what = FCS_to_data, args = c(dots, list(fcs = quote(fcs))))
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
  locale_back <- setloc(c("LC_ALL" = "en_US.UTF-8"))
  enc_back <- options("encoding" = "UTF-8")
  on.exit(suspendInterrupts({setloc(locale_back); options(enc_back)}), add = TRUE)
  now = format(Sys.time(), format = "%d-%b-%y %H:%M:%S")
  
  old_enc <- options("encoding")
  on.exit(options(old_enc), add = TRUE)
  options("encoding" = "UTF-8")
  
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
  if((raw_delimiter == 0x00) || (raw_delimiter > 0x7E)) stop("'delimiter' should be an [0x01-0x7E (1-126)] ASCII character")
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
      fcs = readFCS(fileName = write_to, text_only = TRUE)
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
  pkg_ver = paste0(unlist(recursive = FALSE, use.names = FALSE, packageVersion("IFC")), collapse = ".")
  # is_fcs = length(obj$description$FCS)!=0
  
  # init header
  header = list(version  = "FCS3.0",
                space    = "    ",
                text_beg = "      58",
                text_end = "        ",
                data_beg = "        ",
                data_end = "        ",
                anal_beg = "       0",
                anal_end = "       0")
  
  # we modify features to add populations
  features = fastCbind(obj$features[, setdiff(names(obj$features), names(obj$pops)), drop = FALSE],
                       sapply(names(obj$pops), simplify = FALSE, FUN = function(p) obj$pops[[p]]$obj))
  # need to replace non finite values by something; IDEAS is using 0 so we use 0 also
  # TODO maybe replace -Inf by features min and +Inf by features max ?
  features = as.data.frame(apply(features, 2, cpp_replace_non_finite), stringsAsFactors = TRUE)
  
  # determines length of text_segment2
  # comma (ASCII 0x2C) is not allowed in features names according to fcs specifications so it is replaced by a space
  feat_names = parseFCSname(gsub(","," ",names(features),fixed=TRUE))
  N = feat_names$PnN
  tmp = duplicated(toupper(N))
  if(any(tmp)) stop("$PnN should be unique\n\t- ", paste0(N[tmp], collpase="\n\t- "))
  S = feat_names$PnS
  # S[S == ""] <- names(features)[S == ""] # should we add PnS when not here ? IDEAS does not export $PnS
  L = length(N)
  text2_length = 0
  text_segment2 = lapply(seq_len(L), FUN = function(i) {
    # each time we write /$PnN/<feature_name>/$PnB/32/$PnE/0, 0/$PnG/1/$PnR/<feature_max_value> and /$PnS/<feature_alt-name> if PnS is not empty
    # bar = 1 + diff(range(features[, i], na.rm = TRUE))
    bar = ceiling(max(features[, i], na.rm = TRUE)) # FIXME since we write type "F" this has no importance, shall we use 262144, as it is commonly used ?
    foo = c(paste0("$P",i,"N"), N[i],
            paste0("$P",i,"B"), "32",
            paste0("$P",i,"E"), "0, 0",
            paste0("$P",i,"G"), "1",
            paste0("$P",i,"R"), bar)
    if(S[i] != "") foo = c(foo, paste0("$P",i,"S"), S[i])
    if(any(sapply(foo, FUN = function(x) substr(x,1,1) == delimiter))) stop("keyword-value pairs should not start with 'delimiter'[",delimiter,"]:\n\t- ",
                                                                            paste0(paste0(foo[rep_len(c(TRUE, FALSE), length(foo))], "[", foo[rep_len(c(FALSE, TRUE), length(foo))], "]"), collapse="\n\t- "))
    foo = gsub(pattern = delimiter, x = foo, replacement = delimiter_esc, fixed=TRUE)
    foo = charToRaw(paste(c("", foo), collapse = delimiter))
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
                       "$TOT" = num_to_string(obj$description$ID$objcount),                      #*
                       "$PAR" = L,                                                               #*
                       #* BEGINDATA and ENDDATA are also mandatory and will be added afterwards
                       #* PnB, PnE, PnN, and PnR are also mandatory and are part of text_segment2
                       #* PnG is not mandatory but will be filled in part 2
                       "$CYT" = cyt,
                       "@IFC_fileName" = obj$fileName,
                       "@IFC_fileName_image" = obj$fileName_image,
                       "@IFC_version" = pkg_ver,
                       "@IFC_date" = now
  )
  
  # gather keywords in priority order 
  text_segment1 = c(text_segment1, dots, obj$Keywords)
  # removes keywords whose values are NULL
  text_segment1 = text_segment1[sapply(text_segment1, length) != 0]
  # removes duplicated keywords (priority order is important here)
  text_segment1 = text_segment1[!duplicated(toupper(names(text_segment1)))]
  # removes not allowed keywords (e.g. in text_segment2 ($PnN, $PnS, $PnB, $PnE, $PnG, $PnR) or "")
  text_segment1 = text_segment1[setdiff(names(text_segment1),c(""))]
  text_segment1 = text_segment1[!grepl("^\\$P\\d+[N|S|B|E|G|R]$", names(text_segment1), ignore.case = TRUE)]
  
  # determines length of data
  data_length = 4 * L * nrow(features)
  
  # determines length of mandatory param
  N = names(text_segment1)
  text1_length = sum(c(nchar("$BEGINDATA"), 2, # 2 for additional delimiters
                       nchar("$ENDDATA"),   2, # 2 for additional delimiters, there is already one at the beg of text2
                       sapply(seq_along(text_segment1), FUN = function(i) {
                         foo = c(N[i], text_segment1[[i]])
                         if(any(sapply(foo, FUN = function(x) substr(x,1,1) == delimiter))) stop("keyword-value pairs should not start with 'delimiter' [",delimiter,"]:\n\t- ",
                                                                                                 paste0(foo[1],"[",foo[2],"]"))
                         v = charToRaw(gsub(delimiter, delimiter_esc, N[i], fixed=TRUE))
                         if(any(v < 0x20 | v >= 0x7F)) stop("keyword contains invalid ASCII character, valid are [0x20-0x7E (32-126)]\n\t- ",
                                                            N[i], "[",paste0(paste0("0x",v), collapse = ","),"]")
                         length(charToRaw(gsub(delimiter, delimiter_esc, text_segment1[[i]], fixed=TRUE))) +
                           length(v) + 2 # 2 for additional delimiters
                       }), text2_length,
                       nchar(paste0(header, collapse = ""))
  ))
  
  # compute missing offsets 
  # ENDSTEXT
  # determining text_end is tricky since it depends on its own length
  # so we use a function to optimize it
  f = function(x, text_length) {
    data_beg = x + 1
    data_end = x + data_beg + data_length - 1
    ans = text_length + nchar(num_to_string(data_beg)) + nchar(num_to_string(data_end))
    if(ans != x) ans = f(x = ans, text_length = text_length)
    return(ans)
  }
  text_end = f(x = text1_length, text_length = text1_length)
  if(text_end >= 1e8) stop("primary TEXT segment is too big")
  header$text_end = sprintf("%8i", text_end)
  
  # BEGINDATA / ENDDATA
  data_beg = text_end + 1               # +1 because data start just after text segment end
  data_end = data_beg + data_length - 1 # -1 because last data byte is at minus one from total length
  if((data_beg >= 1e8) || (data_end >= 1e8)) {
    header$data_beg = sprintf("%8i", 0)
    header$data_end = sprintf("%8i", 0)
  } else {
    header$data_end = sprintf("%8i", data_end)
    header$data_beg = sprintf("%8i", data_beg)
  }
  text_segment1 = c(text_segment1, list("$BEGINDATA" = num_to_string(data_beg)))
  text_segment1 = c(text_segment1, list("$ENDDATA"   = num_to_string(data_end)))
  
  towrite = file(description = file_w, open = "wb")
  tryCatch({
    # writes header
    writeBin(object = charToRaw(paste0(header, collapse="")), con = towrite)
    
    # writes text segment1
    N = names(text_segment1)
    lapply(seq_along(text_segment1), FUN = function(i) {
      writeBin(object = charToRaw(paste(c("", 
                                          gsub(delimiter, delimiter_esc, N[i], fixed=TRUE),
                                          gsub(delimiter, delimiter_esc, text_segment1[i], fixed=TRUE)),
                                        collapse = delimiter)), con = towrite)
    })
    
    # writes text segment2
    lapply(seq_along(text_segment2), FUN = function(i) {
      writeBin(object = text_segment2[[i]], con = towrite)
    })
    writeBin(object = charToRaw(delimiter), con = towrite) # we add final delimiter after the last keyword-value
    
    # export features values in data segment
    apply(features, 1, FUN = function(x) writeBin(x, con = towrite, size = 4))
    
    # FIXME write CRC
    writeBin(object = rep(as.raw(0x30), 8), con = towrite) 
    
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
