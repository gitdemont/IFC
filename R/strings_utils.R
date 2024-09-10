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

#' @title Get Locale
#' @description Retrieves current R locale
#' @return a named vector whose members are current locale value.
#' @keywords internal
getloc <- function() {
  sapply(setdiff(.LC.categories, "LC_ALL"), Sys.getlocale)
}

#' @title Set Locale
#' @description Sets aspects of locale
#' @param locale a named vector with locale values to modify. Default is getloc().
#' @return a named vector of the locale values used before the application of the function.
#' @keywords internal
setloc <- function(locale = getloc()) {
  cur_loc <- getloc()
  suppressWarnings(sapply(names(locale), FUN = function(x) Sys.setlocale(x, locale[x])))
  return(cur_loc)
}

#' @title String Truncation
#' @description Truncates character strings
#' @param x a string
#' @param n desired length
#' @details x will be truncated according to 'n' parameter. If x is longer than n '...' are appended.
#' @keywords internal
trunc_string = function(x, n=22) {
  if(length(x) == 0) return("")
  x = as.character(x)
  L = nchar(x,allowNA=TRUE,keepNA=FALSE,type="chars")
  if(L > n) return(paste0(substr(x, 1, n),"..."))
  return(x)
}

#' @title String Sequence Replacement
#' @description Replaces a sequence of strings
#' @param x,pattern,replacement non empty character vectors.
#' @param all whether to replace all instances of 'pattern' or only the 1st one. Default is TRUE, to replace all instances.
#' @details if 'pattern' is found within 'x', 'pattern' will be removed from 'x' and replace by 'replacement'.\cr
#' It looks like gsub but it is different e.g.:\cr
#' x=c("ABD","A","B")\cr
#' pattern=c("A","B")\cr
#' replacement=c("C")\cr
#' - gsub(x=paste0(x,collapse=""),pattern=paste0(pattern,collapse=""),replacement=paste0(replacement,collapse=""),fixed=TRUE) will give "CDC",\cr
#' - gseq(x=x,pattern=pattern,replacement=replacement) will give "ABD","C".
#' @return 'x' where 'pattern' is replaced by replacement.
#' @keywords internal
gseq <- function(x, pattern = "", replacement = character(), all = TRUE) {
  LX <- length(x)
  LM <- length(pattern)
  if(all(c(pattern %in% replacement, replacement %in% pattern))) return(x)
  if(LX >= LM) {
    pos <- cpp_seqmatch(x, pattern)                               # return the position of 1st element found
    if(pos) {                                                     # if something found (i.e. pos != 0)
      a = x[seq_along(integer(pos - 1))];                         # take beg
      b = x[LM + pos - 1 + seq_along(integer(1 + LX - LM - pos))] # take end
      x = c(a, replacement, b)                                    # replacement is inserted
      if(all && (length(x) > 0)) x = gseq(x, pattern)             # another run is launched to remove other match
    }
  }
  x
}

#' @title File Extension Removal
#' @description Removes file extension from file path
#' @param x a file path
#' @details file extension will be removed
#' @keywords internal
remove_ext <- function(x) {
  if(length(x) == 0) return(x)
  xx = as.character(x)
  sapply(1:length(xx), FUN = function(i) gsub(paste0("\\.",getFileExt(xx[i]),"$"), "", xx[i]))
}

#' @title Reverse String
#' @description Reverses string.
#' @param x a character vector.
#' @keywords internal
rev_string <- function(x) {
  if(length(x) == 0) return(x)
  sapply(strsplit(as.character(x), split = "", fixed = TRUE), FUN = function(i) paste0(rev(i), collapse=""))
}

#' @title FCS Name Parser
#' @description 
#' Separates names and alt-names from FCS features names.
#' @param x a character vector of FCS features names.
#' @details when created FCS features names are formatted as 'PnN < PnS >'.
#' The current function allows the separation between PnN and PnS.
#' @return a 2 columns data.frame of names (PnN) and alt-names (PnS).
#' @keywords internal
parseFCSname <- function(x) {
  if(length(x) == 0) return(structure(as.data.frame(matrix(character(), ncol=2, nrow=0)), names = c("PnN", "PnS")))
  foo = strsplit(rev_string(x), split = " < ", fixed = TRUE)
  structure(as.data.frame(t(sapply(foo, FUN = function(x) {
    if((length(x) < 2) || (substr(x[1], 1, 2) != "> ")) return(c(rev_string(paste0(x, collapse = " < ")), ""))
    rev_string(c(paste0(x[-1], collapse=" < "), substr(x[1], 3, nchar(x[1]))))
  })), stringsAsFactors = FALSE), names = c("PnN", "PnS"))
}

#' @title Special Character Replacement
#' @description
#' Helper to replace special character.
#' @param string string where specials will be replaced if found.
#' @param replacement string replacement. Default is "_".
#' @param specials Default is '[\\|\\/|\\:|\\*|\\?|\\"|\'|\\<|\\>|\\|]'.
#' @keywords internal
specialr <- function(string = "", replacement = "_", specials = '[\\|\\/|\\:|\\*|\\?|\\"|\'|\\<|\\>|\\|]') {
  assert(replacement, len = 1, typ = "character")
  assert(specials, len = 1, typ = "character")
  if(grepl(pattern = specials, x = replacement, perl = FALSE)) stop("'replacement' can't contain 'specials'")
  return(gsub(pattern = specials, replacement = replacement, x = string, perl = FALSE))
}

#' @title Name Protection
#' @description
#' Helper to protect population/region name.
#' @param name population/region names
#' @keywords internal
protectn <- function(name) {
  assert(name, typ="character")
  foo = gsub("(.)", "\\\\\\1", name, perl=FALSE)
  paste0("(",paste0(sapply(foo, FUN = function(i) {
    return(paste0("[", i, "]"))
  }), collapse = "|"),")")
}

#' @title Attributes Keeper
#' @description 
#' Helper to keep attributes while using IFC build functions.
#' @param x list to pass to \code{'what'}.
#' @param what a build function.
#' @details all attributes of \code{'x'} will be used with the exception of "names", "dims", and "dimnames".
#' @return result of applying \code{'what'} on \code{'x'} with preserved attributes from \code{'x'}.
#' @keywords internal
keep_attributes <- function(x, what) {
  ans <- do.call(what, x)
  for(k in setdiff(names(attributes(x)),c("names","dim","dimnames"))) { attr(ans, k) <- attr(x, k) }
  ans
}

#' @title Sync Name Decomposition
#' @description
#' Helper that will split text into chunks.
#' @param x string
#' @param name whether to return name only. Default is \code{TRUE}.
#' @param split string used for splitting. Default is \code{"|"}.
#' @details \code{'x'} is supposed to be a string composed of a \code{'name'} followed by \code{'split'} and ending with a string that does not contain \code{'split'}.
#' @return when \code{'x'} is of length 0, \code{""} is returned. Otherwise, the \code{'name'} when \code{'name'} is \code{TRUE} or the \code{'split'} non-containing remaining string when \code{FALSE}.
#' @keywords internal
splits <- function(x, name = TRUE, split = "|") {
  if(length(x) == 0) return("")
  if(length(x) != 1) stop("'x' should be scalar")
  if(length(split) != 1) stop("'split' should be scalar")
  if(!grepl(split, x, fixed = TRUE)) stop("can't find 'split' in 'x'")
  s = strsplit(x, split, fixed = TRUE)[[1]]
  if(name) {
    if(substr(x, nchar(x), nchar(x)) == split) return(substr(x, 1, nchar(x) - 1))
    return(paste0(head(s, -1), collapse = split))
  }
  if(substr(x, nchar(x), nchar(x)) == split) return(character())
  return(tail(s, 1))
}

#' @title String Decomposition with Operators
#' @description
#' Helper that will split definition into chunks of names and operators.
#' @param definition definition to be split
#' @param all_names the names of all allowed names.
#' @param alt_names vector of same length as 'all_names' to use for substitution. It can be used to speed up the process.
#' @param operators operators used. Default is c("And", "Or", "Not", "(", ")").
#' @param split string used for splitting. Default is "|".
#' @param scalar whether to allow presence of scalar or not. Default is FALSE.
#' @param dsplit whether to allow presence of duplicated occurrences of 'split' or not. Default is FALSE.
#' @keywords internal
splitn <- function(definition, all_names, alt_names, operators = c("And", "Or", "Not", "(", ")"),
                   split = "|", scalar = FALSE, dsplit = FALSE) {
  assert(definition, len=1, typ="character")
  assert(all_names, typ="character")
  assert(operators, typ="character")
  assert(split, len=1, typ="character")
  assert(scalar, alw=c(TRUE,FALSE))
  assert(dsplit, alw=c(TRUE,FALSE))
  
  # we create a mapping between all_names and random names
  # we also ensure that random names will not contain any specials and 
  # will not match with themselves nor with names or operators or split to substitute
  # we also order to_substitute by number of character (decreasing) to
  # be sure that longer words will be substitute first
  dups = duplicated(all_names)
  to_substitute = all_names[!dups]
  to_substitute = to_substitute[order(nchar(to_substitute), decreasing = TRUE)]
  if((length(to_substitute) == 0) || (definition == "")) return(definition)
  if(!missing(alt_names)) {
    if(length(all_names) != length(alt_names)) stop("when provided 'alt_names' should be of same length as 'all_names'")
    replace_with = unname(alt_names[!dups])
  } else {
    replace_with = gen_altnames(x = to_substitute, n = max(9,nchar(operators))+1, forbidden = c(to_substitute, operators, split))
  }
  
  # we substitute names and operators with random names
  ans = definition
  for(i in seq_along(to_substitute)) { 
    ans = gsub(pattern = to_substitute[i], replacement = replace_with[i], x = ans, fixed = TRUE)
  }
  
  # we can now split the definition since random names we use
  # do not contain split
  ans = strsplit(ans, split = split, fixed = TRUE)[[1]]
  
  # finally, can replace random names with their corresponding initial names
  # from to_substitute (i.e. all_names + operators)
  ans = sapply(ans, USE.NAMES = FALSE, FUN = function(x) {
    if(x == "") {
      if(dsplit) return("")
      stop(call. = FALSE, 'definition ["',definition,'"] contains duplicated \'split\' [',split,']',
           ifelse(!dsplit, ". Consider the use 'dsplit' = TRUE", ""))
    }
    foo = to_substitute[x == replace_with]
    if(length(foo) == 0) {
      foo = operators[x == operators]
      if(length(foo) == 0) {
        if(scalar && length(na.omit(suppressWarnings(as.numeric(x))))!= 0) return(x)
        stop(call. = FALSE, 'definition ["',definition,'"] contains unexpected name [',x,']',
             ifelse(!scalar && length(na.omit(suppressWarnings(as.numeric(x)))) != 0, ". Consider the use 'scalar' = TRUE", ""))
      }
    }
    return(foo)
  })
  return(ans)
}

#' @title String Decomposition with Placeholders
#' @description
#' Helper aiming to detect placeholder pattern
#' @param write_to string. Default is "\%d/\%s_fromR.\%e"
#' @details 
#' -\%s: shortname (i.e. basename without extension)\cr
#' -\%p: first parent directory\cr
#' -\%d: full path directory\cr
#' -\%e: file extension\cr
#' -\%o: object id\cr
#' -\%c: channel
#' @keywords internal
splitp = function(write_to = "%d/%s_fromR.%e") {
  assert(write_to, len = 1, typ = "character")
  foo = strsplit(write_to, split = "%", fixed = TRUE)[[1]]
  if(length(foo) > 1) {
    pre = foo[1]
    foo = gsub("^(s|p|d|e|o|c)(.*)$", "\\1%\\2", foo[-1])
    foo = unlist(strsplit(foo, split = "%", fixed = TRUE))
    out = lapply(c("d","p","s","e","o","c"), FUN=function(char) {
      which(foo == char)
    })
  } else {
    pre = ""
    out = list(0,0,0,0,0,0)
  }
  names(out) = c("dir", "parent", "short", "ext", "object", "channel")
  names(foo) <- rep("", length(foo))
  names(foo)[out[["dir"]]] <- "dir"
  names(foo)[out[["parent"]]] <- "parent"
  names(foo)[out[["short"]]] <- "short"
  names(foo)[out[["ext"]]] <- "ext"
  names(foo)[out[["object"]]] <- "object"
  names(foo)[out[["channel"]]] <- "channel"
  foo = c(pre, foo)
  out = c(out, decomp = list(foo))
  attr(out, "class") <- "splitp_obj"
  return(out)
}

#' @title File Path Decomposition
#' @description
#' Helper that will split file name into chunks
#' @param file path to file
#' @return a named vector with chunks of 'file'\cr
#' dir: full path directory of 'file'\cr
#' parent: first parent directory of 'file'\cr
#' ext: 'file' extension without leading dot\cr
#' short: 'file' with no extension nor dir\cr
#' input: 'file' path as it was provided.
#' @keywords internal
splitf <- function(file = NULL) {
  b_name = basename(file)
  dir = dirname(file)
  if(dir == "") {
    dir = suppressWarnings(normalizePath(file, mustWork = FALSE, winslash = "/"))
  } else {
    dir = suppressWarnings(normalizePath(dir, mustWork = FALSE, winslash = "/"))
  }
  ext = getFileExt(file)
  short = gsub(paste0("\\.", ext, "$"), "", b_name, ignore.case = TRUE)
  out = c("dir" = dir, "parent" = basename(dir), "ext" = ext, "short" = short, "input" = file)
  class(out) <- "splitf_obj"
  return(out)
}

#' @title File Path Placeholders Formatting
#' @description
#' Helper to format splitp_obj using splitf_obj, channel and object information.
#' @param splitp_obj object returned by \code{\link{splitp}}. 
#' @param splitf_obj object returned by \code{\link{splitf}}. It will be used to substitute \%d, \%p, \%s and \%e.
#' @param channel string to be used to substitute \%c
#' @param object string to be used to substitute \%o
#' @keywords internal
formatn <- function(splitp_obj, splitf_obj, channel = "", object = "") {
  if(missing(splitf_obj)) {
    splitf_obj = list(dir = "", parent = "", file = "", ext = "")
    class(splitf_obj) <- c("splitf_obj", oldClass(splitf_obj))
  }
  N = names(splitp_obj$decomp)
  splitp_obj$decomp[N == "dir"] <- splitf_obj["dir"]
  splitp_obj$decomp[N == "parent"] <- splitf_obj["parent"]
  splitp_obj$decomp[N == "short"] <- splitf_obj["short"]
  splitp_obj$decomp[N == "ext"] <- splitf_obj["ext"]
  splitp_obj$decomp[N == "object"] <- object
  splitp_obj$decomp[N == "channel"] <- channel
  return(paste0(splitp_obj$decomp, collapse=""))
}

#' @title Date Converter
#' @description Formats date.
#' @param x a string.
#' @param iso whether returned value should be ISO8601. Default is \code{FALSE}.
#' @param tryFormats format to try. See \code{\link[base]{format.Date}}. Default is \code{c("\%m/\%d/\%Y \%I:\%M:\%S \%p", "\%d/\%m/\%Y \%H:\%M:\%S", "\%Y:\%m:\%d \%H:\%M:\%S", "\%Y-\%m-\%dT\%H:\%M:\%SZ")}.
#' @param tz time zone name. See \code{\link[base]{format.Date}}. Default is "UTC".
#' @return formatted date
#' @keywords internal
formatdate <- function(x, iso = FALSE, tryFormats = c("%m/%d/%Y %I:%M:%S %p", "%d/%m/%Y %H:%M:%S", "%Y:%m:%d %H:%M:%S", "%Y-%m-%dT%H:%M:%SZ"), tz="UTC") {
  xx = try(as.POSIXlt(x, tryFormats = tryFormats, tz = tz), silent = TRUE)
  if(inherits(xx, "try-error")) return(x)
  if(iso) {
    xxx = try(format.POSIXlt(xx, format = "%Y-%m-%dT%H:%M:%SZ"), silent = TRUE)
  } else {
    xxx = try(format.POSIXlt(xx, format = "%Y:%m:%d %H:%M:%S"), silent = TRUE)
  }
  if(inherits(xxx, "try-error")) return(x)
  return(xxx)
}

#' @title First Letter Only Capitalization
#' @description
#' Helper to capitalize the first letter of strings and leave the rest to lower case
#' @param text a string
#' @keywords internal
toCapFirstOnly <- function(text) {
  x = as.character(text)
  paste0(toupper(substr(x,1,1)), tolower(substr(x,2,nchar(x))))
}

#' @title Color Mapping
#' @name map_color
#' @description
#' Converts IDEAS/INSPIRE colors toR and inversely
#' @param color a character vector Default is missing.
#' @param toR whether to convert color toR or back. Default is TRUE.
#' @return a character vector.
#' @keywords internal
map_color <- function(color, toR = TRUE) {
  set1 = c("Teal", "Green", "Lime", "Control")
  set2 = c("Cyan4", "Green4", "Chartreuse", "Gray81")
  if(toR) {
    foo = color %in% set1
    if(any(foo)) {
      bar = sapply(color[foo], FUN = function(x) which(set1 %in% x))
      color[foo] <- set2[bar]
    }
  } else {
    foo = color %in% set2
    if(any(foo)) {
      bar = sapply(color[foo], FUN = function(x) which(set2 %in% x))
      color[foo] <- set1[bar]
    }
  }
  return(color)
}

#' @title Style Mapping
#' @name map_style
#' @description
#' Converts IDEAS/INSPIRE style toR and inversely
#' @param style a pch (converted to integer) or a character vector. Default is missing.
#' @param toR whether to convert color toR or back. Default is FALSE.
#' @return an integer vector when toR is TRUE or a character vector.
#' @keywords internal
map_style <- function(style, toR = FALSE) {
  set1 = c("Simple Dot", "Cross", "Plus", 
                       "Empty Circle", "Empty Diamond", "Empty Square", 
                       "Empty Triangle", "Solid Diamond", "Solid Square", 
                       "Solid Triangle")
  set2 = c(20, 4, 3, 1, 5, 0, 2, 18, 15, 17)
  if(toR) {
    foo = style %in% set1
    if(any(foo)) {
      bar = sapply(style[foo], FUN = function(x) which(set1 %in% x))
      style[foo] <- set2[bar]
    }
    if(!all(style %in% set2)) stop("not supported 'style'")
    style = as.integer(style)
  } else {
    style_ = suppressWarnings(as.integer(style))
    foo = style_ %in% set2
    if(any(foo)) {
      bar = sapply(style_[foo], FUN = function(x) which(set2 %in% x))
      style[foo] <- set1[bar]
    }
    if(!all(style %in% set1)) stop("not supported 'style'")
  }
  return(style)
}

#' @title Random Name Generator
#' @name random_name
#' @description
#' Generates random name
#' @param n number of characters of the desired return name. Default is \code{10}.
#' @param ALPHA upper case letters. Default is \code{LETTERS}.
#' @param alpha lower case letters. Default is \code{letters}.
#' @param num integer to use. Default is \code{0L:9L} 
#' @param special characters. Default is \code{c("#", "@@", "?", "!", "&", "\%", "$")}.
#' @param forbidden forbidden character vector. Default is \code{character()}.
#' @param prefix string to be prepended. Default is \code{""}.
#' @details 'forbidden' should not encompass all possible returned value otherwise the function will never end.
#' @return a character string.
#' @keywords internal
random_name <- function(n = 10, ALPHA = LETTERS, alpha = letters, num = 0L:9L, special = c("#", "@", "?", "!", "&", "%", "$"), forbidden = character(), prefix = "") {
  if(length(ALPHA)!=0) assert(ALPHA, alw = LETTERS)
  if(length(alpha)!=0) assert(alpha, alw = letters)
  if(length(num)!=0) assert(num, cla="integer", alw = 0L:9L)
  if(length(forbidden)!=0) forbidden = unique(forbidden); assert(forbidden, typ = "character")
  assert(prefix, len = 1, typ = "character")
  id = paste0(prefix, paste0(sample(x = c(ALPHA, alpha, num, special), size = n, replace = TRUE), collapse = ""))
  while(id %in% forbidden) {
    id = paste0(prefix, paste0(sample(x = c(ALPHA, alpha, num, special), size = n, replace = TRUE), collapse = ""))
  }
  return(id)
}

#' @title Alternative Names Generator
#' @name gen_altnames
#' @description
#' Generates unique non matching alternative names
#' @param x a character vector.
#' @param n number of characters of the desired returned name. Default is 10.
#' @param forbidden forbidden character vector. Default is character().
#' @param random_seed a list of elements to pass to \link[base]{set.seed} or a single value, interpreted as an integer, or NULL.
#' Default is list(seed = 0xFC, "Mersenne-Twister", "Inversion", "Rounding").
#' Note that NA_integer_ or list(seed = NA_integer_) can be used to not call \link[base]{set.seed} at all.
#' Note also that the default is chosen because it is compatible with old R version.
#' @details 'forbidden' should not encompass all possible returned value otherwise the function will never end.
#' @return a character vector.
#' @keywords internal.
gen_altnames <- function(x, n = 10, forbidden = character(),
                         random_seed = list(seed = 0xFC, "Mersenne-Twister", "Inversion", "Rounding")) {
  if(n < 5) stop("can't generate altnames with n < 5")
  x = as.character(x); assert(x, typ = "character")
  forbidden = unique(as.character(forbidden)); assert(forbidden, typ = "character")
  xx = unique(x)
  xx = unname(xx[order(nchar(xx), decreasing = TRUE)])
  pat = setdiff(c(xx, forbidden),"")
  spat = pat[nchar(pat) <= n]
  let = setdiff(letters, pat)
  LET = setdiff(LETTERS, pat)
  NUM = setdiff(0:9, suppressWarnings(na.omit(as.integer(pat))))
  SEED = fetch_seed(random_seed)
  # the suppressWarning is here to silence the "Rounding" stuff because gen_altnames function 
  # is used to generates alternative names internally for splitting definitions or else where
  # there is no need be statistically relevant
  suppressWarnings(with_seed({
    ans = c()
    for(i in seq_along(x)) {
      foo = paste0(sample(x = c(LET, let, NUM), size = n, replace = TRUE), collapse = "")
      while(cpp_mpfmatch(x = foo, pattern = spat)) foo = paste0(sample(x = c(LET, let, NUM), size = n, replace = TRUE), collapse = "")
      ans <- c(ans, foo)
    }
    ans
  }, SEED$seed, SEED$kind, SEED$normal.kind, SEED$sample.kind))
}

#' @title Numeric to String Formatting
#' @name num_to_string
#' @description
#' Formats numeric to string used for features, images, ... values conversion when exporting to xml.
#' @param x a numeric vector.
#' @param precision number of significant decimal digits to keep. Default is 22.
#' @return a string vector.
#' @keywords internal
num_to_string <- function(x, precision = 22) {
  old <- options("scipen")
  on.exit(options(old))
  options("scipen" = 18)
  xx = toupper(as.character(round(x, precision)))
  xx[is.na(x)] <- "NaN"
  return(xx)
}

#' @title Next Component Prediction
#' @description
#' Helper to define next allowed component in a boolean vector.
#' @param x a string, current component. Default is character(0).
#' @param count an integer, representing current number of opened/closed bracket.
#' @param obj_alias a character vector, alias used for name(s). Default is character(0).
#' @return a vector of next allowed components.
#' @keywords internal
next_bool = function(x = character(0), count = 0L, obj_alias = character(0)) {
  if(length(x) == 0) return(c("Not", "(", obj_alias))
  return(switch(x,
                "And" = {
                  c("Not", "(", obj_alias)
                },
                "Or" = {
                  c("Not", "(", obj_alias)
                },
                "Not" = {
                  c("(", obj_alias)
                },
                "(" = {
                  c("Not", "(", obj_alias)
                },
                ")" = {
                  tmp = c("And", "Or")
                  if(count > 0) tmp = c(tmp, ")")
                  tmp
                },
                {  
                  tmp = c("And", "Or")
                  if(count > 0) tmp = c(tmp, ")")
                  tmp
                }))
}

#' @title Boolean Expression Validation
#' @description
#' Helper to check if a boolean vector is valid.
#' @param x a string vector representing the boolean expression to be validated. Default is ""
#' @param all_names a character vector of scalars which are allowed to be part of the the boolean expression.
#' @return x is returned if no exception is raised during validation process.
#' @keywords internal
validate_bool = function(x = "", all_names = "") {
  operators = c("And", "Or", "Not", "(", ")")
  all_names = setdiff(all_names, c(operators, ""))
  if(length(x) == 0) stop("object definition is not possible: empty")
  if(!any(all_names %in% x)) stop("object definition is not possible: no match found")
  obj_alias = gen_altnames("foo", forbidden = c("", operators, all_names, x))
  count = 0L
  alw = c("Not", "(", obj_alias)
  lapply(1:length(x), FUN = function(i) {
    if(x[i] == "(") count <<- count + 1L
    if(x[i] == ")") count <<- count - 1L
    if(x[i] %in% alw) {
      alw <<- next_bool(x[i], count, obj_alias)
      return(NULL)
    } else {
      if((obj_alias %in% alw) && (x[i] %in% all_names)) {
        alw <<- next_bool(x[i], count, obj_alias)
        return(NULL)
      }
    }
    stop("object definition is not possible: '",x[i],"' is not allowed at position [",i,"] in\n",x)
  })
  if(count != 0) stop("object definition is not possible: invalid number of opened bracket")
  if(x[length(x)] %in% c("And", "Or", "Not", "(")) stop("object definition is not possible: it should not end with '",x[length(x)],"'")
  return(x)
}
