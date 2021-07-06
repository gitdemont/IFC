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

#' @title String Truncation
#' @description Truncs character strings
#' @param x a string
#' @param n desired length
#' @details x will be truncated according to 'n' parameter. If x is longer than n '...' are appended.
#' @keywords internal
trunc_string = function(x, n=22) {
  x = as.character(x)
  L = nchar(x)
  if(L > n) return(paste0(substr(x, 1, n),"..."))
  return(x)
}

#' @title File Extenstion Removal
#' @description Remove file extension from file path
#' @param x a file path
#' @details file extension will be removed
#' @keywords internal
remove_ext <- function(x) {
  x = as.character(x)
  ext = getFileExt(x)
  gsub(paste0("\\.",ext,"$"), "", x)
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

#' @title String Decomposition with Operators
#' @description
#' Helper that will split population definition into chunks of names and operators.
#' @param definition population definition to be splitted
#' @param all_names the names of all allowed populations
#' @param operators operators used. Default is c("And", "Or", "Not", "(", ")").
#' @keywords internal
splitn <- function(definition, all_names, operators = c("And", "Or", "Not", "(", ")")) {
  assert(definition, len=1, typ="character")
  assert(all_names, typ="character")
  assert(operators, typ="character")
  
  # we create a mapping between all_names + operators and random names
  # we also ensure that random names will bot contain any specials and 
  # will not match with themselves nor with names or operators to substitute
  # we also order to_substitute by number of character (decreasing) to
  # be sure that longer words will be substitute first
  to_substitute = c(all_names, operators)
  if((length(to_substitute) == 0) || (definition == "")) return(definition)
  to_substitute = to_substitute[order(nchar(to_substitute), decreasing = TRUE)]
  replace_with = c()
  for(i in 1:length(to_substitute)) { 
    replace_with = c(replace_with, random_name(n=10, special = NULL, forbidden = c(replace_with, to_substitute)))
  }
  
  # we substitute names and operators with random names
  ans = definition
  for(i in 1:length(to_substitute)) { 
    ans = gsub(pattern = to_substitute[i], replacement = replace_with[i], x = ans, fixed = TRUE)
  }
  
  # we can now split the definition with "|" since random names we use
  # do not contain this special character
  ans = strsplit(ans, split = "|", fixed = TRUE)[[1]]
  
  # finally, can replace random names with their corresponding initial names
  # from to_substitute (i.e. all_names + operators)
  ans = sapply(ans, USE.NAMES = FALSE, FUN = function(x) {
    foo = to_substitute[x == replace_with]
    if(length(foo) == 0) stop("definition contains unexpected name")
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

#' @title Color Mapping
#' @name map_color
#' @description
#' Converts IDEAS/INSPIRE colors toR and inversely
#' @param color a character vector Default is missing.
#' @param toR whether to convert color toR or back. Default is TRUE.
#' @return a character vector
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
#' @param n number of characters of the desired return name. Default is 10.
#' @param ALPHA upper case letters. Default is LETTERS.
#' @param alpha lower case letters. Default is letters.
#' @param num integer to use. Default is 0:9 
#' @param special characters. Default is c("#", "@@", "?", "!", "&", "\%", "$").
#' @param forbidden forbidden character strings. Default is character().
#' @return a character string
#' @keywords internal
random_name <- function(n = 10, ALPHA = LETTERS, alpha = letters, num = 0L:9L, special = c("#", "@", "?", "!", "&", "%", "$"), forbidden = character()) {
  if(length(ALPHA)!=0) assert(ALPHA, alw = LETTERS)
  if(length(alpha)!=0) assert(alpha, alw = letters)
  if(length(num)!=0) assert(num, cla="integer", alw = 0L:9L)
  id = paste0(sample(x = c(ALPHA, alpha, num, special), size = n), collapse = "")
  while(id %in% forbidden) { id <- random_name(n = n, ALPHA = ALPHA, alpha = alpha, num = num, special = special, forbidden = forbidden) }
  return(id)
}

#' @title Numeric to String Formatting
#' @name num_to_string
#' @description
#' Formats numeric to string used for features, images, ... values conversion when exporting to xml.
#' @param x a numeric vector.
#' @param precision number of significant decimal digits to keep when abs(x) < 1. Default is 15.
#' @return a string vector.
#' @keywords internal
num_to_string <- function(x, precision = 15) {
  return(formatC(x, digits = precision, width = -1, drop0trailing = TRUE))
  # return(cpp_num_to_string(x, precision))
}

#' @title Next Component Prediction
#' @description
#' Helper to define next allowed component in a boolean vector.
#' @param x a string, current component.
#' @param count an integer, representing current number of opened/closed bracket.
#' @return a vector of next allowed components.
#' @keywords internal
next_bool = function(x = "", count = 0L) {
  return(switch(x,
                "And" = {
                  c("Not", "(", "obj")
                },
                "Or" = {
                  c("Not", "(", "obj")
                },
                "Not" = {
                  c("(", "obj")
                },
                "(" = {
                  c("Not", "(", "obj")
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
#' @param x a string vector representing the boolean expression to be validated.
#' @param all_names a character vector of scalars which are allowed to be part of the the boolean expression.
#' @return x is returned if no exception is raised during validation process.
#' @keywords internal
validate_bool = function(x = "", all_names = "") {
  operators = c("And", "Or", "Not", "(", ")")
  all_names = setdiff(all_names, c(operators, ""))
  if(!any(all_names %in% x)) stop("object definition is not possible: no match found")
  count = 0L
  alw = c("Not", "(", "obj")
  lapply(1:length(x), FUN = function(i) {
    if(x[i] == "(") count <<- count + 1L
    if(x[i] == ")") count <<- count - 1L
    if(x[i] %in% alw) {
      alw <<- next_bool(x[i], count)
      return(NULL)
    } else {
      if(("obj" %in% alw) && (x[i] %in% all_names)) {
        alw <<- next_bool(x[i], count)
        return(NULL)
      }
    }
    stop("object definition is not possible: '",x[i],"' is not allowed at position [",i,"] in\n",x)
  })
  if(count != 0) stop("object definition is not possible: invalid number of opened bracket")
  if(x[length(x)] %in% c("And", "Or", "Not", "(")) stop("object definition is not possible: it should not end with '",x[length(x)],"'")
  return(x)
}
