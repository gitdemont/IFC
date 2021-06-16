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

#' @title Parse Transformation Information
#' @name parseTrans
#' @description Helper to parse transformation parameter.
#' @param string character string describing transformation used and its hyperparameter.
#' @details -If string is "P" no transformation will be applied.\cr
#' -If string is of length 1 and coercible to a numeric, smoothLinLog will be applied with hyper parameter set with as.numeric(string) value.\cr
#' -Otherwise, string will be split with "|", 1st element will be considered as function to call and all other elements will be passed after being coerced to numeric to this function in the order they are provided with the exception of 1st one; if coercion results in NA, the argument will be skipped.
#' @examples 
#' ## this will use smoothLinLog, skip argument hyper and fill base with 2
#' parseTrans("smoothLinLog||2")
#' @return a list with 2 members what (=the transformation function) and args (=the parameters to pass to this function) 
#' @keywords internal
parseTrans <- function(string) {
  if(length(string) == 0) return(list(what = "return", args = list()))
  foo = strsplit(string, split = "|", fixed = TRUE)[[1]]
  if(length(foo) == 1 && foo == "P") return(list(what = "return", args = list()))
  if(length(foo) != 1 || is.na(suppressWarnings(as.numeric(foo[1])))) {
    fun = foo[1]
    fun_args <- formalArgs(fun)
    args <- suppressWarnings(as.numeric(foo)) # eliminates 1st arg (e.g 'x')
    names(args) <- fun_args[seq_along(foo)]
    args <- args[!is.na(args)]
    return(list(what = fun, args = as.list(args)))
  } else {
    return(list(what = "smoothLinLog", args = list(hyper = as.numeric(foo[1]))))
  }
}

# version to include ...
# parseTrans <- function(string) {
#   foo = strsplit(string, split = "|", fixed = TRUE)[[1]]
#   if(length(foo) == 1 && foo == "P") return(list(what = "return", args = list()))
#   if(length(foo) != 1 || is.na(suppressWarnings(as.numeric(foo[1])))) {
#     fun = foo[1]
#     fun_args <- formals(fun)
#     args <- lapply(foo, FUN = function(x) suppressWarnings(as.numeric(x)))
#     names(args) <- names(fun_args[seq_along(foo)])
#     for(i in na.omit(names(args[is.na(unlist(args))]))) {
#       try({args[[i]] <- eval(parse(text=fun_args[[i]])) }, silent = TRUE) # try is used to not include args that evaluates other args 
#     }
#     if(anyNA(names(args)) && "..." %in% names(args)) { args[["..."]] <- c(args[["..."]], na.omit(unname(unlist(args[is.na(names(args))])))) }
#     args <- args[!is.na(unlist(args))]
#     args <- args[!is.na(names(args))]
#     return(list(what = fun, args = args))
#   } else {
#     return(list(what = "smoothLinLog", args = list(hyper = as.numeric(foo[1]))))
#   }
# }

#' @title Apply Transformation
#' @name applyTrans
#' @description Helper apply transformation.
#' @param x a numeric vector.
#' @param trans the object returned by parseTrans().
#' @param inverse whether or not to apply the inverse transformation. Default is FALSE.
#' @details for the moment, in addition to no transformation, only smoothLinLog and asinh are supported.
#' @examples
#' x <- 10^(1:10)
#' ## parse the transformation instruction
#' trans <- parseTrans("smoothLinLog|100|4")
#' ## transform x
#' xx <- applyTrans(x, trans)
#' ## inverse transform the transformed values
#' x_back <- applyTrans(xx, trans, inverse = TRUE)
#' stopifnot(all.equal(x, x_back))
#' @return the transformation of the input
#' @keywords internal
applyTrans <- function(x, trans, inverse = FALSE) {
  fun = trans$what
  if(inverse) fun = switch(fun,
                           "return" = "return" ,
                           "smoothLinLog" = "inv_smoothLinLog",
                           "smoothAsinh" = "inv_smoothAsinh",
                           "asinh" = "sinh",
                           stop("can't find inverse transformation for: '",fun,"'" ))
  do.call(what = fun, args = c(list(x), trans$args))
}

#' @title Smooth Asinh Transformation
#' @description Transforms values in asinh
#' @keywords internal
smoothAsinh <- function(x, hyper = 1000) {
  return(asinh(x/hyper))
}

#' @title Inverse Smooth Asinh Transformation
#' @description Gets values back just to their original values before applying smoothAsinh
#' @keywords internal
inv_smoothAsinh <- function(x, hyper = 1000) {
  return(sinh(x)*hyper)
}
