################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2022 Yohann Demont                                             #
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

#' @title Combine by Columns
#' @name fastCbind
#' @description
#' Helper to combine by columns
#' @param obj1 an object either a data.frame, a list or something coercible to numeric matrix.
#' @param obj2 an object either a data.frame, a list or something coercible to numeric matrix.
#' @param add_id a bool determining if 1st column of returned object should be given 1 to nrow integers
#' @details if obj1 or obj2 is a data.frame returned object will inherit class of obj1 or obj2 respectively in this order.
#' /!\ if obj1 or obj2 needs to be coerced take care that it can not handle factor during object coercion
#' @return a combined object
#' @keywords internal
fastCbind <- function(obj1, obj2, add_id = FALSE) {
  # fun to check and create columns names
  get_colnames <- function(obj1, obj2, add_id = FALSE) {
    N1 = names(obj1); N1 = N1[N1 != ""]
    N2 = names(obj2); N2 = N2[N2 != ""]
    s1 = ncol(obj1); if(length(s1) == 0) s1 = length(obj1)
    s2 = ncol(obj2); if(length(s2) == 0) s2 = length(obj2)
    while(length(N1) != s1) N1 = c(N1, random_name(special = NULL, forbidden = c(N1, N2)))
    while(length(N2) != s2) N2 = c(N2, random_name(special = NULL, forbidden = c(N1, N2))) 
    N = NULL
    if(add_id) N = random_name(special = NULL, forbidden = c(N1, N2))
    all_names = c(N, N1, N2)
    if(anyDuplicated(all_names)) stop("names should be unique")
    return(all_names)
  }
  if(inherits(obj1, "data.frame") && inherits(obj2, "data.frame")) {
    all_names = get_colnames(obj1, obj2, add_id)
    return(structure(cpp_fast_cbind_DF_DF(obj1, obj2, add_id = add_id), names = all_names, class = class(obj1)))
  }
  if(inherits(obj1, "data.frame")) {
    if(inherits(obj2, "list")) {
      ans = cpp_fast_cbind_DF_L(obj1, obj2, add_id = add_id)
      class(ans) = class(obj1)
      colnames(ans) = get_colnames(obj1, obj2, add_id)
      return(ans)
    }
    if(inherits(obj2, "matrix")) {
      ans = cpp_fast_cbind_DF_M(obj1, obj2, add_id = add_id)
      class(ans) = class(obj1)
      colnames(ans) = get_colnames(obj1, obj2, add_id)
      if(typeof(obj2) == "logical") lapply(colnames(obj2), FUN = function(x) {ans[, x] <<- as.logical(ans[, x])})
      if(typeof(obj2) == "integer") lapply(colnames(obj2), FUN = function(x) {ans[, x] <<- as.integer(ans[, x])})
      return(ans)
    } else {
      M2 = as.matrix(obj2)
      ans = cpp_fast_cbind_DF_M(obj1, M2, add_id = add_id)
      class(ans) = class(obj1)
      colnames(ans) = get_colnames(obj1, M2, add_id)
      if(typeof(M2) == "logical") lapply(colnames(M2), FUN = function(x) {ans[, x] <<- as.logical(ans[, x])})
      if(typeof(M2) == "integer") lapply(colnames(M2), FUN = function(x) {ans[, x] <<- as.integer(ans[, x])})
      return(ans) 
    }
  }
  if(inherits(obj2, "data.frame")) {
    if(inherits(obj2, "list")) {
      ans = cpp_fast_cbind_L_DF(obj1, obj2, add_id = add_id)
      class(ans) = class(obj2)
      colnames(ans) = get_colnames(obj1, obj2, add_id)
      return(ans)
    }
    if(inherits(obj1, "matrix")) {
      ans = cpp_fast_cbind_M_DF(obj1, obj2, add_id = add_id)
      class(ans) = class(obj2)
      colnames(ans) = get_colnames(obj1, obj2, add_id)
      if(typeof(obj1) == "logical") lapply(colnames(obj1), FUN = function(x) {ans[, x] <<- as.logical(ans[, x])})
      if(typeof(obj1) == "integer") lapply(colnames(obj1), FUN = function(x) {ans[, x] <<- as.integer(ans[, x])})
      return(ans)
    } else {
      M1 = as.matrix(obj1)
      ans = cpp_fast_cbind_M_DF(M1, obj2, add_id = add_id)
      class(ans) = class(obj2)
      colnames(ans) = get_colnames(M1, obj2, add_id)
      if(typeof(M1) == "logical") lapply(colnames(M1), FUN = function(x) {ans[, x] <<- as.logical(ans[, x])})
      if(typeof(M1) == "integer") lapply(colnames(M1), FUN = function(x) {ans[, x] <<- as.integer(ans[, x])})
      return(ans)
    }
  }
  M1 = as.matrix(obj1)
  M2 = as.matrix(obj2)
  all_names = get_colnames(M1, M2, add_id)
  return(structure(cpp_fast_cbind_M_M(M1, M2, add_id = add_id), colnames = all_names))
}
