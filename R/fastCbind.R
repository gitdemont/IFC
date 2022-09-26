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
  old <- options("stringsAsFactors")
  on.exit(options(old))
  options("stringsAsFactors" = FALSE)
  d1 = dim(obj1); l1 = length(d1); t1 = typeof(obj1)
  d2 = dim(obj2); l2 = length(d2); t2 = typeof(obj2)
  if(typeof(obj1) == "list") {
    N1 = names(obj1)
  } else {
    N1 = dimnames(obj1)[[2]]
  }
  if(typeof(obj2) == "list") {
    N2 = names(obj2)
  } else {
    N2 = dimnames(obj2)[[2]]
  }
  if((l1 == 2) || (l2 == 2)) {
    if((l1 == 0) && (typeof(obj1) == "list")) obj1 = do.call(what = ifelse(t2 == "list", "data.frame", "cbind"), args = obj1)
    if((l2 == 0) && (typeof(obj2) == "list")) obj2 = do.call(what = ifelse(t1 == "list", "data.frame", "cbind"), args = obj2)
  }
  if(l1 != 2) {
    if(typeof(obj1) == "list") {
      obj1 = do.call(what = ifelse(t2 == "list", "data.frame", "cbind"), args = obj1)
    } else {
      obj1 = do.call(what = ifelse(t2 == "list", "data.frame", "cbind"), args = list(obj1))
    }
  }
  if(l2 != 2) {
    if(typeof(obj2) == "list") {
      obj2 = do.call(what = ifelse(t1 == "list", "data.frame", "cbind"), args = obj2)
    } else {
      obj2 = do.call(what = ifelse(t1 == "list", "data.frame", "cbind"), args = list(obj2))
    }
  }
  d1 = dim(obj1); l1 = length(d1)
  d2 = dim(obj2); l2 = length(d2)
  if((t1 != "list") && (t2 != "list")) {
    ans = matrix(ncol = 0, nrow = 0)
  } else {
    ans = data.frame() 
  }
  if((l1 == 2) && (l2 == 2)) {
    if(d1[1] == 0) ans = cbind(obj2, deparse.level = 0)
    if(d2[1] == 0) ans = cbind(obj1, deparse.level = 0)
    if((d1[1] != 0) && (d2[1] != 0)) ans = cbind(obj1, obj2, deparse.level = 0)
    if((d1[1] == 0) && (d2[1] == 0)) ans = cbind(obj1, obj2, deparse.level = 0)
    s1 = ncol(obj1)
    s2 = ncol(obj2)
  } else {
    if(l1 == 2) {
      ans = cbind(obj1, deparse.level = 0)
      s1 = ncol(ans); s2 = 0
    }
    if(l2 == 2) {
      ans = cbind(obj2, deparse.level = 0)
      s2 = ncol(ans); s1 = 0
    }
  }
  if(length(N1) == 0) replicate(s1, N1 <<- c(N1, random_name(special = NULL, forbidden = as.character(c(N1, N2)))))
  if(length(N2) == 0) replicate(s2, N2 <<- c(N2, random_name(special = NULL, forbidden = as.character(c(N1, N2)))))
  while(length(N1 == "") != s1) {
    tmp = which(N1=="")[1]
    N1[ifelse(is.na(tmp),1,tmp)] <- random_name(special = NULL, forbidden = as.character(c(N1, N2)))
  }
  while(length(N2 == "") != s2) {
    tmp = which(N2=="")[1]
    N2[ifelse(is.na(tmp),1,tmp)] <- random_name(special = NULL, forbidden = as.character(c(N1, N2)))
  }
  N = NULL
  if(add_id) {
    ans = cbind(seq_along(integer(nrow(ans))), ans, deparse.level = 0)
    N = random_name(special = NULL, forbidden = as.character(c(N1, N2)))
  }
  all_names = c(N, N1, N2)
  if(anyDuplicated(all_names)) stop("names should be unique")
  colnames(ans) = all_names
  return(ans)
}
