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

#' @title Seed Instructions Extraction
#' @description
#' Retrieve seed instructions from input
#' @param seed a list of elements to pass to \link[base]{set.seed} or a single value, interpreted as an integer, or NULL. NA_integer_ or list(seed = NA_integer_) can be used to prevent 'seed' argument from being passed to \link[base]{set.seed}. Default is NA_integer_.
#' @keywords internal
fetch_seed <- function(seed = NA_integer_) {
  ans = list(seed = NA_integer_, kind = NULL, normal.kind = NULL, sample.kind = NULL)
  old <- options("warn")
  on.exit(options(old))
  options("warn" = 2)
  if(typeof(seed) == "list") {
    Ns = names(seed)
    Na = names(ans) 
    if(anyDuplicated(Ns)) stop("found duplicated names in 'seed'")
    a = c()
    b = c()
    for(ii in 1:4) {            # pass by name
      tmp = which(Ns == Na[ii])
      if(length(tmp) != 0) {  
        a = c(a, tmp)
        b = ii
        if(is.null(seed[[tmp]])) {
          ans[ii] <- list(NULL)
        } else{
          ans[[ii]] <- seed[[tmp]]
        }
      }
    }
    seed = seed[setdiff(seq_along(seed), a)]
    for(ii in setdiff(1:4, b)) { # pass by index
      if(length(seed) != 0) {
        if(is.null(seed[[1]])) {
          ans[ii] <- list(NULL)
        } else {
          ans[[ii]] <- seed[[1]]  
        }
        seed = seed[-1]
      }
    }
  } else {
    if(is.null(seed)) {
      ans["seed"] = list(NULL)
    } else {
      ans$seed = seed
    }
  }
  v = ans$seed
  if(length(v) != 0) ans$seed = try(as.integer(v), silent = TRUE)
  if((length(v) != 0) && (length(v) > 1 ||
     is.nan(v) ||
     (!is.na(v) && inherits(ans$seed, "try-error")))) stop("'seed' should be either NULL, NA or a single integer")
  return(ans)
}

#' @title With Seed Evaluation
#' @description
#' Evaluates expression with a seed and resets to initial seed state on exit
#' @param expr expression to evaluate.
#' @param seed a single value, interpreted as an integer, or NULL, with the exception that NA can be provided to prevent passing 'seed' argument. Default is NA_integer_.
#' @param kind character or NULL. If kind is a character string, set R's RNG to the kind desired. Use "default" to return to the R default.
#' @param normal.kind	character string or NULL. If it is a character string, set the method of Normal generation. Use "default" to return to the R default. NULL makes no change.
#' @param sample.kind	character string or NULL. If it is a character string, set the method of discrete uniform generation (used in sample, for instance). Use "default" to return to the R default. NULL makes no change.
#' @details see ‘Details’, from  \link[base]{set.seed}, with the exception of 'seed'
#' @keywords internal
with_seed <- function(expr, seed = NA_integer_, kind = NULL, normal.kind = NULL, sample.kind = NULL) {
  cur_kind = RNGkind()
  e = globalenv()
  cur_seed <- e$.Random.seed
  on.exit(suspendInterrupts({
    do.call(args = list(kind = cur_kind[1], normal.kind = cur_kind[2], sample.kind = cur_kind[3]), what = RNGkind)
    if (is.null(cur_seed)) {
      rm(".Random.seed", envir = e, inherits = FALSE)
    } else {
      assign(".Random.seed", value = cur_seed, envir = e, inherits = FALSE)
    }
  }))
  RNGkind(kind = kind, normal.kind = normal.kind, sample.kind = sample.kind)
  if(!is.na(seed)) set.seed(seed = seed)
  force(expr)
}
