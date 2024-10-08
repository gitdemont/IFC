################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2020 Yohann Demont                                             #
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

#' @title Add Population to IFC_data Object
#' @description
#' Adds populations to an already existing `IFC_data` object.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param pops a list of population(s) to add to 'obj'. Each element of this list will be coerced by \code{\link{buildPopulation}}.
#' @param pnt_in_poly_algorithm algorithm used to determine if object belongs to a polygon region or not. Default is 1.\cr
#' Note that for the moment only 1(Trigonometry) is available.
#' @param pnt_in_poly_epsilon epsilon to determine if object belongs to a polygon region or not. It only applies when algorithm is 1. Default is 1e-12.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @details A warning will be thrown if a provided population is already existing in 'obj'.\cr
#' In such a case this population will not be added to 'obj'.\cr
#' If any input population is not well defined and can't be created then an error will occur.
#' @param ... Other arguments to be passed.
#' @source For pnt_in_poly_algorithm, Trigonometry, is an adaptation of Jeremy VanDerWal's code \url{https://github.com/jjvanderwal/SDMTools}
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#'   ## copy 1st population from existing daf
#'   pop <- daf$pops[[1]]
#'   if(length(pop) != 0) {
#'     pop_copy <- pop
#'     ## modify name, obj and type of copied population
#'     pop_copy$name <- paste0(pop_copy$name,"_copy")
#'     pop_copy$obj <- (which(pop_copy$obj)-1)[1]
#'     pop_copy$type <- "T"
#'     ## create new object with this new population
#'     dafnew <- data_add_pops(obj = daf, pops = list(pop_copy))
#'   }
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return an IFC_data object with pops added.
#' @export
data_add_pops <- function(obj, pops, pnt_in_poly_algorithm = 1, pnt_in_poly_epsilon = 1e-12, display_progress = TRUE, ...) {
  dots = list(...)
  assert(obj, cla = "IFC_data")
  obj_number = obj$description$ID$objcount
  
  # try to coerce inputs to compatible daf format
  pops = lapply(pops, keep_attributes, what=buildPopulation)
  names(pops) = sapply(pops, FUN=function(x) x$name)

  # removes duplicated inputs
  tmp = duplicated(names(pops))
  if(any(tmp)) {
    warning("duplicated pops automatically removed:\n\t-", paste0(names(pops)[tmp], collapse = "\n\t-"), immediate. = TRUE, call. = FALSE)
    pops = pops[!tmp]
  }
  
  exported_pops = sapply(1:length(pops), FUN=function(i_pop) {
    pop = pops[[i_pop]]
    if(pop$name%in%names(obj$pops)) {
      warning(pop$name, ", not added: trying to add an already defined population", immediate. = TRUE, call. = FALSE)
      return(FALSE)
    }
    if(pop$type=="T") {
      if(anyNA(pop$obj)) stop(pop$name, ", trying to add a tagged population containing NA/NaN")
      K = typeof(pop$obj)
      if(length(pop$obj)==0) {
        warning(pop$name, ", not exported: trying to add a tagged population of length = 0", immediate. = TRUE, call. = FALSE)
        return(FALSE)
      }
      if(K%in%"logical") {
        if(sum(pop$obj) ==0) {
          warning(pop$name, ", not exported: trying to add a tagged population of length = 0", immediate. = TRUE, call. = FALSE)
          return(FALSE)
        }
        if(obj_number != length(pop$obj)) stop(pop$name, ", trying to add a tagged population with more element(s) than total number of objects acquired")
      }
      if(K%in% c("double","integer")) {
        if((obj_number <= max(pop$obj)) | (min(pop$obj) < 0) | any(duplicated(pop$obj))) stop(pop$name, ", trying to add a tagged population with element(s) outside of objects acquired")
        log_obj = rep(FALSE, obj_number)
        log_obj[pop$obj + 1] = TRUE
        pops[[i_pop]]$obj <<- log_obj
      }
    }
    return(TRUE)
  })
  exported_pops = pops[exported_pops]
  names(exported_pops) = sapply(exported_pops, FUN=function(x) x$name)
  
  exported_pops = popsCompute(pops = c(obj$pops, exported_pops), 
                              regions = obj$regions, 
                              features = obj$features, 
                              pnt_in_poly_algorithm = pnt_in_poly_algorithm, 
                              pnt_in_poly_epsilon = pnt_in_poly_epsilon, 
                              display_progress = display_progress, 
                              title_progress = basename(obj$fileName), ...)

  obj$pops = exported_pops
  if(nrow(obj$stats)!=0) obj$stats = get_pops_stats(exported_pops, as.integer(obj$description$ID$objcount))
  return(obj)
}
