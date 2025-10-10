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

#' @title Add Region to IFC_data Object
#' @description
#' Adds regions to an already existing `IFC_data` object.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param regions a list of region(s) to add to obj. Each element of this list will be coerced by \code{\link{buildRegion}}.
#' @param create_pops whether population(s) corresponding to \code{'regions'} should be also created. Default is \code{FALSE}.
#' If \code{TRUE}, at least \code{'base'}, \code{'fx'} and, eventually, \code{'fy'} should be members of sub elements of \code{'regions'}.
#' Other members, see \code{\link{buildPopulation}}, are not mandatory but can also be provided.
#' @details A warning will be thrown if a provided region is duplicated or already existing in 'obj'.\cr
#' In such a case this region, and resulting population when \code{'create_pops'} is \code{TRUE}, will not be added to 'obj'.\cr
#' If any input region is not well defined, or resulting population when \code{'create_pops'} is \code{TRUE}, and can't be created then an error will occur.
#' @param ... Other arguments to be passed.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#'   ## copy 1st region found in daf
#'   reg <- daf$regions[[1]]
#'   if(length(reg) != 0) {
#'     reg_copy <- reg
#'     ## modify region label and x boundaries
#'     reg_copy$label <- paste0(reg_copy$label,"_copy")
#'     reg_copy$x <- reg_copy$x*0.9
#'     ## create new object with this new region
#'     dafnew <- data_add_regions(obj = daf, regions = list(reg_copy))
#'   }
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return an IFC_data object with regions added.
#' @export
data_add_regions <- function(obj, regions, create_pops = FALSE, ...) {
  assert(obj, cla = "IFC_data")
  
  # extract regions names
  names(regions) = sapply(regions, FUN=function(x) x$label)
  
  # removes duplicated inputs
  tmp = duplicated(names(regions))
  if(any(tmp)) {
    warning("duplicated regions automatically removed:\n\t-", paste0(unique(names(regions)[tmp]),collapse="\n\t-"), immediate. = TRUE, call. = FALSE)
    regions = regions[!tmp]
  }
  
  # removes already defined regions
  exported_regions = sapply(regions, FUN=function(reg) {
    if(reg$label%in%names(obj$regions)) {
      warning(paste0(reg$label, "\nnot exported: trying to export an already defined region"), immediate. = TRUE, call. = FALSE)
      return(FALSE)
    }
    return(TRUE)
  })
  regions = regions[exported_regions]
  
  # build pops if create_pops is TRUE
  pops = list()
  if(create_pops) {
    pops = lapply(regions, FUN = function(r) {
      tryCatch({
        args = list(type = "G",
                    color = r$color,
                    lightModeColor = r$lightcolor,
                    region = r$label)
        if(any("base" %in% names(r))) {
          args$base = r$base
          args$name = ifelse(identical(r$base, "All"), r$label, paste(r$label, r$base, sep = " & "))
        }
        if(any("name" %in% names(r))) args$name = r$name
        if(any(arg$name %in% names(obj$pops))) stop("population ['",arg$name,"' is already defined in 'obj'")
        if(any("style" %in% names(r))) args$style = r$style
        if(!any(r$fx %in% names(obj$features))) stop("'fx' is not a feature of 'obj'")
        if(any("fx" %in% names(r))) args$fx = r$fx
        if(!identical(r$type, "line")) {
          if(!any(r$fy %in% names(obj$features))) stop("'fy' is not a feature of 'obj'")
          args$fy = r$fy
        }
        do.call(buildPopulation, args)
      }, error = function(e) {
        stop("can't create population for region['",r$label,"']:\n", 
                e$message, call. = FALSE)
      })
    })
    names(pops) = sapply(pops, FUN = function(p) p$name)
    tmp = duplicated(c(names(pops), names(obj$pops))) # throw error if creating pop results in already existing pop
    if(any(tmp)) stop("creating populations will lead to duplicated names (provide another 'name' to solve this):\n\t-", paste0(unique(c(names(pops), names(obj$pops))[tmp]),collapse="\n\t-"))
  }
  
  # try to coerce regions inputs
  regions = lapply(regions, keep_attributes, what=buildRegion)
  names(regions) = sapply(regions, FUN=function(x) x$label)
  tmp = duplicated(c(names(regions), names(obj$regions))) # normally, this should not happen
  if(any(tmp)) stop("creating regions will lead to duplicated names:\n\t-", paste0(unique(c(names(regions), names(obj$regions))[tmp]),collapse="\n\t-"))
  
  # change colors to R compatible
  for(i in 1:length(regions)) {
    regions[[i]]$color = map_color(regions[[i]]$color)
    regions[[i]]$lightcolor = map_color(regions[[i]]$lightcolor)
  }
  
  K = class(obj$regions)
  obj$regions = c(obj$regions, regions)
  class(obj$regions) = K
  
  if(length(pops) != 0) return(data_add_pops(obj, pops, ...))
  return(obj)
}
