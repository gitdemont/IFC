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

#' @title Modify Populations within IFC_data Object
#' @description
#' Modify populations in an already existing `IFC_data` object.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param regions a list of region(s) to modify in 'obj'. Each element of this list will be coerced by \code{\link{buildRegion}}.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @details regions names should be present in names(obj$regions), otherwise an error will be raised.\cr
#' Note that If you want to rename regions, you should do it by changing 'label' member,
#' e.g. regions[[1]]$label <- "bar" while names(regions[[1]]) is "foo" and "foo" is part of names(obj$regions).
#' @param ... Other arguments to be passed.
#' @return an IFC_data object with regions modified.
#' @keywords internal
data_modify_regions <- function(obj, regions, display_progress = TRUE, ...) {
  dots = list(...)
  assert(obj, cla = "IFC_data")
  mutation = names(regions)
  if(!all(mutation %in% names(obj$regions))) stop("can't find regions to modify in 'obj'", call. = FALSE)
  R = lapply(regions, FUN = function(x) {
    reg = do.call(what = buildRegion, args = x)
    reg$color = map_color(reg$color)
    reg$lightcolor = map_color(reg$lightcolor)
    reg
  })
  names(R) = sapply(R, FUN = function(x) x$label)
  tmp = duplicated(names(R))
  if(any(tmp)) stop(paste0("duplicated regions found: ", unique(names(regions)[tmp])), call. = FALSE)
  names(mutation) = names(R)
  type1 = sapply(obj$regions[mutation], FUN = function(r) r$type)
  type2 = sapply(R, FUN = function(r) r$type)
  tmp = type1 == type2
  if(!all(tmp)) stop(paste0("'type' modification is not allowed:\n\t- ",
                            paste0(paste(paste0(mutation[!tmp], " [", type1[!tmp], "]"),
                                         paste0(names(mutation)[!tmp], " [", type2[!tmp], "]"), sep = " -> "),
                                   collapse = "\n\t- ")), call. = FALSE)
  
  tmp = names(mutation) %in% names(obj$region[setdiff(names(obj$region), mutation)])
  if(any(tmp)) stop("trying to rename region",ifelse(sum(tmp) > 1, "s", ""),
                    " with already existing name",ifelse(sum(tmp) > 1, "s", ""),":\n\t-",
                    paste0(names(mutation)[tmp], collapse="\n\t-"))
  
  # regions modification, it can be everything
  ans = obj
  N = names(ans$regions)
  K = class(ans$regions)
  ans$regions[mutation] = R
  names(ans$regions) = sapply(ans$regions, FUN = function(r) r$label)
  class(ans$regions) <- K
  
  # keep track of graphical population names that could be modified due to region renaming
  old_pop_names = c()
  new_pop_names = c()
  
  # now we can only modify names within pops and graphs
  if(length(mutation) != 0) {
    # populations modification
    for(i_pop in seq_along(ans$pops)) {
      found = mutation %in% ans$pops[[i_pop]]$region
      if(any(found)) {
        foo = sub(paste0(ans$pops[[i_pop]]$region, " & "), "", ans$pops[[i_pop]]$name, fixed = TRUE)
        if(foo != ans$pops[[i_pop]]$region) {
          if(names(mutation[found]) != mutation[found]) {
            old_pop_names <- c(old_pop_names, ans$pops[[i_pop]]$name)
            new_pop_names <- c(new_pop_names, paste(names(mutation[found]), foo, sep = " & "))
          }
        } else {
          if(names(mutation[found]) != mutation[found]) {
            old_pop_names <- c(old_pop_names, ans$pops[[i_pop]]$name)
            new_pop_names <- c(new_pop_names, names(mutation[found]))
          }
        }
        ans$pops[[i_pop]]$region <- names(mutation[found])
      }
    }
    
    # graphs modification
    for(i_graph in seq_along(ans$graphs)) {
      for(i_g in seq_along(ans$graphs[[i_graph]]$GraphRegion)) {
        found = ans$graphs[[i_graph]]$GraphRegion[[i_g]]$name == mutation
        if(any(found)) ans$graphs[[i_graph]]$GraphRegion[[i_g]]$name <- names(mutation[found])
      }
    }
  }
  
  # now we can change name of resulting populations if any
  if(length(old_pop_names) != 0) {
    new_pop_names_back = new_pop_names
    found = TRUE
    n = 1
    while(found) {
      tmp = new_pop_names %in% setdiff(names(ans$pops),old_pop_names) 
      found = any(tmp)
      if(found) new_pop_names[tmp] = paste0(new_pop_names_back[tmp], n)
      n = n + 1
    }
    ans = do.call(what = popsRename,
                  args = c(list(obj = ans,
                                old_names = old_pop_names,
                                new_names = new_pop_names),
                           dots[setdiff(names(dots), c("obj","old_names","new_pop_names"))]))
  }
  
  # since regions have been changed we need to recompute objects
  ans$pops <- do.call(what = popsCompute,
                      args = c(list(pops = ans$pops,
                                    regions = ans$regions,
                                    features = ans$features,
                                    display_progress = display_progress),
                               dots[setdiff(names(dots),c("pops","features"))]))
  
  # modify stats
  if(nrow(obj$stats)!=0) ans$stats = get_pops_stats(ans$pops, as.integer(obj$description$ID$objcount))
  return(ans)
}
