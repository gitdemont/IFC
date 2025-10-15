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
    xx = keep_attributes(x, what = buildRegion)
    xx$color = map_color(xx$color)
    xx$lightcolor = map_color(xx$lightcolor)
    xx
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
  ans$regions[mutation] = lapply(seq_along(mutation), FUN = function(i) {
    x = obj$regions[mutation][[i]]
    for(k in setdiff(names(attributes(x)),c("names","dim","dimnames"))) {
      if(!identical(attr(x, k), attr(R[[i]], k))) {
        if(!is.null(attr(R[[i]], k))) warning("modification of attribute['",k,"'] is not allowed for ", x$label, call. = FALSE, immediate. = TRUE)
        attr(R[[i]], k) <- attr(x, k)
      }
    }
    for(k in setdiff(names(attributes(R[[i]])),c("names","dim","dimnames"))) { 
      if(is.null(attr(x, k)) && !is.null(attr(R[[i]], k))) {
        warning("creation of attribute['",k,"'] is not allowed for ", x$label, call. = FALSE, immediate. = TRUE) 
        attr(R[[i]], k) <- NULL
      }
    }
    R[[i]]
  })
  sync_names = sync_parse(ans$regions)
  sync_mut = structure(lapply(obj$regions[mutation], sync_name), names = mutation)
  for(s in unlist(setdiff(na.omit(sync_mut), ""))) {
    tmp = sapply(sync_names, identical, s)
    if(any(tmp))
    ans$regions[tmp] = lapply(ans$regions[tmp], FUN = function(r) {
      found = names(na.omit(which(s == sync_mut)))
      if(length(found) == 0) return(r)
      sync_typ = sync_type(r)
      sync_num = sync_part(r, sync_typ)
      if(length(found) != 1) {
        x = sapply(obj$regions[found], FUN = function(r) r$x[1])
        y = sapply(obj$regions[found], FUN = function(r) r$y[1])
        xx = sapply(regions[found], FUN = function(r) r$x[1])
        yy = sapply(regions[found], FUN = function(r) r$y[1])
        if((identical(sync_typ, "dual") && (length(unique(xx[xx != x])) > 1)) ||
           (identical(sync_typ, "quad") && (length(unique(xx[xx != x | yy != y])) > 1 || length(unique(yy[xx != x| yy != y])) > 1)) ) {
          stop("data_modify_regions: trying to set disparate coordinates for 'sync' regions:\n\t-",
               paste0(paste(found, sapply(regions[found], FUN = function(r) r$label), sep = " -> "), collapse="\n\t-"))
        }
      }
      if(sync_num != "") {
        coords = list(x = rep(regions[[found[1]]]$x[1], 2), y = rep(regions[[found[1]]]$y[1], 2))
        coords = switch(sync_typ,
                        "dual" = {
                          coords$y = rep(r$y[1], 2)
                          dual_coords(coords, sync_num)
                        },
                        "quad" = quad_coords(coords, sync_num),
                        coords)
        r$x <- coords$x
        r$y <- coords$y
      }
      return(r)
    })
  }
  tryCatch({
    ans$regions = sync_validate(ans$regions, quietly = FALSE, caller = "data_modify_regions")
  }, warning = function(w) stop(w$message))
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
                  args = c(list(obj = quote(ans),
                                old_names = old_pop_names,
                                new_names = new_pop_names),
                           dots[setdiff(names(dots), c("obj","old_names","new_pop_names"))]))
  }
  
  # since regions have been changed we need to recompute objects
  ans$pops <- do.call(what = popsCompute,
                      args = c(list(pops = quote(ans$pops),
                                    regions = quote(ans$regions),
                                    features = quote(ans$features),
                                    display_progress = display_progress),
                               dots[setdiff(names(dots),c("pops","features"))]))
  
  # modify stats
  if(nrow(obj$stats)!=0) ans$stats = get_pops_stats(ans$pops, as.integer(obj$description$ID$objcount))
  return(ans)
}
