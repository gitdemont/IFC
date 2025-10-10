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

#' @title Populations Renaming
#' @description
#' Renames populations in an `IFC_data` object
#' @param obj an `IFC_data`.
#' @param old_names character vector of name(s) of population(s) to rename inside 'obj'. Default is character().
#' @param new_names character vector of desired new population(s) name(s). Default is character().
#' @param loops a positive integer specifying the maximum number of recursive loops before raising an error. Default is 10L.
#' @param verbose whether to show a final message about the renaming. Default is TRUE.
#' @param ... other arguments to be passed.
#' @return an object of class `IFC_data`.
#' @export
popsRename <- function(obj, old_names = character(), new_names = character(), loops = 10L, verbose = TRUE, ...) {
  assert(obj, cla = "IFC_data")
  assert(old_names, typ="character")
  assert(new_names, typ="character")
  if(length(old_names) != length(new_names)) stop("'old_names' and 'new_names' should be character vectors of same length")
  sam = seq_along(old_names, FUN = function(i) identical(old_names[i], new_names[i]))
  old_names = old_names[!sam]
  new_names = new_names[!sam]
  if(length(old_names) == 0) return(obj)
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  loops = na.omit(as.integer(loops)); loops = loops[loops > 0]; assert(loops, len = 1, typ="integer")
  if(any(old_names %in% c(character(), "", NA_character_, "All"))) stop("'old_names' should not be NA, \"\", nor \"All\"")
  if(any(new_names %in% c(character(), "", NA_character_, "All"))) stop("'new_names' should not be NA, \"\", nor \"All\"")
  ori_names <- names(obj$pops)
  tmp = old_names %in% ori_names
  if(!all(tmp)) warning("can't find population",ifelse(sum(!tmp) > 1, "s", "")," to rename:\n",
                         paste0("\t- ", old_names[!tmp], collapse="\n"))
  if(!any(tmp)) return(obj)
  alt_names = gen_altnames(new_names[tmp], forbidden = c(ori_names, new_names[tmp]))
  mutations = list(structure(old_names[tmp], names = alt_names),
                   structure(alt_names, names = new_names[tmp]))
  P = obj$pops[old_names[tmp]]
  
  # create function for sorting unique values (and keeping names)
  sort_unique <- function(x) {
    N = names(x)
    foo = duplicated(x)
    xx = structure(x[!foo], names = N[!foo])
    sort(xx)
  }
  
  Kp = class(obj$pops)
  Kg = class(obj$graphs)
  # start population(s) renaming
  for(i in 1:2) {
    mutation = sort_unique(mutations[[i]])
    mutation_back = structure(character(), names = character())
    count = 0L 
    K <- class(obj$pops)
    while(!identical(mutation, mutation_back) || (count > loops)) {
      count = count + 1L
      mutation_back = mutation
      N = names(mutation)
      for(i_pop in seq_along(obj$pops)) {
        foo = mutation %in% obj$pops[[i_pop]]$name
        if(any(foo)) obj$pops[[i_pop]]$name <- N[foo]
        if(obj$pops[[i_pop]]$type == "C"){
          foo <- which(mutation %in% obj$pops[[i_pop]]$names)
          for(i in foo) {
            obj$pops[[i_pop]]$names[obj$pops[[i_pop]]$names == mutation[i]] <- N[i]
            obj$pops[[i_pop]]$split[obj$pops[[i_pop]]$split == mutation[i]] <- N[i]
            obj$pops[[i_pop]]$definition = paste0(obj$pops[[i_pop]]$split, collapse = "|")
          }
        }
        foo <- mutation %in% obj$pops[[i_pop]]$base
        if(any(foo)) {
          obj$pops[[i_pop]]$base <- N[foo]
          if(obj$pops[[i_pop]]$type == "G") {
            bar = sub(paste0(obj$pops[[i_pop]]$region, " & "), "", obj$pops[[i_pop]]$name, fixed = TRUE)
            foo = mutation %in% bar
            if(any(foo)) {
              g_name <- paste0(obj$pops[[i_pop]]$region, " & ", N[foo])
              mutation <- sort_unique(c(mutation, structure(obj$pops[[i_pop]]$name, names = g_name)))
              N <- names(mutation)
              obj$pops[[i_pop]]$name <- g_name
            }
          }
        }
      }
    }
    all_names = names(obj$pops)
    alt_names = gen_altnames(names(obj$pops), forbidden = c(names(obj$pops), c(ori_names[tmp], new_names[tmp])))
    names(obj$pops) = sapply(obj$pops, FUN = function(p) p$name)
    N = names(mutation)
    
    # now modify graphs
    K = class(obj$graphs)
    obj$graphs <- sapply(obj$graphs, simplify = FALSE, USE.NAMES = TRUE, FUN = function(g) {
      # modify parameters that depend on populations
      g$BasePop = sapply(g$BasePop, simplify = FALSE, USE.NAMES = TRUE, FUN = function(gg) {
        foo <- mutation %in% gg$name
        if(any(foo)) gg$name <- N[foo]
        gg
      })
      g$GraphRegion = sapply(g$GraphRegion, simplify = FALSE, USE.NAMES = TRUE, FUN = function(gg) {
        foo <- mutation %in% gg$def
        if(any(foo)) gg$def <- N[foo]
        gg
      })
      g$ShownPop = sapply(g$ShownPop, simplify = FALSE, USE.NAMES = TRUE, FUN = function(gg) {
        foo <- mutation %in% gg$name
        if(any(foo)) gg$name <- N[foo]
        gg
      })
      # check if title is default or has been given by user, when default try to modify it
      bar = try(trimws(splitn(definition = g$title, all_names = all_names, alt_names = alt_names, split = " , ")), silent = TRUE)
      if(!inherits(bar, "try-error")) {
        g$title = paste0(sapply(bar, FUN = function(x) {
          foo = mutation %in% x
          if(any(foo)) return(N[foo])
          x
        }), collapse = " , ")
      } 
      # change population in order
      bar = try(splitn(definition = g$order, all_names = all_names, alt_names = alt_names), silent = TRUE)
      if(!inherits(bar, "try-error")) {
        g$order = paste0(sapply(bar, FUN = function(x) {
          foo = mutation %in% x
          if(any(foo)) return(N[foo])
          x
        }), collapse = "|")
      } else { # something went wrong: g$order is reset and will be recomputed by buildGraph
        g$order <- NULL
      }
      g$xstatsorder <- NULL
      do.call(buildGraph, args = g)
    })
  }
  # throw error if max allowed number of recursions has been reached
  if(count >= loops) stop("can't rename population(s), max number of recursive loops reached")
  
  # modify attributes
  for(i in seq_along(old_names[tmp])) {
    attributes(obj$pops[[new_names[tmp][i]]])[setdiff(names(attributes(obj$pops[[new_names[tmp][i]]])), "names")] <- NULL
    for(k in setdiff(names(attributes(P[[old_names[tmp][i]]])),"names")) {
      attr(obj$pops[[new_names[tmp][i]]], k) <- attr(P[[old_names[tmp][i]]], k)
    }
  }
  # throw error if resulting 'obj' has duplicated pops names
  names(obj$pops) = sapply(obj$pops, FUN = function(p) p$name)
  tmp = duplicated(names(obj$pops)) 
  if(any(tmp)) stop("population renaming results in duplicated names:\n\t-", paste0(unique(names(obj$pops)[tmp]),collapse="\n\t-"))
  class(obj$pops) <- Kp
  class(obj$graphs) <- Kg
  
  if(verbose) {
    if(length(mutation) > 0) {
      message("population",ifelse(length(mutation) > 1,"s have"," has")," been successfully renamed:\n\t- ", paste(old_names, new_names, sep = " -> ", collapse = "\n\t- "))
    } else {
      message("no population was renamed")
    }
  } 
  return(obj)
}
