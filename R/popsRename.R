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
  verbose = as.logical(verbose); assert(verbose, len = 1, alw = c(TRUE, FALSE))
  loops = na.omit(as.integer(loops)); loops = loops[loops > 0]; assert(loops, len = 1, typ="integer")
  if(any(old_names %in% c(character(), "", NA_character_, "All"))) stop("'old_names' should not be NA, \"\", nor \"All\"")
  if(any(new_names %in% c(character(), "", NA_character_, "All"))) stop("'new_names' should not be NA, \"\", nor \"All\"")
  if(length(old_names) != length(new_names)) stop("'old_names' and 'new_names' should be character vectors of same length")
  ori_names <- names(obj$pops)
  tmp1 = old_names %in% ori_names
  tmp2 = new_names %in% ori_names
  if(!all(tmp1)) warning("can't find population",ifelse(sum(!tmp1) > 1, "s", "")," to rename:\n",
                         paste0("\t- ", old_names[!tmp1], collapse="\n"))
  if(any(tmp2)) warning("trying to rename population",ifelse(sum(tmp2) > 1, "s", ""),"with already existing name",ifelse(sum(tmp2) > 1, "s", ""),":\n",
                        paste0("\t- ", new_names[tmp2], collapse="\n"))
  mutation <- old_names[tmp1 & !tmp2]
  names(mutation) <- new_names[tmp1 & !tmp2]
  
  # create function for sorting unique values (and keeping names)
  sort_unique <- function(x) {
    N = names(x)
    foo = duplicated(x)
    xx = structure(x[!foo], names = N[!foo])
    sort(xx)
  }
  
  # start population(s) renaming
  mutation = sort_unique(mutation)
  mutation_back = structure(character(), names = character())
  count = 0L 
  K <- class(obj$pops)
  while(!identical(mutation, mutation_back) || (count > loops)) {
    count = count + 1L
    mutation_back = mutation
    N = names(mutation)
    obj$pops <- lapply(obj$pops, FUN = function(p) {
      foo <- mutation %in% p$name
      if(any(foo)) p$name <- N[foo]
      if(p$type == "C"){
        foo <- mutation %in% p$names
        if(any(foo)) {
          p$names[p$names == mutation[foo]] <- N[foo]
          p$split[p$split == mutation[foo]] <- N[foo]
          p$definition = paste0(p$split, collapse = "|")
        }
      }
      foo <- mutation %in% p$base
      if(any(foo)) {
        p$base <- N[foo]
        if(p$type == "G") {
          bar = gsub(paste0(p$region, " & "), "", p$name, fixed = TRUE)
          foo = mutation %in% bar
          if(any(foo)) {
            g_name <- paste0(p$region, " & ", N[foo])
            mutation <<- sort_unique(c(mutation, structure(p$name, names = g_name)))
            N <<- names(mutation)
            p$name <- g_name
            
          }
        }
      }
      p
    })
    names(obj$pops) = sapply(obj$pops, FUN = function(p) p$name)
  }
  class(obj$pops) <- K
  if(count >= loops) stop("can't rename population(s), max number of recursive loops reached")
  N = names(mutation)
  
  # now modify graphs
  obj$graphs <- sapply(obj$graphs, simplify = FALSE, USE.NAMES = TRUE, FUN = function(g) {
    # modify parameters that depend on populations
    g$BasePop = sapply(g$BasePop, simplify = FALSE, USE.NAMES = TRUE, FUN = function(gg) {
      foo <- mutation %in% gg$name
      if(any(foo)) gg$name <- N[foo]
      gg
    })
    g$GraphRegion = sapply(g$GraphRegion, simplify = FALSE, USE.NAMES = TRUE, FUN = function(gg) {
      foo <- mutation %in% gg$name
      if(any(foo)) gg$name <- N[foo]
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
    bar = try(splitn(trimws(g$title), all_names = ori_names, split = " , "), silent = TRUE)
    if(!inherits(bar, "try-error")) {
      g$title = paste0(sapply(bar, FUN = function(x) {
        foo = mutation %in% x
        if(any(foo)) return(N[foo])
        x
      }), collapse = " , ")
    } 
    # change population in order
    bar = try(splitn(g$order, all_names = ori_names), silent = TRUE)
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
  if(verbose) {
    if(length(mutation) > 0) {
      message("population",ifelse(length(mutation) > 1,"s have"," has")," been successfully renamed:\n\t- ", paste(mutation, names(mutation), sep = " -> ", collapse = "\n\t- "))
    } else {
      message("no population was renamed")
    }
  } 
  return(obj)
}
################################# test
# library(IFC)
# f = system.file(package = "IFCdata", "extdata", "example.daf")
# daf = readIFC(f)
# ndaf = popsRename(daf, c("foo","Events","Cells"), c("bar", "Events2","Focused"))
