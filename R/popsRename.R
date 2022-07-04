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
#' @param old_names character vector of name(s) of population(s) to rename. Default is character().
#' @param new_names character vector of new population(s) name(s). Default is character().
#' @param loops an positive integer specifying the maximum number of recursive loops before raising an error. Default is 10L
#' @param ... other arguments to be passed.
#' @return an object of class `IFC_data`.
#' @keywords internal
popsRename <- function(obj, old_names = character(), new_names = character(), loops = 10L, ...) {
  assert(obj, cla = "IFC_data")
  assert(old_names, typ="character")
  assert(new_names, typ="character")
  loops = na.omit(as.integer(loops)); loops = loops[loops > 0]; assert(loops, len = 1, typ="integer")
  if(any(old_names %in% c(character(), "", "All"))) stop("'old_names' should not be \"\", nor \"All\"")
  if(any(new_names %in% c(character(), "", "All"))) stop("'new_names' should not be \"\", nor \"All\"")
  if(length(old_names) != length(new_names)) stop("'old_names' and 'new_names' should be character vectors of same length")
  ori_names <- names(obj$pops)
  foo = old_names %in% ori_names
  if(!all(foo)) warning("can't find population(s) to rename:\n", paste0("\t- ", old_names[!foo], collapse="\n"))
  mutation <- old_names[foo]
  names(mutation) <- new_names[foo]
  foo = new_names %in% ori_names
  if(any(foo)) stop("trying to rename population with an already existing name:\n", paste0("\t- ", new_names[foo], collapse="\n"))

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
    # check if title is default or has been given by user
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
    } else { # something went wrong order is reset and will be recomputed by buildGraph
      g$order <- NULL
    }
    g$xstatsorder <- NULL
    do.call(buildGraph, args = g)
  })
  return(obj)
}
