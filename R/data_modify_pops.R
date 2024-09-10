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
#' @param pops a list of population(s) to modify in 'obj'. Each element of this list will be coerced by \code{\link{buildPopulation}}.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @details pops names should be present in names(obj$pops), otherwise an error will be raised.\cr
#' Note that If you want to rename pops, you should do it by changing 'name' member,
#' e.g. pops[[1]]$name <- "bar" while names(pops[[1]]) is "foo" and "foo" is part of names(obj$pops).
#' However, new names should not be present in 'obj'.
#' @param ... Other arguments to be passed.
#' @return an IFC_data object with pops modified.
#' @keywords internal
data_modify_pops <- function(obj, pops, display_progress = TRUE, ...){
  dots = list(...)
  assert(obj, cla = "IFC_data")
  mutation = names(pops)
  if(!all(mutation %in% names(obj$pops))) stop("can't find pops to modify in 'obj'", call. = FALSE)
  P = lapply(pops, keep_attributes, what = buildPopulation)
  names(P) = sapply(P, FUN = function(x) x$name)
  tmp = duplicated(names(P))
  if(any(tmp)) stop(paste0("duplicated pops found: ", unique(names(pops)[tmp])), call. = FALSE)
  names(mutation) = names(P)
  type1 = sapply(obj$pops[mutation], FUN = function(p) p$type)
  type2 = sapply(P, FUN = function(p) p$type)
  tmp = type1 == type2
  if(!all(tmp)) stop(paste0("'type' modification is not allowed:\n\t- ",
                            paste0(paste(paste0(mutation[!tmp], " [", type1[!tmp], "]"),
                                         paste0(names(mutation)[!tmp], " [", type2[!tmp], "]"), sep = " -> "),
                                   collapse = "\n\t- ")), call. = FALSE)
  
  # not any mutation should be found in already existing names
  # it is safer to stop because what to do if user rename a population but does not change its definition ?
  tmp = sapply(seq_along(mutation), FUN = function(i) any(names(mutation[i]) == names(obj$pops[setdiff(names(obj$pops), mutation[i])])))
  if(any(tmp)) stop("trying to rename population",ifelse(sum(tmp) > 1, "s", ""),
                    " with already existing name",ifelse(sum(tmp) > 1, "s", ""),
                    " use popsRename() instead:\n\t-",
                    paste0(names(mutation)[tmp], collapse="\n\t-"))
  
  # rename pops
  ans = do.call(what = popsRename,
                args = c(list(obj = quote(obj),
                              old_names = mutation,
                              new_names = names(mutation)),
                         dots[setdiff(names(dots), c("obj","old_names", "new_names"))]))
  
  # apply other modifications
  K = class(ans$pops)
  if(any(sapply(P, FUN = function(p) p$type == "C"))) {
    all_names = names(ans$pops)
    alt_names = gen_altnames(all_names) 
  }
  ans$pops[names(mutation)] = lapply(names(mutation), FUN = function(i_p) {
    p = ans$pops[[i_p]]
    p$color <- P[[i_p]]$color
    p$lightModeColor <- P[[i_p]]$lightModeColor
    p$style <- P[[i_p]]$style
    p$base <- P[[i_p]]$base
    if(p$type == "T") p$obj = P[[i_p]]$obj
    if(p$type == "C") {
      operators = c("And", "Or", "Not", "(", ")")
      p$definition = P[[i_p]]$definition
      p$split = splitn(definition = p$definition, all_names = all_names, alt_names = alt_names, operators = operators)
      p$names = setdiff(p$split, operators)
    }
    return(p)
  })
  names(ans$pops) = sapply(ans$pops, FUN = function(p) p$name)
  class(ans$pops) <- K
  ans$pops <- do.call(what = popsCompute,
                      args = c(list(pops = quote(ans$pops),
                                    regions = quote(ans$regions),
                                    features = quote(ans$features),
                                    display_progress = display_progress),
                               dots[setdiff(names(dots),c("pops","features"))]))
  
  # modify graphs (needed for example when base from a graphical population is changed)
  for(i_graph in seq_along(ans$graphs)) ans$graphs[[i_graph]] = adjustGraph(ans, ans$graphs[[i_graph]])

  # modify attributes
  for(i in names(mutation)) {
    attributes(ans$pops[[i]])[setdiff(names(attributes(ans$pops[[i]])), "names")] <- NULL
    for(k in setdiff(names(attributes(pops[[mutation[i]]])),"names")) {
      attr(ans$pops[[i]], k) <- attr(pops[[mutation[i]]], k)
    }
  }
  
  # modify stats
  if(nrow(obj$stats)!=0) ans$stats = get_pops_stats(ans$pops, as.integer(obj$description$ID$objcount))
  return(ans)
}
