################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2021 Yohann Demont                                             #
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

################################################################################
#             functions described hereunder are experimental                   #
#              inputs and outputs may change in the future                     #
################################################################################

#' @title Apply Gating Strategy
#' @description
#' Applies Gating Strategy to an `IFC_data` object
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param gating an `IFC_gating` object extracted by \code{\link{readGatingStrategy}}.
#' @param keep names of population(s) present in 'obj' that should not be overwritten using 'gating'.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @param ... other arguments to be passed.
#' @details /!\ Please note that all former gating strategy (i.e. regions, pops, graphs and stats) will be removed from returned object, with the exception of population(s) described in 'keep'.\cr
#' An error will be thrown if a feature is required to create a population or a graph but can't be found in 'obj'.\cr
#' When tagged population(s) is(are) imported, objects from this(these) population(s) outside 'obj' will be discarded.\cr
#' If this results in NULL, then all objects will be tagged.
#' @return A named list of class `IFC_data` with new regions, pops and graphs
#' @keywords internal
applyGatingStrategy = function(obj, gating, keep, display_progress = TRUE, verbose = FALSE, ...) {
  dots = list(...)
  # due to large population and feature names error returned if any can be long
  # so we may need to modify options if error happens, as a consequence we store and return to default on exit
  old_opt = getOption("warning.length")
  on.exit(options("warning.length" = old_opt))
  
  # check mandatory param
  assert(obj, cla = "IFC_data")
  assert(gating, cla = "IFC_gating")
  if(missing(keep)) keep = NULL
  keep = as.character(keep); assert(keep, typ = "character")
  to_remove = grepl("^regions$|^pops$|^graphs$|^stats$", names(obj))
  if(!all(to_remove)) {
    ans = obj[!to_remove]
  } else {
    stop("'obj' is not of correct format", call. = FALSE)
  }
  class(ans) = class(obj)
  tryCatch({
    if(length(gating$spillover$spillover) != 0) {
      compensate <- function(raw, spillover) {
        ans = as.matrix(raw) %*% solve(t(spillover))
        colnames(ans) <- colnames(raw)
        return(ans)
      }
      comp_table = gating$spillover$spillover
      tmp1 = rownames(comp_table) %in% names(ans$features)
      if(all(tmp1)) {
        tmp2 = names(ans$features) %in% rownames(comp_table)
        ans$features[,  tmp2] <- compensate(ans$features[,  tmp2], spillover = comp_table) 
      } else {
        warning("spillover can not be applied on 'obj', some feature(s) are missing:\n\t-",
                paste0(rownames(comp_table)[!tmp1], collapse="\n\t-"), call. = FALSE, immediate. = TRUE)
      }
    }
    
    ans$regions = gating$regions
    ans$graphs = gating$graphs
    ans$pops = gating$pops
    if(length(ans$pops) == 0) {
      ans$pops = obj$pops["All"]
    }
    names(ans$pops) = sapply(ans$pops, FUN = function(p) p$name)
    
    if(length(keep) != 0) {
      tmp = keep %in% names(obj$pops)
      if(!all(tmp)) {
        warning("population(s) listed in 'keep' can't be found in 'obj':\n\t-",
                paste0(keep[!tmp], collapse="\n\t-"), call. = FALSE, immediate. = TRUE)
        keep = keep[tmp]
      }
    }
    ans$pops = ans$pops[!(names(ans$pops) %in% setdiff(keep, c("All", "")))]
    
    tmp1 = NULL
    tmp2 = NULL
    if(length(ans$pops) != 0) {
      feat_for_pops = lapply(ans$pops, FUN = function(g) g[c("fx", "fy")])
      tmp1 = which(!sapply(feat_for_pops, FUN = function(g) length(unlist(g))==0 || all(unlist(g) %in% names(ans$features))))
    }
    if(length(ans$graphs) != 0) {
      feat_for_graphs = lapply(ans$graphs, FUN = function(g) g[c("f1", "f2")])
      tmp2 = which(!sapply(feat_for_graphs, FUN = function(g) all(unlist(g) %in% names(ans$features))))
    }
    
    msg = NULL
    if(length(tmp1)!=0 || length(tmp2)!=0) {
      if(length(tmp1)!=0) msg = c(msg, paste0(paste0("missing feature",ifelse(length(unlist(feat_for_pops[tmp1]))>1,"s", "")," to create pop",ifelse(length(tmp1)>1,"s", ""),":\n"),
                                              paste0(paste0("\t- ",
                                                            sapply(tmp1, FUN = function(x) paste0("[", names(ans$pops)[x],"] ", paste0(unlist(feat_for_pops[x]), collapse = " - "))),
                                                            collapse = "\n"))))
      if(length(tmp2)!=0) msg = c(msg, paste0(paste0("missing feature",ifelse(length(unlist(feat_for_graphs[tmp2]))>1,"s", "")," to create graph",ifelse(length(tmp2)>1,"s", ""),":\n"), 
                                              paste0(paste0("\t- ", 
                                                            sapply(tmp2, FUN = function(x) paste0("[", x,"] ", paste0(unlist(feat_for_graphs[x]), collapse = " - "))),
                                                            collapse = "\n"))))
      options("warning.length" = (length(tmp1) + length(tmp2)) * 14 + nchar(msg))
      stop(paste0(msg, collapse = "\n"), call. = FALSE)
    } 
    
    ##### forces tagged population import
    ans$pops = lapply(ans$pops, FUN = function(p) {
      if(p$type != "T") return(p)
      obj = obj$features$`Object Number` %in% p$obj
      S = sum(obj)
      if(S == 0) {
        warning("importing tagged population [",p$name,"] resulted in NULL and has been replaced with \"All\" objects import", call. = FALSE, immediate. = TRUE)
        p$obj = rep(TRUE, length(obj))
      } else {
        if(S != length(p$obj)) warning("mismatch found between imported tagged population [",p$name,"] and objects stored within 'obj'", call. = FALSE, immediate. = TRUE)
        p$obj = obj
      }
      return(p)
    })
    
    ##### uses pops from obj that are listed in keep
    if(length(keep) != 0) {
      ans$pops <- c(ans$pops, obj$pops[keep])
    }
    class(ans$pops) <- "IFC_pops"
    
    ##### determines which object belongs to each population and changes styles and colors
    ans$pops = popsCompute(pops = ans$pops, regions = ans$regions, features = ans$features,
                           display_progress = display_progress, ...)
    
    ##### retrieves name(s) of graphical population created by region applied in graph
    if(length(ans$graphs) > 0) {
      ans$graphs = lapply(ans$graphs, FUN = function(g) {
        if(length(g$GraphRegion) != 0) {
          g$GraphRegion = lapply(g$GraphRegion, FUN = function(r) {
            foo = sapply(ans$pops,
                         FUN = function(p) {
                           bar = (p$type == "G") && 
                             (p$region == r$name) && 
                             (p$base %in% unique(unlist(lapply(g$BasePop, FUN = function(b) b$name)))) &&
                             (g$f1 == p$fx)
                           if(ans$regions[[r$name]]$type != "line") bar = bar && (g$f2 == p$fy)
                           return(bar)
                         })
            return(c(r, list(def = names(which(foo)))))
          })
        }
        return(do.call(buildGraph, args = g))
      })
    }
    class(ans$graphs) <- "IFC_graphs"
    
    ##### checks that graphs are correctly defined
    tryCatch({
      lapply(ans$graphs, FUN = function(g) { p = plot_lattice(plotGraph(obj = ans, graph = g, draw = FALSE, stats_print = FALSE))})
    }, error = function(e) stop("something went wrong while creating graphs"))
    
    ##### computes stats
    ans$stats = data.frame(stringsAsFactors = FALSE, check.rows = FALSE, check.names = FALSE, t(sapply(names(ans$pops), FUN=function(p) {
      count = sum(ans$pops[[p]]$obj)
      base = ans$pops[[p]]$base
      type = ans$pops[[p]]$type
      if(base=="") base = "All"
      parent = sum(ans$pops[[base]]$obj)
      c("type" = type, "parent" = base, "count" = count, "perc_parent" = count/parent*100, "perc_tot" = count/ncol(ans$features)*100)
    })))
    ans$stats[,3] = as.numeric(ans$stats[,3])
    ans$stats[,4] = as.numeric(ans$stats[,4])
    ans$stats[,5] = as.numeric(ans$stats[,5])
  },
  error = function(e) {
    stop("can't apply gating strategy:\n", e$message, call. = FALSE)
  })
  return(ans)
}