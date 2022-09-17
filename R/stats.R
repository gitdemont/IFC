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

#' @title IFC Statistics Coercion
#' @description
#' Helper to build a list to allow statistics export.
#' @param type stats's type. Allowed are "COUNT","PERCENT_GATED","CONCENTRATION","PERCENT","MEAN","MEDIAN","STDDEV","MAD","CV","MINIMUM","MAXIMUM","GEOMETRIC_MEAN","MODE","VARIANCE","NAN","MEAN_RD","MEDIAN_RD".
#' @param title stats's title. If missing, it will be determined thanks to def.
#' @param def definition of the stats.
#' @details when 'type' is:\cr
#' - "COUNT","PERCENT_GATED","CONCENTRATION", 'def' will be "",
#' - "PERCENT", 'def' has to be a population name,
#' - "MEAN_RD","MEDIAN_RD" 'def' has to be the concatenation of a feature name and a population name collapse with "|". e.g. "Area_MC|All",
#' - otherwise, 'def' has to be a feature name.
#' @param ... Other arguments to be passed.
#' @return a list containing all statistics information.
#' @keywords internal
buildStats <- function(type, title, def, ...) {
  alw_type = c("COUNT","PERCENT_GATED","CONCENTRATION","PERCENT","MEAN","MEDIAN","STDDEV","MAD",
               "CV","MINIMUM","MAXIMUM","GEOMETRIC_MEAN","MODE","VARIANCE","NAN","MEAN_RD","MEDIAN_RD")
  type <- match.arg(type, choices = alw_type, several.ok = FALSE)
  if(missing(title)) {
    title = ""
  } else {
    title=as.character(title); assert(title, len=1)
  }
  if(type %in% alw_type[1:3]) {
    def = "";
  } else {
    def=as.character(def); def=setdiff(def, ""); assert(def, len=1)   
  }
  return(list("type"=type, "title"=title, "def"=def))
}

#' @title Statistics Extraction from Populations
#' @description
#' Extracts populations statistics
#' @param pops list of populations
#' @param objcount total number of objects.
#' @return a data.frame with computed statistics.
#' @keywords internal
get_pops_stats <- function(pops, objcount = 0) {
  names(pops) = sapply(pops, FUN = function(p) p$name)
  stats = data.frame(stringsAsFactors = FALSE, check.rows = FALSE, check.names = FALSE, t(sapply(names(pops), FUN=function(p) {
    count = sum(pops[[p]]$obj, na.rm = TRUE)
    base = pops[[p]]$base
    type = pops[[p]]$type
    if(base=="") base = "All"
    parent = sum(pops[[base]]$obj, na.rm = TRUE)
    c("type" = type, "parent" = base, "count" = count, "perc_parent" = count/parent*100, "perc_tot" = count/objcount*100)
  })))
  stats[,3] = as.numeric(stats[,3])
  stats[,4] = as.numeric(stats[,4])
  stats[,5] = as.numeric(stats[,5])
  stats
}

#' @title Statistics Extraction
#' @description
#' Extracts statistics from `IFC_data` object
#' @param obj an `IFC_data` object
#' @param feat_name a feature name.
#' @param trans character string describing transformation used and its parameters. See parseTrans().
#' @return a data.frame with computed statistics.
#' @keywords internal
extractStats <- function(obj, feat_name, trans="P") {
  assert(obj, cla="IFC_data")
  stats = get_pops_stats(obj$pops, as.integer(obj$description$ID$objcount))
  if(!missing(feat_name)) {
    msg = try(assert(feat_name, len=1, alw=names(obj$features)), silent=TRUE)
    if(inherits(msg, what="try-error")) {
      warning(attr(msg, "condition")$message, call.=FALSE, immediate.=TRUE)
    } else {
      trans_s = parseTrans(trans)
      val = sapply(rownames(stats), FUN = function(p) {
        foo = obj$features[obj$pops[[p]]$obj, feat_name]
        foo = na.omit(foo)
        if(length(foo) == 0) return(structure(rep(NA, 6), names = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")))
        foo = applyTrans(foo, trans_s)
        return(summary(foo))
      })
      stats = cbind(stats, t(val))
    }
  }
  return(stats)
}

#' @title IFC_stats Computation
#' @description
#' Function used to compute `IFC_stats` object\cr
#' @param obj an `IFC_data` object
#' @param stats list of statistics. It will be coerced by buildStats()
#' @param width desired width. Default is 80 * (1 + length(stats)).
#' @param height desired height Default is 400.
#' @param xlocation x location in analysis workspace. Default is 0.
#' @param ylocation y location in analysis workspace. Default is 0.
#' @param ... Other arguments to be passed.
#' @return an object of class `IFC_pops`.
#' @keywords internal
statsCompute <- function(obj, stats, width=80*length(stats), height=240, xlocation=0, ylocation=0, ...) {
  dots=list(...)
  dots=dots[!grepl("stats", names(dots))]
  assert(obj, cla="IFC_data")
  width=as.integer(na.omit(width)); width=width[width>=0]; assert(width, len=1)
  height=as.integer(na.omit(height)); height=height[height>=0]; assert(height, len=1)
  xlocation=as.integer(na.omit(xlocation)); xlocation=xlocation[xlocation>=0]; assert(xlocation, len=1)
  ylocation=as.integer(na.omit(ylocation)); ylocation=ylocation[ylocation>=0]; assert(ylocation, len=1)
  alw_type = c("COUNT","PERCENT_GATED","CONCENTRATION", # no def
               "PERCENT", # needs pop
               "MEAN","MEDIAN","STDDEV","MAD","CV","MINIMUM","MAXIMUM","GEOMETRIC_MEAN","MODE","VARIANCE","NAN", # needs feature
               "MEAN_RD","MEDIAN_RD") # needs feature and pop
  stats = lapply(stats, FUN = function(i_stats) {
    i_stats = do.call(what = buildStats, args = i_stats)
    if(i_stats$type == "PERCENT") {
      if(!i_stats$def %in% names(obj$pops)) stop("can't find [",i_stats$def,"] in obj$pops")
      if(i_stats$title == "") i_stats$title = paste0("%",i_stats$def)
    } else {
      if(i_stats$type %in% alw_type[5:15]) {
        if(!i_stats$def %in% names(obj$features)) stop("can't find [",i_stats$def,"] in obj$features")
        if(i_stats$title == "") i_stats$title = paste(i_stats$def, toCapFirstOnly(i_stats$type), sep = ", ")
      } else {
        if(i_stats$type %in% alw_type[16:17]) {
          foo = try(splitn(definition = i_stats$def, all_names = c(names(obj$pops), names(obj$features)), operators = ""), silent = TRUE)
          if(inherits(foo, what = "try-error")) stop("can't find [",i_stats$def,"] in obj$features + obj$pops")
          if(length(foo) != 2) stop("bad stats definition length")
          if(!all(c(foo[1] %in% names(obj$features), foo[2] %in% names(obj$pops)))) stop("bad stats definition")
          if(i_stats$title == "") i_stats$title = paste(paste0(foo, collapse = ", "),
                                                         ifelse(i_stats$type == "MEAN_RD", "Mean RD", "Median RD"), sep = ", ")
        } else {
          if(i_stats$title == "") i_stats$title = switch(i_stats$type, "COUNT"="Count","PERCENT_GATED"="%Gated","CONCENTRATION"="Objects/mL")
        }
      }
    }
    i_stats
  })
  ans = list(useShortNames="true", width=width, heigth=height, xlocation=xlocation, ylocation=ylocation, stats=stats, dots)
  class(ans) = "IFC_stats"
  return(ans)
}

#' @title IFC_stats XML Conversion
#' @description 
#' Helper to convert stats (`IFC_stats` object) to XML nodes.
#' @param stats an `IFC_stats` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a xml_node.
#' @keywords internal
toXML2_stats = function(stats, verbose = FALSE) {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating stats node")
  assert(stats, cla = "IFC_stats")
  if(length(stats$stats)==0) return(xml_new_node(name = "StatisticsControl", text = ""))
  xml_new_node(name = "StatisticsControl", attrs = stats[!(grepl("stats", names(stats)) | sapply(stats, length) == 0)],
               .children = lapply(stats$stats, FUN=function(i_stats) {
                 xml_new_node(name = "StatisticColumn",
                              attrs = i_stats)
               }))
}
