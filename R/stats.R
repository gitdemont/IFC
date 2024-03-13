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

################################################################################
#             functions described hereunder are experimental                   #
#              inputs and outputs may change in the future                     #
################################################################################

################################################################################
#       functions to compute features / pops stats from IFC_data object        #
################################################################################

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
      col_stats = c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","mad","sd","var","NaN")
      trans_s = parseTrans(trans)
      val = sapply(rownames(stats), FUN = function(p) {
        fv = as.numeric(obj$features[obj$pops[[p]]$obj, feat_name])
        no_nas = na.omit(fv)
        foo = fv[is.finite(fv)]
        if(length(foo) == 0) return(structure(rep(NaN, length(col_stats)), names = col_stats))
        foo = applyTrans(foo, trans_s)
        vv = c(summary(foo)[1:6], 
               mad(foo),
               ifelse(length(no_nas) == 1, 0, sd(foo)),
               ifelse(length(no_nas) == 1, 0, var(foo)),
               sum(!is.finite(fv))) # ifelse(any(is.finite(fv)), sum(!is.finite(fv)), NaN))
        vv[!is.finite(vv)] <- NaN
        return(structure(vv, names = col_stats))
      })
      stats = cbind(stats, t(val))
    }
  }
  return(stats)
}

################################################################################
#       functions to build ad export statistics panel inside daf file          #
################################################################################

#' @title IFC Statistics Coercion
#' @description
#' Helper to build a list to allow statistics export.
#' @param obj an `IFC_data` object
#' @param stats list of statistics instructions, whose members are list containing 3 instructions:\cr
#' - 'type' stats's type. Allowed are "COUNT","PERCENT_GATED","CONCENTRATION","PERCENT","MEAN","MEDIAN","STDDEV","MAD","CV","MINIMUM","MAXIMUM","GEOMETRIC_MEAN","MODE","VARIANCE","NAN","MEAN_RD","MEDIAN_RD".\cr
#' - 'title' stats's title. If missing, it will be determined thanks to def.\cr
#' - 'def' definition of the stats.\cr
#' Default is:\cr
#' list(list(type="COUNT",title="Count",def=""),\cr
#'      list(type="PERCENT_GATED",title="\%Gated",def=""))
#' @details when stats$type is:\cr
#' - "COUNT","PERCENT_GATED","CONCENTRATION", stats$def will be "",\cr
#' - "PERCENT", stats$def has to be a population name,\cr
#' - "MEAN_RD","MEDIAN_RD" stats$def has to be the concatenation of a feature name and a population name collapse with "|". e.g. "Area_MC|All",\cr
#' - otherwise, stats$def has to be a feature name.
#' @param width desired width. Default is 80 * (1 + length(stats)).
#' @param height desired height Default is 400.
#' @param xlocation x location in analysis workspace. Default is 0.
#' @param ylocation y location in analysis workspace. Default is 0.
#' @param ... Other arguments to be passed.
#' @return an object of class `IFC_stats`.
#' @keywords internal
buildStats <- function(obj, stats = list(list(type="COUNT",title="Count",def=""), list(type="PERCENT_GATED",title="%Gated",def="")), 
                       width=80*(1 + length(stats)), height=240, xlocation=0, ylocation=0, ...) {
  dots=list(...)
  dots=dots[!grepl("stats", names(dots))]
  assert(obj, cla="IFC_data")
  width=as.integer(na.omit(width)); width=width[width>=0]; assert(width, len=1)
  height=as.integer(na.omit(height)); height=height[height>=0]; assert(height, len=1)
  xlocation=as.integer(na.omit(xlocation)); xlocation=xlocation[xlocation>=0]; assert(xlocation, len=1)
  ylocation=as.integer(na.omit(ylocation)); ylocation=ylocation[ylocation>=0]; assert(ylocation, len=1)
  alw_type = c("COUNT","PERCENT_GATED","CONCENTRATION",                                                          # no def
               "PERCENT",                                                                                        # needs pop
               "MEAN","MEDIAN","STDDEV","MAD","CV","MINIMUM","MAXIMUM","GEOMETRIC_MEAN","MODE","VARIANCE","NAN", # needs feature
               "MEAN_RD","MEDIAN_RD")                                                                            # needs feature and pop
  all_names = c(names(obj$pops), names(obj$features))
  alt_names = gen_altnames(all_names)
  check_stats <- function(type, title, def, ...) {
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
  stats = lapply(stats, FUN = function(i_stats) {
    i_stats = do.call(what = check_stats, args = i_stats)
    if(i_stats$type == "PERCENT") {
      if(!i_stats$def %in% names(obj$pops)) stop("can't find [",i_stats$def,"] in obj$pops")
      if(i_stats$title == "") i_stats$title = paste0("%",i_stats$def)
    } else {
      if(i_stats$type %in% alw_type[5:15]) {
        if(!i_stats$def %in% names(obj$features)) stop("can't find [",i_stats$def,"] in obj$features")
        if(i_stats$title == "") i_stats$title = paste(i_stats$def, toCapFirstOnly(i_stats$type), sep = ", ")
      } else {
        if(i_stats$type %in% alw_type[16:17]) {
          foo = try(splitn(definition = i_stats$def, all_names = all_names, alt_names = alt_names, operators = ""), silent = TRUE)
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

################################################################################
#            functions to extract and generate statistics report               #
################################################################################

#' @title Statistics Report Template Extraction
#' @description
#' Retrieves statistics report template from .ast / .daf files.
#' @param fileName path to file.
#' @details Allowed \code{'statistics'} names are: \code{"Count"},\code{"Mean"},\code{"\%Total"},\code{"\%Gated"},\code{"\%"},\code{"Objects/mL"},\code{"RD - Mean"},
#' \code{"Median"},\code{"CV"},\code{"stddev"},\code{"NaN"},\code{"MAD"},\code{"min"},\code{"RD - Median"},\code{"Variance"},\code{"max"},\code{"geomean"},\code{"Mode"}.\cr
#' For \code{"\%Total"},\code{"\%Gated"},\code{"\%"},\code{"RD - Mean"},\code{"RD - Median"}, \code{'type'} has to be \code{"ratio"} and both \code{'population1'}
#' and \code{'population2'} should be provided. Otherwise, \code{'type'} is \code{"value"} and only \code{'population1'} is mandatory.\cr
#' \strong{/!\\} Note that \code{"Mode"} and \code{"Objects/mL"} can't be determined and will result in \code{NA}.
#' @return a 6-columns character matrix describing report instructions:\cr
#' - \code{'name'}, for the desired name of exported \code{'statistics'},\cr
#' - \code{'type'}, for the type of stats to return (either \code{"value"} or \code{"ratio"}),\cr
#' - \code{'population1'}, determines the population on which \code{'statistics'} will be performed,\cr
#' - \code{'population2'}, determines the reference population (when \code{'type'} is \code{"ratio"}, see \strong{Details}),\cr
#' - \code{'feature'}, determines the feature's name on which \code{'statistics'} will be computed,\cr
#' - \code{'statistics'}, controls the mathematical function that will be applied (see \strong{Details}).
#' @keywords internal
getSTATSREPORT <- function(fileName) {
  if(length(fileName) != 1) stop("'fileName' should be of length 1")
  if(!file.exists(fileName)) stop(paste0("can't find ",fileName))
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf","ast"))
  fileName = enc2native(normalizePath(fileName, winslash = "/", mustWork = FALSE))
  end_node = ifelse(file_extension == "daf", '</Assay>', '</AssayTemplate>')
  toskip = cpp_scanFirst(fileName, charToRaw(end_node), start = 0, end = 0)
  if(toskip == 0) stop(fileName, "\ndoes not seem to be well formatted: ", end_node, " not found")
  toskip = toskip + nchar(end_node) - 1
  tmp_daf = read_xml(readBin(con = fileName, what = "raw", n = toskip), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  alw_names = c("name","type","population1","population2","feature","statistics")
  stats_node = xml_find_all(tmp_daf, "//StatisticsReports/StatisticsReportColumn")
  if(length(stats_node) == 0) return(data.frame(matrix(nrow = 0, ncol = length(alw_names), dimnames = c(list(NULL),list(alw_names)))))
  do.call(what = rbind,
          args = sapply(stats_node,
                        USE.NAMES = FALSE, simplify = FALSE,
                        FUN = function(x) {
                          foo = xml_attrs(x)
                          names(foo)[names(foo) == "population"] <- "population1"
                          structure(sapply(alw_names, USE.NAMES = FALSE, FUN = function(n) foo[n]), names = alw_names)
                        }))
}

#' @title Statistical Report Generation
#' @description
#' Generates stats report from `IFC_data` object.
#' @param obj an `IFC_data` object.
#' @param stats template defining stats to compute as extracted by getSTATSREPORT().
#' @return a named vector with extracted statistics.
#' @keywords internal
StatsReport <- function(obj, stats) {
  if(!any(nrow(stats))) return(numeric(0))
  structure(apply(stats, 1, FUN = function(x) {
    n1 = x["population1"]
    n2 = x["population2"]
    n3 = ""
    fn = x["feature"]
    p1 = NA_real_; p2 = NA_real_; p3 = NA_real_; fv1 = NA_real_; fv2 = NA_real_
    if(n1 %in% names(obj$pops)) {
      p1 = obj$pops[[n1]]$obj
      n3 = obj$pops[[x["population1"]]]$base
    }
    if(n2 %in% names(obj$pops)) {
      p2 = obj$pops[[n2]]$obj
    }
    if(n3 %in% names(obj$pops)) p3 = obj$pops[[n3]]$obj
    if(fn %in% names(obj$features)) {
      fv1 = na.omit(as.numeric(obj$features[p1, fn, drop = TRUE]))
      if(n2 %in% names(obj$pops)) fv2 = na.omit(as.numeric(obj$features[p2, fn, drop = TRUE]))
    }
    # /!\ in R sd(1) == NA_real_, whereas with IDEAS stddev of a feature with 1 unique value gives 0
    # /!\ in R sd(numeric()) == NA_real_, whereas with IDEAS stddev of a NULL feature gives NaN
    # same for var
    ans = suppressWarnings(
      switch(
        x["statistics"],
        "Count"      = sum(p1, na.rm = TRUE),
        "Mean"       = mean(fv1),
        "%Total"     = 100 * sum(p1, na.rm = TRUE) / length(p1),
        "%Gated"     = 100 * sum(p1, na.rm = TRUE) / sum(p3, na.rm = TRUE),
        "%"          = 100 * sum(p1, na.rm = TRUE) / sum(p2, na.rm = TRUE),
        "Objects/mL" = {
          # How to compute concentration ?
          return(NA_real_)
        },
        "RD - Mean"   = abs(mean(fv1) - mean(fv2)) / (ifelse(length(fv1) == 1, 0, sd(fv1)) + ifelse(length(fv2) == 1, 0, sd(fv2))),
        "Median"      = median(fv1),
        "CV"          = 100 * ifelse(length(fv1) == 1, 0, sd(fv1)) / mean(fv1),
        "stddev"      = ifelse(length(fv1) == 1, 0, sd(fv1)),
        "NaN"         = {
          if(n1 %in% names(obj$pops) && fn %in% names(obj$features)) {
            fv = as.numeric(obj$features[p1, fn, drop = TRUE])
            v = sum(!is.finite(fv))
            return(ifelse((length(p1) - v == 0), NaN, ifelse(any(is.finite(fv)), v, NaN)))
          } else {
            return(0)
          }
          ifelse(n1 %in% names(obj$pops) && fn %in% names(obj$features),
                               sum(is.na(as.numeric(obj$features[p1, fn, drop = TRUE]))),
                               0)
        },
        "MAD"         = mad(fv1),
        "min"         = min(fv1),
        "RD - Median" = abs(median(fv1) - median(fv2)) / (mad(fv1) + mad(fv2)),
        "Variance"    = ifelse(length(fv1) == 1, 0, var(fv1)),
        "max"         = max(fv1),
        "geomean"     = exp(mean(log(fv1), na.rm=TRUE)),
        "Mode" = {
          # How to compute concentration ?
          # What to do here ? the following does not work
          # 1/ as.numeric(names(which.max(table(fv1, useNA = "no"))))
          # 2/ mean(as.numeric(strsplit(gsub("\\(|\\]","",names(which.max(table(cut(fv1, breaks = 128))))),",",fixed = TRUE)[[1]]))
          return(NA_real_)
        },
        {
          stop("statistics [" , x["statistics"], "] is not allowed")
        }))
    if(!is.finite(ans)) return(NaN)
    ans
  }), names = stats[,"name"])
}

#' @title Batch Generation of Statistic Report
#' @description
#' Generates statistics report on batch of files or `IFC_data` objects.
#' @param fileName,obj either one or the other. Path to file(s) to read from for 'fileName' or list of `IFC_data` objects for obj.
#' @param stats template defining stats to compute as extracted by getSTATSREPORT().
#' @param gating an `IFC_gating` object as extracted by readGatingStrategy(). Default is missing.
#' If not missing, each `IFC_data` provided in 'obj' or read from 'fileName' will be passed to applyGatingStrategy() before creating the report.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other parameters to be passed.
#' @return a data.frame of statistics
#' @keywords internal
BatchStatsReport <- function(fileName, obj, stats, gating, display_progress = FALSE, ...) {
  # check fileName or obj
  dots = list(...)
  if(missing(obj) && missing(fileName)) stop("you should provide either 'fileName' or 'obj'")
  if(!missing(fileName) && !missing(obj)) stop("you should provide either 'fileName' or 'obj' not both")
  if(!missing(fileName) && is.list(fileName) && inherits(fileName[[1]], what="IFC_data")) {
    warning("'fileName' will be treated as a list of `IFC_data` objects")
    obj = fileName
  }
  is_obj = FALSE
  if(!missing(obj)) {
    is_obj = TRUE
    fileName = sapply(seq_along(obj), FUN = function(i_obj) {
      assert(obj[[i_obj]], cla="IFC_data")
      obj[[i_obj]]$fileName
    })
  }
  assert(fileName, typ = "character")
  display_progress = as.logical(display_progress); assert(display_progress, alw = c(TRUE,FALSE))
  apply_gating = FALSE
  if(!missing(gating)) apply_gating = TRUE
  if(display_progress) {
    pb = newPB(min = 0, max = length(fileName), initial = 0, style = 3)
    on.exit(endPB(pb))
  }
  ans = lapply(seq_along(fileName), FUN = function(i_file) {
    tryCatch({
      if(is_obj) {
        i_obj = obj[[i_file]]
      } else {
        i_obj = readIFC(fileName=fileName[i_file], display_progress=FALSE,
                        extract_features=TRUE, extract_images=FALSE,
                        extract_offsets=FALSE, extract_stats=FALSE)
      }
      if(apply_gating) i_obj = applyGatingStrategy(obj=i_obj, gating=gating, display_progress=FALSE)
      foo = StatsReport(i_obj, stats)
      if(length(foo) == 0) return(list())
      foo
    }, error = function(e) {
      warning("can't extract stats for fileName [", fileName[i_file],"]\n", e$message, call. = FALSE, immediate. = TRUE)
      return(list())
    }, finally = {
      if(display_progress) setPB(pb, value = i_file, title = "Computing Batch Stats", label = basename(fileName[i_file]))
    })
  })
  no_err = sapply(ans, length) != 0
  if(!any(no_err)) no_err = rep(TRUE, length(no_err))
  cbind("File" = basename(fileName[no_err]), data.frame(do.call(rbind, args = ans[no_err]), check.names = FALSE, check.rows = FALSE), deparse.level = 0)
}
