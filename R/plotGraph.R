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

#' @title Plot and Stats Computation for IFC Graph
#' @description
#' Computes plot and stats from a IFC graph
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param graph a graph from 'obj' or a list that can be coerced by \code{\link{buildGraph}}.
#' @param draw whether to draw plot or not. Default is FALSE.
#' @param stats_print whether to print stats or not. Default is given by 'draw' argument.
#' @param color_mode whether to extract colors from 'obj' in white or black mode. Default is "white".
#' @param add_key whether to draw a "global" key under title or in the first "panel" or "both". Default is "panel".\cr
#' Accepted values are either: FALSE, "panel", "global", "both" or c("panel", "global").\cr
#' Note that it only applies when display is seen as overlaying populations.
#' @param precision when graphs is a 2D scatter with population overlay, this argument controls amount of information displayed. Default is "light".\cr
#' -"light", the default, will only display points of same coordinates that are among the other layers.\cr
#' -"full" will display all the layers.
#' @param trunc_labels maximum number of characters to display for labels. Default is 38.
#' @param trans the name of a transformation function for density graphs. If missing the default, the BasePop[[1]]$densitytrans, if any, will be retrieved, otherwise "asinh" will be used.
#' @param bin number of bin used for histogram / density. Default is missing.
#' @param viewport either "ideas", "data" or "max" defining limits used for the graph. Default is "ideas".\cr
#' -"ideas" will use same limits as the one defined in ideas.\cr
#' -"data" will use data to define limits.\cr
#' -"max" will use data and regions drawn to define limits.
#' @param ... other arguments to be passed.
#' @return it invisibly returns a list whose members are:\cr
#' -plot, "trellis" object that can be displayed using plot, if 'draw' was TRUE,\cr
#' -stats, a table of statistics computed for the graph, if 'stats_print' was TRUE,\cr
#' -input, a list with input parameters.
#' @export
plotGraph = function(obj, graph, draw = FALSE, stats_print = draw,
                     color_mode = c("white","black")[1], add_key = "panel", precision = c("light","full")[1],
                     trunc_labels = 38, trans = "asinh", bin, viewport = "ideas", ...) {
  dots = list(...)
  ret <- list()
  tryCatch({
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))

  # check mandatory param
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  assert(color_mode, len=1, alw=c("white","black"))
  color_mode=c(2,1)[c("white","black")==color_mode]
  assert(precision, len=1, alw=c("light","full"))
  if(!all(add_key%in%c("panel","global","both",FALSE))) stop("Accepted values for add_key are either: FALSE, 'panel', 'global', 'both' or c('panel', 'global')")
  trunc_labels = na.omit(as.integer(trunc_labels)); trunc_labels=trunc_labels[trunc_labels > 0]
  assert(trunc_labels, len=1, typ="integer")
  draw = as.logical(draw); assert(draw, len=1, alw=c(TRUE,FALSE))
  stats_print = as.logical(stats_print); assert(stats_print, len=1, alw=c(TRUE,FALSE))
  assert(viewport, len = 1, alw = c("ideas","data","max"))
  
  # shortcuts
  normalize = FALSE
  P = obj$pops
  R = obj$regions
  g = do.call(what=buildGraph, args=graph)
  foo = c(g$f1, g$f2)
  tmp = foo %in% names(obj$features)
  if(!all(tmp)) stop(paste0("trying to plot feature(s) not found in obj$features: ",  paste0(foo[!tmp], collapse=", ")))
  if(missing(trans)) trans = g$BasePop[[1]]$densitytrans
  if(length(trans) == 0 || (trans == "")) trans = "asinh"
  is_fun = !inherits(x = try(parseTrans(trans), silent = TRUE), what = "try-error")
  if((length(g$BasePop[[1]][["densitylevel"]]) != 0) && (g$BasePop[[1]][["densitylevel"]] != "")) trans = "return"
  
    # defines binning (for histogram and density)
  if(missing(bin)) {
    if(g$type=="histogram") {
      nbin = g$bincount
    } else {
      nbin = g$BasePop[[1]]$densitybincount
    }
    nbin = na.omit(as.integer(nbin))
    if(length(nbin)==0) nbin=ifelse(g$type=="histogram", 520, 128)
    if(nbin==0) nbin=ifelse(g$type=="histogram", 520, 128)
  } else {
    nbin = na.omit(as.integer(bin)); assert(nbin, len=1, typ="integer")
  }
  
  # extracts graph information
  if(g$type == "histogram") {
    D = obj$features[,c("Object Number",g$f1)]
    colnames(D) = c("Object Number","x1")
  } else {
    D = obj$features[,c("Object Number",g$f1,g$f2)]
    colnames(D) = c("Object Number","x1","y1")
  }
  Xlim = c(g$xmin, g$xmax)
  Ylim = c(g$ymin, g$ymax)
  Xtrans = g$xtrans; if(length(Xtrans) == 0) Xtrans = g$xlogrange
  Ytrans = g$ytrans; if(length(Ytrans) == 0) Ytrans = g$ylogrange
  trans_x <- parseTrans(Xtrans)
  trans_y <- parseTrans(Ytrans)
  # Xlim = Xlim + c(-0.07,0.07)*diff(Xlim) # fix, this should not be here error
  
  base_n = unlist(lapply(g$BasePop, FUN=function(x) x$name))
  reg_n = unlist(lapply(g$GraphRegion, FUN=function(x) x$name))
  shown_n = unlist(lapply(g$ShownPop, FUN=function(x) x$name))
  graph_n = unlist(lapply(g$GraphRegion, FUN=function(x) x$def))
  
  operators = c("And","Or","Not","(",")")
  displayed_n = unique(splitn(definition = g$order, all_names = c(base_n, graph_n, shown_n, "Selected Bin"), operators = operators))
  displayed_n = setdiff(displayed_n, "Selected Bin")
  displayed_r = rev(displayed_n)
  tmp = displayed_n %in% names(P)
  if(!all(tmp)) stop(paste0("trying to display a population not found in obj$pops: ",  paste0(displayed_n[!tmp], collapse=", ")))
  L = length(displayed_n)

  base_o = sapply(base_n, FUN=function(x) which(displayed_n%in%x))
  base_n = base_n[order(base_o)]
  
  if(length(shown_n) == 0) {
    shown_o = NULL
    shown_n = NULL
  } else {
    shown_o = unlist(sapply(shown_n, FUN=function(x) which(displayed_n%in%x)))
    shown_n = names(shown_o)[order(shown_o)]
  }
  
  # subset data
  # base = sapply(base_n, FUN=function(x) P[[x]]$obj)
  data_sub = fastAny(lapply(base_n, FUN=function(x) P[[x]]$obj))
  displayed_o = c(base_o, shown_o)
  Dall = fastCbind(D, sapply(displayed_n, simplify = FALSE, FUN=function(x) P[[x]]$obj), TRUE)
  D = Dall[data_sub, , drop = FALSE]
  dens_feat = numeric()
  
  xy_subset = rep(TRUE, nrow(D))
  if(length(xy_subset) == 0) xy_subset = TRUE
  D[,"x2"] = applyTrans(D[,"x1"], trans_x)
  Xlim = applyTrans(Xlim, trans_x)
  # computes limits / stats
  if(g$type == "histogram") {
    if(viewport == "data") {
      Xlim = cpp_fast_range(D[, "x1", drop=TRUE])
      Xlim = applyTrans(Xlim, trans_x)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
    }
    if(viewport == "max") {
      regx = sapply(reg_n, FUN=function(r) {
        reg = R[[r]] 
        coords = reg[["x"]]
        return(c(reg$cx, coords))
      })
      Xlim = cpp_fast_range(c(D[,"x1", drop=TRUE], unlist(regx)))
      Xlim = applyTrans(Xlim, trans_x)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
    }
    if(viewport == "ideas") {
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
      no_nas = !is.na(D[,"x2"])
      D[no_nas & (D[,"x2"] < Xlim[1]), "x2"] <- Xlim[1] # D = D[(D[,"x2"] >= Xlim[1]) & (D[,"x2"] <= Xlim[2]), ]
      D[no_nas & (D[,"x2"] > Xlim[2]), "x2"] <- Xlim[2] #
    }
    if(!all(is.finite(Xlim))) {
      Xlim = c(g$xmin, g$xmax)
      Xlim = applyTrans(Xlim, trans_x)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      if(!all(is.finite(Xlim))) Xlim = c(-1, 1)
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
      no_nas = !is.na(D[,"x2"])
      D[no_nas & (D[,"x2"] < Xlim[1]), "x2"] <- Xlim[1]
      D[no_nas & (D[,"x2"] > Xlim[2]), "x2"] <- Xlim[2]
    }
    smooth = (g$histogramsmoothingfactor != 0)
    br = do.breaks(Xlim, nbin)
    type = "count"
    if(as.logical(g$freq)) type = "percent"
    if(nrow(D) > 0) {
      if(viewport == "ideas") {
        Ylim = c(g$ymin, g$ymax)
        if(Ylim[1] == Ylim[2]) Ylim = Ylim[1] + c(0,0.07)
      } else {
        Ylim = c(0,max(sapply(displayed_n, FUN=function(x) get_ylim(x=D[D[,x],"x2"], type=type, br=br))*1.07))
        if(Ylim[1] == Ylim[2]) Ylim = Ylim[1] + c(0,0.07)
      } 
    } else {
      Ylim = c(g$ymin, g$ymax)
    }
    coln_stats = c("count","perc","Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
    stats = structure(matrix(numeric(), ncol = length(coln_stats), nrow = 0), dimnames = list(character(), coln_stats))
    colnames(stats) = c(coln_stats[1:2], paste0("x-",coln_stats[3:8]))
  } else {
    D[,"y2"] = applyTrans(D[,"y1"], trans_y)
    Ylim = applyTrans(Ylim, trans_y)
    if(viewport == "data") {
      Xlim = cpp_fast_range(D[,"x1", drop=TRUE])
      Xlim = applyTrans(Xlim, trans_x)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      Ylim = cpp_fast_range(D[,"y1", drop=TRUE])
      Ylim = applyTrans(Ylim, trans_y)
      Ylim = Ylim + c(-0.07,0.07)*diff(Ylim)
    }
    if(viewport == "max") {
      regx = sapply(reg_n, FUN=function(r) {
        reg = R[[r]] 
        coords = reg[["x"]]
        return(c(reg$cx, coords))
      })
      Xlim = cpp_fast_range(c(D[,"x1", drop=TRUE], unlist(regx)))
      Xlim = applyTrans(Xlim, trans_x)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      regy = sapply(reg_n, FUN=function(r) {
        reg = R[[r]] 
        coords = reg[["y"]]
        return(c(reg$cy, coords))
      })
      Ylim = cpp_fast_range(c(D[,"y1", drop=TRUE], unlist(regy)))
      Ylim = applyTrans(Ylim, trans_y)
      Ylim = Ylim + c(-0.07,0.07)*diff(Ylim)
    }
    if(!all(is.finite(Xlim))) {
      Xlim = c(g$xmin, g$xmax)
      Xlim = applyTrans(Xlim, trans_x)
      Xlim = Xlim + c(-0.07,0.07)*diff(Xlim)
      if(!all(is.finite(Xlim))) Xlim = c(-1, 1)
      if(Xlim[1] == Xlim[2]) Xlim = Xlim[1] + c(-0.07,0.07)
    }
    if(!all(is.finite(Ylim))) {
      Ylim = c(g$ymin, g$ymax)
      Ylim = applyTrans(Ylim, trans_y)
      Ylim = Ylim + c(-0.07,0.07)*diff(Ylim)
      if(!all(is.finite(Ylim))) Ylim = c(-1, 1)
      if(Ylim[1] == Ylim[2]) Ylim = Ylim[1] + c(-0.07,0.07)
    }
    if(nrow(D) > 0) {
      xy_subset = rep(FALSE, nrow(D))
      if(g$maxpoints <= 1) {
        xy_subset[cpp_fast_sample(n = nrow(D), size = g$maxpoints * nrow(D), replace = FALSE)] <- TRUE
      } else {
        xy_subset[cpp_fast_sample(n = nrow(D), size = min(g$maxpoints,nrow(D)), replace = FALSE)] <- TRUE
      }
    }
    if(!is_fun && (trans!="return")) {
      if(!any(names(obj$features) %in% trans)) stop(paste0("trying to plot a feature not found in obj$features: ", trans))
      dens_feat = obj$features[D[xy_subset, 1, drop = TRUE], trans, drop = TRUE]
      dens_ran = cpp_fast_range(dens_feat)
      dens_feat = (dens_feat-dens_ran[1])/diff(dens_ran)
    }
    coln_stats = c("count","perc","Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
    stats = structure(matrix(numeric(), ncol = length(coln_stats), nrow = 0), dimnames = list(character(), coln_stats))
    colnames(stats) = c(coln_stats[1:2], paste0("x-",coln_stats[3:8]), paste0("y-",coln_stats[9:14]))
  }
  
  # define text/points size
  dv <- dev.list()
  tryCatch({
    lt <- custom.theme(bg=c("black","white")[color_mode], fg=c("white","black")[color_mode])
    lt$grid.pars <- get.gpar()
  }, finally = while(!identical(dv, dev.list())) {
    dev.off(which = rev(dev.list())[1])
  })
  lt$grid.pars$fontfamily <- "serif"
  lt$fontsize$text <- lt$grid.pars$fontsize
  lt$fontsize$points <- 4
  for(i in c("xlab","xlab","zlab","main")) {
    lt[[paste0("par.",i,".text")]]$fontfamily <- "serif"
    lt[[paste0("par.",i,".text")]]$cex = g$axislabelsfontsize/lt$grid.pars$fontsize
  }
  lt[["axis.text"]]$fontfamily <- "serif"
  lt[["axis.text"]]$cex = g$axistickmarklabelsfontsize/lt$grid.pars$fontsize
  lt[["add.text"]]$fontfamily <- "serif"
  lt[["add.text"]]$cex = g$regionlabelsfontsize/lt$grid.pars$fontsize
  foo = list(par.settings = lt)
  if(g$type == "histogram") {
    ret_order = c("Object Number","x1","x2",displayed_n)
  } else {
    ret_order = c("Object Number","x1","x2","y1","y2",displayed_n)
  }
  displayed = lapply(obj$pops[displayed_n], FUN = function(p) {
    return(p[!(names(p) %in% "obj")])
  })
  suball = rep(FALSE, nrow(D))
  suball[D[xy_subset, 1, drop = TRUE]] <- TRUE
  ret = list("plot" = foo,
             "stats" = stats,
             "input" = list("data" = structure(D[ ,ret_order], features=dens_feat), 
                            "trunc_labels" = trunc_labels,
                            "title" = g$title,
                            "xlab" = g$xlab, "ylab" = g$ylab,
                            "xlim" = Xlim, "ylim" = Ylim, 
                            "trans_x" = Xtrans, "trans_y" = Ytrans,
                            "trans" = trans,
                            "order" = displayed_o,
                            "base" = g$BasePop,
                            "displayed" = displayed,
                            "graphical" = graph_n,
                            "regions" = obj$regions[reg_n],
                            "viewport" = viewport,
                            "bin" = nbin,
                            "type" = ifelse(g$type=="histogram", type, g$type),
                            "histogramsmoothingfactor" = g$histogramsmoothingfactor,
                            "normalize" = normalize,
                            "precision" = precision,
                            "add_key" = add_key,
                            "subset" = xy_subset,
                            "suball" = suball,
                            "mode" = color_mode))
  class(ret) <- "IFC_plot"
  invisible(ret)
  },
  finally = {
    if(stats_print) {
      ret$stats = plot_stats(ret)
      print(ret$stats)
    }
    if(draw) {
      tryCatch({
        ret$plot = plot_lattice(ret)
        plot(ret$plot)
      })
    }
    return(invisible(ret))
  })
}