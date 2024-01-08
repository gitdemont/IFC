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

#' @title Change Graphs Layout
#' @description
#' Reconstructs `IFC_graphs` object layout.
#' @param graphs an `IFC_graphs` object extracted with features extracted.
#' @param size Integer, graphs' sizes. Default is 320.
#' @param splitterdistance Integer. Default is 120.
#' @param stats Logical. Whether to show stats or not. Default is TRUE.
#' @param byrow Logical. Whether to layout graphs by row or not. Default is FALSE.
#' @param times Integer. Max number of graphs by row/column (depending on `byrow`). Default is 4.
#' @param layout Integer matrix. Desired layout. Default is NULL.
#' When NULL, the default, graphs will be relayout using `times` and `byrow` parameters.
#' Otherwise, layout will be used to determine the position of the graphs (NA value can be used for empty place).
#' Note that when not NULL layout should contain all indices of graphs only once.
#' @return an  `IFC_graphs` object
#' @keywords internal
relayout <- function(graphs, size = 320, splitterdistance = 120, stats = TRUE, byrow = FALSE, times = 4, layout = NULL) {
  assert(graphs, cla="IFC_graphs")
  size = suppressWarnings(na.omit(as.integer(size[size>0]))); assert(size,len=1,typ="integer")
  splitterdistance = suppressWarnings(na.omit(as.integer(splitterdistance[splitterdistance>0]))); assert(splitterdistance,len=1,typ="integer")
  assert(stats, len=1, alw=c(TRUE,FALSE))
  xsize = size
  ysize = size
  if(stats) ysize = ysize + splitterdistance
  if(length(layout) != 0) {
    for(i_graph in seq_along(graphs)) {
      tmp = which(layout == i_graph, arr.ind = TRUE, useNames = FALSE)
      if(!(1 %in% nrow(tmp))) stop("'layout' should contain all indices of graphs at most once")
      graphs[[i_graph]]$xsize <- xsize
      graphs[[i_graph]]$ysize <- ysize
      graphs[[i_graph]]$splitterdistance <- splitterdistance
      graphs[[i_graph]]$stats <- ifelse(stats, "true", "false")
      graphs[[i_graph]]$xlocation <- tmp[1] * xsize
      graphs[[i_graph]]$ylocation <- tmp[2] * ysize
    }
  } else {
    assert(byrow, len=1, alw=c(TRUE,FALSE))
    times = suppressWarnings(na.omit(as.integer(times[times>0]))); assert(times,len=1,typ="integer")
    for(i_graph in seq_along(graphs)) {
      xloc = ((i_graph - 1) %% times)
      yloc = ((i_graph - 1) %/% times)
      graphs[[i_graph]]$xsize <- xsize
      graphs[[i_graph]]$ysize <- ysize
      graphs[[i_graph]]$splitterdistance <- splitterdistance
      graphs[[i_graph]]$stats <- ifelse(stats, "true", "false")
      if(byrow) {
        graphs[[i_graph]]$xlocation <- yloc * xsize
        graphs[[i_graph]]$ylocation <- xloc * ysize
      } else {
        graphs[[i_graph]]$xlocation <- xloc * xsize
        graphs[[i_graph]]$ylocation <- yloc * ysize
      }
    } 
  }
  return(graphs)
}

#' @title Report Layout Extraction
#' @description
#' Extracts report layout from `IFC_graphs` object.
#' @param graphs an `IFC_graphs` object extracted with features extracted.
#' @return a list containing:\cr
#' -lay, a 3 columns (N, x, y) data.frame, where N is the graph index and x and y its coordinates on the layout,\cr
#' -mat, a matrix describing the layout.
#' @keywords internal
layoutReport <- function(graphs) {
  assert(graphs, cla="IFC_graphs")
  lay=structure(as.data.frame(matrix(integer(),nrow=0, ncol=3)), names=c("N","x","y"))
  lay_mat=structure(matrix(integer(),nrow=0,ncol=0),"row.vars"=list(y=character()), "col.vars"=list(x=character()))
  if(length(graphs) != 0) {
    lay=lapply(graphs, FUN=function(g) with(g, c(x=xlocation,y=ylocation)))
    lay=as.data.frame(do.call(rbind, lay))
    row.names(lay)=seq_along(graphs)
    lay=by(lay, lay$y, FUN=function(d) {
      d=d[order(d$x), ]
      d$x=seq_along(d$x)
      cbind("N"=as.integer(row.names(d)), d, stringsAsFactors=FALSE)
    })
    lay=lapply(seq_along(c(lay)), FUN=function(i) {
      d=lay[[i]]
      d$y=i
      return(d)
    })
    lay=do.call("rbind", c(lay, make.row.names=FALSE))
    lay_mat=ftable(by(lay$N, lay[,c("y","x")], FUN=function(x) x))
  }
  return(list(lay=lay, mat=unclass(lay_mat)))
}

#' @title Report File Creation
#' @description
#' Checks if report files can be created.
#' @param fileName default fileName path.
#' @param write_to pattern used to export file(s).
#' Placeholders, like c("\%d/\%s_fromR.pdf", "\%d/\%s_fromR.csv"), will be substituted:\cr
#' -\%d: with full path directory of 'fileName'\cr
#' -\%p: with first parent directory of 'fileName'\cr
#' -\%e: with extension of 'fileName' (without leading .)\cr
#' -\%s: with shortname from 'fileName' (i.e. basename without extension).\cr
#' Exported file(s) extension(s) will be deduced from this pattern. Note that has to be a .pdf and/or .csv.
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' Note that if TRUE, it will overwrite file. In addition a warning message will be sent.
#' @return a list with path(s) to pdf and csv file and overwritten status
#' @keywords internal
tryReportFileCreation <- function(fileName, write_to, overwrite = FALSE) {
  if(missing(write_to)) stop("'write_to' can't be missing")
  assert(write_to, typ="character")
  file_extension = getFileExt(write_to)
  assert(file_extension, alw=c("pdf", "csv"))
  if(any(duplicated(file_extension))) stop("'write_to' has to be only one .pdf and / or only one .csv")
  
  create_pdf = FALSE
  create_csv = FALSE
  export_to_pdf = NULL
  export_to_csv = NULL
  
  splitf_obj = splitf(fileName)
  if(any(file_extension%in%"pdf")) {
    create_pdf = TRUE
    splitp_obj_pdf = splitp(write_to[file_extension=="pdf"])
    if(any(splitp_obj_pdf$channel > 0)) message("'write_to' (pdf part) has %c argument but channel information can't be retrieved with tryReportFileCreation()")
    if(any(splitp_obj_pdf$object > 0)) message("'write_to' (pdf part) has %o argument but channel information can't be retrieved with tryReportFileCreation()")
    export_to_pdf = formatn(splitp_obj_pdf, splitf_obj)
  }
  if(any(file_extension%in%"csv")) {
    create_csv = TRUE
    splitp_obj_csv = splitp(write_to[file_extension=="csv"])
    if(any(splitp_obj_csv$channel > 0)) message("'write_to', (csv part) has %c argument but channel information can't be retrieved with tryReportFileCreation()")
    if(any(splitp_obj_csv$object > 0)) message("'write_to' (csv part) has %o argument but channel information can't be retrieved with tryReportFileCreation()")
    export_to_csv = formatn(splitp_obj_csv, splitf_obj)
  }
  write_to = c(export_to_pdf, export_to_csv)
  assert(overwrite, len=1, alw=c(TRUE,FALSE))
  overwritten = FALSE
  tmp = file.exists(write_to)
  if(any(tmp)) {
    if(!overwrite) stop(paste0("file ",paste0(write_to[tmp]," already exists"), collapse="\n"))
    if(create_pdf) if(file.exists(export_to_pdf)) export_to_pdf = normalizePath(export_to_pdf, winslash = "/")
    if(create_csv) if(file.exists(export_to_csv)) export_to_csv = normalizePath(export_to_csv, winslash = "/")
    overwritten = TRUE
  }
  
  if(create_pdf) {
    if(!dir.exists(dirname(export_to_pdf))) {
      if(!dir.create(dirname(export_to_pdf), showWarnings = FALSE, recursive = TRUE)) {
        stop(paste0(export_to_pdf,"\ncan't be created: check dirname"), call. = FALSE)
      }
    }
    tryCatch({
      pdf(file=export_to_pdf)
    }, error = function(e) stop(paste0(export_to_pdf,"\ncan't be created: check name / currently opened ?"), call. = FALSE))
    dev.off(dev.cur())
    export_to_pdf = normalizePath(export_to_pdf, winslash = "/")
  }
  if(create_csv) {
    if(!dir.exists(dirname(export_to_csv))) {
      if(!dir.create(dirname(export_to_csv), showWarnings = FALSE, recursive = TRUE)) {
        stop(paste0(export_to_csv,"\ncan't be created: check dirname"), call. = FALSE)
      }
    }
    tryCatch({
      write.table(x=rbind(c("pop","count","perc","x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.",
                            "y-Min.","y-1st Qu.","y-Median","y-Mean","y-3rd Qu.","y-Max.")),
                  sep=",", row.names = FALSE, col.names = FALSE, file=export_to_csv)
    }, error = function(e) stop(paste0(export_to_csv,"\ncan't be created: check name / currently opened ?"), call. = FALSE))
    export_to_csv = normalizePath(export_to_csv, winslash = "/")
  }
  write_to = c(export_to_pdf,export_to_csv)
  message(paste0(ifelse(length(write_to)==2, "files", "file")," will be exported in :\n",paste0(normalizePath(dirname(write_to), winslash = "/"),collapse="\n")))
  return(list(pdf = export_to_pdf, csv = export_to_csv, overwritten = overwritten))
}

#' @title Graph Report Generation
#' @description
#' Generates graph report (plot + statistics) from `IFC_data` object.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param selection indices of desired graphs. It can be provided as an integer vector or as a matrix.\cr
#' In such case, the layout of the matrix will reflect the layout of the extracted graphs.\cr
#' NA value will result in an empty place.\cr
#' Otherwise, when 'selection' is provided as a vector not identical to seq_along(obj$graphs), 'onepage' parameter will be set to FALSE.\cr
#' Note that indices are read from left to right, from top to bottom. Default is missing for extracting all graphs.
#' @param onepage whether to generate a pdf with all graphs on one page or not. Default is TRUE.
#' @param color_mode Whether to extract colors from obj in white or black mode. Default is "white".
#' @param add_key whether to draw a "global" key under title or in the first "panel" or "both". Default is "panel".\cr
#' Accepted values are either: FALSE, "panel", "global", "both" or c("panel", "global").\cr
#' Note that it only applies when display is seen as overlaying populations.
#' @param precision when graphs is a 2D scatter with population overlay, this argument controls amount of information displayed. Default is "light".\cr
#' -"light", the default, will only display points of same coordinates that are among the other layers.\cr
#' -"full" will display all the layers.
#' @param trunc_labels maximum number of characters to display for labels. Default is 38.
#' @param trans name of the transformation function for density graphs. If missing the default, the BasePop[[1]]$densitytrans, if any, will be retrieved, otherwise "asinh" will be used.
#' @param bin default number of bin used for histogram. Default is missing.
#' @param viewport Either "ideas", "data" or "max" defining limits used for the graph. Default is "ideas".\cr
#' -"ideas" will use same limits as the one defined in ideas.\cr
#' -"data" will use data to define limits.\cr
#' -"max" will use data and regions drawn to define limits.
#' @param backend backend used for drawing. Allowed are "lattice", "base", "raster". Default is "lattice".\cr
#' -"lattice" is the original one used in \pkg{IFC} using \pkg{lattice},\cr
#' -"base" will produce the plot using \pkg{base},\cr
#' -"raster" uses "base" for plotting but 2D graphs points will be produced as \code{\link[graphics]{rasterImage}}.\cr
#' This has the main advantage of being super fast allowing for plotting a huge amount of points while generating smaller objects (in bytes).
#' However, plot quality is impacted with "raster" method and resizing can lead to unpleasant looking.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other parameters to be passed.
#' @return a list with onepage, layout, layout_matrix, graphs, grobs, and stats.
#' @keywords internal
CreateGraphReport <- function(obj, selection, onepage=TRUE,
                              color_mode=c("white","black")[1], add_key="panel", precision=c("light","full")[1],    # parameters to pass to plotGraph
                              trunc_labels=38, trans="asinh", bin, viewport="ideas", backend="lattice",             # parameters to pass to plotGraph
                              display_progress=TRUE, ...) {
  dots = list(...)
  
  # backup dev.list
  dv <- dev.list()
  on.exit(dev_close(dv), add = TRUE)
  
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
  # check mandatory param
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
  assert(backend, len=1, alw=c("lattice","base","raster"))
  assert(onepage, len=1, alw=c(TRUE,FALSE))
  if(missing(bin) || (length(bin) == 0)) {
    bin = NULL
  } else {
    bin=na.omit(as.integer(bin)); assert(bin, len=1, typ="integer")
  }
  do_draw = na.omit(dots$draw)
  if(length(do_draw) != 1 || do_draw != FALSE) {
    do_draw = TRUE
  } else {
    do_draw = FALSE
  }
  dots$draw = FALSE
  do_stats = na.omit(dots$stats_print)
  if(length(do_stats) != 1 || do_stats != FALSE) {
    do_stats = TRUE
  } else {
    do_stats = FALSE
  }
  dots$stats_print = FALSE
  plotGraph_args = c(dots, list(obj=quote(obj), color_mode=color_mode, add_key=add_key,
                     precision=precision, trunc_labels=trunc_labels, viewport=viewport))
  if(length(bin) != 0) plotGraph_args = c(plotGraph_args, list(bin=bin))
  if(!missing(trans)) plotGraph_args = c(plotGraph_args, list(trans=trans))
  title_progress = basename(obj$fileName)
  
  # shortcuts
  G = obj$graphs
  if(length(G)==0) stop("there is no graph defined in 'obj'")
  
  # defines layout
  lay=layoutReport(G)
  lay_mat=lay$mat
  lay=lay$lay
  if(missing(selection) || (length(suppressWarnings(na.omit(as.integer(c(selection))))) == 0)) {
    selection = seq_along(G)
  } else {
    if(!all(na.omit(selection)%in%(seq_along(G)))) stop("'selection' refers to graph absent from 'obj'")
  }
  if(is.matrix(selection)) {
    G=G[c(selection)]
    lay_mat = matrix(seq_along(selection), nrow = nrow(selection))
    lay_mat[is.na(selection)] <- NA
    lay = data.frame(N = seq_along(selection), x = (1:ncol(selection)), y = (1:nrow(selection)))
    lay = lay[!is.na(selection),,drop=FALSE]
  } else {
    if(!identical(selection, seq_along(G))) onepage = FALSE
    if(onepage) {
      G=G[order(lay$N)]
    } else {
      G=G[na.omit(selection)]
    }
  }
  gl = length(G)
  
  # backup last state of graphic device
  if(display_progress) {
    pb_gr = newPB(min = 0, max = gl, initial = 0, style = 3)
    on.exit(endPB(pb_gr), add = TRUE)
  }
  default_stats_1D = matrix(NaN, nrow = 1, ncol = 8)
  colnames(default_stats_1D) = c("count","perc",
                                 "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.")
  default_stats_2D = matrix(NaN, nrow = 1, ncol = 14)
  colnames(default_stats_2D) = c("count","perc",
                                 "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.",
                                 "y-Min.","y-1st Qu.","y-Median","y-Mean","y-3rd Qu.","y-Max.")
  default_offscreen = function(w, h) {
    pdf(width=w, height=h, onefile=TRUE, pagecentre=TRUE, useDingbats=FALSE)
    dev.control("enable")
  }
  suppressWarnings({
    stList = list()
    grList = lapply(seq_along(G), FUN=function(i) {
      if(display_progress) setPB(pb = pb_gr, value = i, title = title_progress, label = "computing graphs and stats")
      stats = default_stats_1D
      rownames(stats) = paste0("Error: ", G[[i]]$title)
      if(length(G[[i]]) != 0) {
        g = do.call(what = plotGraph, args = c(plotGraph_args, list(graph = G[[i]])))
        if(G[[i]]$type=="histogram") {
          stats = default_stats_1D
        } else {
          stats = default_stats_2D
        }
        rownames(stats) = paste0("Error: ", G[[i]]$title)
        if(backend == "lattice") {
          foo = grob(p = do.call(what = ifelse(inherits(g, "error"),"plot_lattice_error",ifelse(inherits(g, "empty"),"plot_lattice_empty","plot_lattice")), args = list(obj = g)), vp = viewport(x=0.5, y=unit(0.5,"npc")), cl = "lattice")
        } else {
          foo = grid.grabExpr(grid.echo({plot_backend(g, backend); recordPlot()}, device = default_offscreen),
                              device = default_offscreen, width = 3*2.54-1.5, height = 3*2.54-0.5)
        }
        foo$children$`graphics-background`$width = unit(1,"npc")
        foo$children$`graphics-background`$height = unit(1,"npc")
      } else {
        do_stats = FALSE
        do_draw = FALSE
        foo = grob()
      }
      stats_try = stats
      if(do_stats) try({stats_try = plot_stats(g)}, silent = TRUE)
      if(!inherits(stats_try, "try-error")) stats = stats_try
      stList <<- c(stList, list(stats))
      stats = stats[,!grepl("Qu", colnames(stats)),drop=FALSE]
      foo_lay = matrix(c(rep(1,9), c(NA,2,NA)), nrow=4, ncol=3, byrow = TRUE)
      if(onepage) {
        if(nrow(stats) > 7) {
          stats = stats[1:8, ]
          stats[8, ] = "..."
        }
      }
      if(do_stats) {
        tab = tableGrob(format(stats, scientific=FALSE, digits=5), theme = ttheme_default(base_colour = "black", base_size=4, base_family = "sans"))
        tab$vp <- viewport(x=0.5, y=unit(1,"npc") - 0.5*unit(sum(tab$heights, na.rm=TRUE),"npc"))
        foo = arrangeGrob(foo, tab, layout_matrix = foo_lay, respect = TRUE)
      }
      return(foo)
    })
  })
  return(list(onepage = onepage, graphs = G, layout = lay, matrix = lay_mat, grobs = grList, stats = stList))
}

#' @title Graphical and Statistic Report Generation
#' @description
#' Generates report from `IFC_data` object.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param selection indices of desired graphs. It can be provided as an integer vector or as a matrix.\cr
#' In such case, the layout of the matrix will reflect the layout of the extracted graphs.\cr
#' NA value will result in an empty place.\cr
#' Otherwise, when 'selection' is provided as a vector not identical to seq_along(obj$graphs), 'onepage' parameter will be set to FALSE.\cr
#' Note that indices are read from left to right, from top to bottom. Default is missing for extracting all graphs.
#' @param write_to pattern used to export file(s).
#' Placeholders, like c("\%d/\%s_fromR.pdf", "\%d/\%s_fromR.csv"), will be substituted:\cr
#' -\%d: with full path directory of 'obj$fileName'\cr
#' -\%p: with first parent directory of 'obj$fileName'\cr
#' -\%e: with extension of 'obj$fileName' (without leading .)\cr
#' -\%s: with shortname from 'obj$fileName' (i.e. basename without extension).\cr
#' Exported file(s) extension(s) will be deduced from this pattern. Note that has to be a .pdf and/or .csv.
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' Note that if TRUE, it will overwrite file. In addition a warning message will be sent.
#' @param onepage whether to generate a pdf with all graphs on one page or not. Default is TRUE.
#' @param color_mode Whether to extract colors from obj in white or black mode. Default is "white".
#' @param add_key whether to draw a "global" key under title or in the first "panel" or "both". Default is "panel".\cr
#' Accepted values are either: FALSE, "panel", "global", "both" or c("panel", "global").\cr
#' Note that it only applies when display is seen as overlaying populations.
#' @param precision when graphs is a 2D scatter with population overlay, this argument controls amount of information displayed. Default is "light".\cr
#' -"light", the default, will only display points of same coordinates that are amoung the other layers.\cr
#' -"full" will display all the layers.
#' @param trunc_labels maximum number of characters to display for labels. Default is 38.
#' @param trans name of the transformation function for density graphs. If missing the default, the BasePop[[1]]$densitytrans, if any, will be retrieved, otherwise "asinh" will be used.
#' @param bin default number of bin used for histogram. Default is missing.
#' @param viewport Either "ideas", "data" or "max" defining limits used for the graph. Default is "ideas".\cr
#' -"ideas" will use same limits as the one defined in ideas.\cr
#' -"data" will use data to define limits.\cr
#' -"max" will use data and regions drawn to define limits.
#' @param backend backend used for drawing. Allowed are "lattice", "base", "raster". Default is "lattice".\cr
#' -"lattice" is the original one used in \pkg{IFC} using \pkg{lattice},\cr
#' -"base" will produce the plot using \pkg{base},\cr
#' -"raster" uses "base" for plotting but 2D graphs points will be produced as \code{\link[graphics]{rasterImage}}.\cr
#' This has the main advantage of being super fast allowing for plotting a huge amount of points while generating smaller objects (in bytes).
#' However, plot quality is impacted with "raster" method and resizing can lead to unpleasant looking.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other parameters to be passed.
#' @details depending on 'write_to', function will create .pdf and/or .csv file(s) report with according to graphs found in 'obj'.\cr
#' - csv file if created will contain "Min.","1st Qu.","Median","Mean","3rd Qu.","Max." for each graph found for x and y (if not histogram) for drawn populations and regions.\cr
#' - pdf file if created will contain graphs and to a certain extent some stats "Min.", "Median", "Mean", "Max." (no more than 7 rows).\cr
#' Note that only graphs will be exported (no images, features values, population stats, ...) in the same layout they were created and without sizing.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   tmp <- tempdir(check = TRUE)
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf, extract_images = FALSE,
#'                         extract_offsets = FALSE, display_progress = FALSE)
#'   L = length(daf$graphs)
#'   if(L > 0) { 
#'     ## randomly export at most 5 graphs from daf
#'     sel = sample(1:L, min(5, L))
#'     ExportToReport(obj = daf, selection = sel,
#'                    write_to = paste0(tmp, "\\test.pdf"), overwrite = TRUE)
#'   }
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return It invisibly returns full path of exported .pdf and/or .csv file(s).
#' @export
ExportToReport = function(obj, selection, write_to, overwrite=FALSE, onepage=TRUE,
                          color_mode=c("white","black")[1], add_key="panel", precision=c("light","full")[1],    # parameters to pass to plotGraph
                          trunc_labels=38, trans="asinh", bin, viewport="ideas", backend="lattice",             # parameters to pass to plotGraph
                          display_progress=TRUE, ...) {
  dots = list(...)
  dev_old <- options("device" = "pdf")
  on.exit(options(dev_old), add = TRUE)
  
  # backup dev.list
  dv <- dev.list()
  on.exit(dev_close(dv), add = TRUE)
  
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
  # check mandatory param
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
  title_progress = basename(obj$fileName)
  if(missing(bin)) bin = NULL
  
  # check file(s) can be created
  can_write = tryReportFileCreation(obj$fileName, write_to, overwrite)
  export_to_csv = can_write$csv
  export_to_pdf = can_write$pdf
  overwritten = can_write$overwritten
  create_csv = length(export_to_csv) != 0
  create_pdf = length(export_to_pdf) != 0
  write_to = c(export_to_pdf, export_to_csv)
  
  tryCatch({
    foo = CreateGraphReport(obj, selection, onepage=onepage,
                            color_mode=color_mode, add_key=add_key, precision=precision,
                            trunc_labels=trunc_labels, trans=trans, bin=bin, viewport=viewport, backend=backend,
                            display_progress=display_progress, ...)
    gl = length(foo$graphs)
    if(create_csv) {
      lapply(seq_along(foo$graphs), FUN = function(i) {
        g = foo$graphs[[i]]
        if(length(g)!=0) {
          write.table(x=rbind(c(g$f1,"x")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          if(g$type!="histogram") write.table(x=rbind(c(g$f2,"y")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          write.table(x=foo$stats[[i]], file=export_to_csv, append=TRUE, sep=",", col.names = FALSE)
        }
      })
    }
    if(create_pdf) {
      if(display_progress) {
        pb_pdf = newPB(title = title_progress,
                       label = ifelse(foo$onepage, "writing to pdf (no update but file is being processed)","writing to pdf"),
                       min = 0, max = gl, initial = 0, style = 3)
        on.exit(endPB(pb_pdf), add = TRUE)
      }
      if(foo$onepage) {
        # TODO add a progress bar
        pdf(file=export_to_pdf,
            width=3*max(foo$layout$x)*2.54, height=3*max(foo$layout$y)*2.54, 
            family="sans", onefile=TRUE, pagecentre=TRUE, useDingbats=FALSE)
        # no reason to have newpage = TRUE unless pdf is not open
        grid.arrange(grobs = foo$grobs[foo$layout$N], top = title_progress, newpage = names(dev.cur()) != "pdf", layout_matrix = foo$matrix, as.table = FALSE)
      } else {
        pdf(file=export_to_pdf, paper = "a4",
            family="sans", onefile=TRUE, pagecentre=TRUE, useDingbats=FALSE)
        # no reason to have newpage = TRUE unless pdf is not open
        grid.arrange(foo$grobs[[1]], top = title_progress, newpage = names(dev.cur()) != "pdf")
        if(gl > 1) for(i in 2:gl) {
          grid.arrange(foo$grobs[[i]], top = title_progress, newpage = TRUE)
          if(display_progress) setPB(pb_pdf, value = i, title = title_progress, label = "writing to pdf")
        }
      }
    }
    message(paste0("\n######################\n",ifelse(length(write_to)==2, "files have", "file has"), " been successfully ", ifelse(overwritten, "overwritten", "exported"),"\n"))
    return(invisible(write_to))
  }, error = function(e) {
    message(paste0(ifelse(length(write_to)==2, "files have", "file has"), " been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n"))
    stop(e$message, call. = FALSE)
  })
}

#' @title Graphical and Statistic Report Display
#' @description
#' Displays report from `IFC_data` object.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other parameters to be passed.
#' @return It invisibly returns NULL
#' @keywords internal
DisplayReport = function(obj, display_progress = TRUE, ...) {
  dots = list(...)
  
  # check mandatory param
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
  
  dots$obj <- quote(obj)
  dots$display_progress <- display_progress
  title_progress = basename(as.character(obj$fileName))
  
  foo = do.call(what = CreateGraphReport, args = dots)
  gl = length(foo$graphs)
  if(display_progress) {
    pb_pdf = newPB(title = title_progress,
                   label = ifelse(foo$onepage, "plotting graph merge (no update but file is being processed)","plotting each graph"),
                   min = 0, max = gl, initial = 0, style = 3)
    on.exit(endPB(pb_pdf), add = TRUE)
  }
  if(foo$onepage) {
    # TODO add a progress bar
    grid.arrange(grobs = foo$grobs[foo$layout$N], top = title_progress, newpage = TRUE, layout_matrix = foo$matrix, as.table = FALSE)
  } else {
    for(i in seq_along(foo$grobs)) {
      grid.arrange(foo$grobs[[i]], top = title_progress, newpage = TRUE)
      if(display_progress) setPB(pb_pdf, value = i, title = title_progress, label = "plotting each graph")
    }
  }
}

#' @title Batch Generation of Graphical and Statistic Report
#' @description
#' Batch creates graphical an statistical report.
#' @param fileName,obj either one or the other. Path to file(s) to read from for 'fileName' or list of `IFC_data` objects for obj.
#' @param selection indices of desired graphs. It can be provided as an integer vector or as a matrix.\cr
#' In such case, the layout of the matrix will reflect the layout of the extracted graphs for each 'fileName' or ''obj'.\cr
#' NA value will result in an empty place. When missing, it will be determined by the whole layout of 1st 'fileName' or 'obj' with 'gating' applied when provided
#' @param write_to pattern used to export file(s).
#' Placeholders, like c("\%d/\%s_fromR.pdf", "\%d/\%s_fromR.csv"), will be substituted:\cr
#' -\%d: with full path directory\cr
#' -\%p: with first parent directory\cr
#' -\%e: with extension (without leading .)\cr
#' -\%s: with shortname (i.e. basename without extension).\cr
#' Exported file(s) extension(s) will be deduced from this pattern using either 1st 'fileName' or 'obj'. Note that has to be a .pdf and/or .csv.
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' Note that if TRUE, it will overwrite file. In addition a warning message will be sent.
#' @param gating an `IFC_gating` object as extracted by readGatingStrategy(). Default is missing.
#' If not missing, each `IFC_data` provided in 'obj' or read from 'fileName' will be passed to applyGatingStrategy() before creating the report.
#' @param main the main title of the document. Default is missing.
#' @param byrow whether to add selected graphs for each file by row or not. Default is FALSE.
#' @param times number of files to add before starting a new row or column (depending on 'byrow').
#' @param color_mode Whether to extract colors in white or black mode. Default is "white".
#' @param add_key whether to draw a "global" key under title or in the first "panel" or "both". Default is "panel".\cr
#' Accepted values are either: FALSE, "panel", "global", "both" or c("panel", "global").\cr
#' Note that it only applies when display is seen as overlaying populations.
#' @param precision when graphs is a 2D scatter with population overlay, this argument controls amount of information displayed. Default is "light".\cr
#' -"light", the default, will only display points of same coordinates that are among the other layers.\cr
#' -"full" will display all the layers.
#' @param trunc_labels maximum number of characters to display for labels. Default is 38.
#' @param trans name of the transformation function for density graphs. If missing the default, the BasePop[[1]]$densitytrans, if any, will be retrieved, otherwise "asinh" will be used.
#' @param bin default number of bin used for histogram. Default is missing.
#' @param viewport Either "ideas", "data" or "max" defining limits used for the graph. Default is "ideas".\cr
#' -"ideas" will use same limits as the one defined in ideas.\cr
#' -"data" will use data to define limits.\cr
#' -"max" will use data and regions drawn to define limits.
#' @param backend backend used for drawing. Allowed are "lattice", "base", "raster". Default is "lattice".\cr
#' -"lattice" is the original one used in \pkg{IFC} using \pkg{lattice},\cr
#' -"base" will produce the plot using \pkg{base},\cr
#' -"raster" uses "base" for plotting but 2D graphs points will be produced as \code{\link[graphics]{rasterImage}}.\cr
#' This has the main advantage of being super fast allowing for plotting a huge amount of points while generating smaller objects (in bytes).
#' However, plot quality is impacted with "raster" method and resizing can lead to unpleasant looking.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other parameters to be passed.
#' @return It invisibly returns full path of exported .pdf and/or .csv file(s).
#' @export
BatchReport <- function(fileName, obj, selection, write_to, overwrite=FALSE, 
                        gating, main, byrow = FALSE, times = 5, 
                        color_mode=c("white","black")[1], add_key="panel", precision=c("light","full")[1],    # parameters to pass to plotGraph
                        trunc_labels=38, trans="asinh", bin, viewport="ideas", backend="lattice",             # parameters to pass to plotGraph
                        display_progress=TRUE, ...) {
  dots = list(...)
  create_pdf = FALSE
  create_csv = FALSE
  
  # backup dev.list
  dv <- dev.list()
  on.exit(dev_close(dv), add = TRUE)
  
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
  # check fileName or obj
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
  # check file(s) can be created
  assert(fileName, typ = "character")
  can_write = tryReportFileCreation(fileName[1], write_to, overwrite)
  export_to_csv = can_write$csv
  export_to_pdf = can_write$pdf
  overwritten = can_write$overwritten
  create_csv = length(export_to_csv) != 0
  create_pdf = length(export_to_pdf) != 0
  write_to = c(export_to_pdf, export_to_csv)
  
  # checks & coercion
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
  byrow = as.logical(byrow); assert(byrow, len=1, alw=c(TRUE,FALSE))
  times = na.omit(as.integer(times)); assert(times, len = 1)
  apply_gating = FALSE
  if(!missing(gating)) apply_gating = TRUE
  if(missing(bin)) {
    bin = NULL
  } else {
    bin=na.omit(as.integer(bin)); assert(bin, len=1, typ="integer")
  }
  if(!"draw" %in% names(dots)) dots$draw = FALSE
  if(!"stats_print" %in% names(dots)) dots$stats_print = FALSE
  plotGraph_args=c(dots, list(color_mode=color_mode, add_key=add_key,
                              precision=precision, trunc_labels=trunc_labels,
                              viewport=viewport))
  if(length(bin) != 0) plotGraph_args = c(plotGraph_args, list(bin=bin))
  if(!missing(trans)) plotGraph_args = c(plotGraph_args, list(trans=trans))
  
  # layout
  n_tiles = length(fileName)
  s_tiles = seq(from = 1, to = n_tiles, by = 1)
  lay = split(seq(from = 1, to = n_tiles, by = 1), ceiling(seq(from = 1, to = n_tiles, by = 1) / times) )
  if(length(lay) > 1) lay[[length(lay)]] = c(lay[[length(lay)]], rep(NA, times - length(lay[[length(lay)]])))
  lay = do.call(what = ifelse(byrow, "cbind", "rbind"), args = lay)
  if(missing(selection)) {
    if(apply_gating) {
      selection=gating$graphs
    } else {
      if(is_obj) {
        selection=obj[[1]]$graphs
      } else {
        selection=readIFC(fileName=fileName[1], display_progress=FALSE,
                          extract_features=TRUE, extract_images=FALSE,
                          extract_offsets=FALSE, extract_stats=FALSE)$graphs
      }
    }
    selection=layoutReport(selection)$mat
  }
  if(!is.matrix(selection)) {
    selection = na.omit(selection)
    selection = matrix(selection, nrow=length(selection))
    if(!byrow) selection = t(selection)
  }
  if(length(selection) == 0) selection=matrix(1,ncol=1,nrow=1)
  
  if(display_progress) {
    pb = newPB(min = 0, max = length(fileName), initial = 0, style = 3)
    on.exit(endPB(pb), add = TRUE)
  }
  tryCatch({
    grobs = lapply(seq_along(fileName), FUN = function(i_file) {
      if(display_progress) setPB(pb = pb, value = i_file, title = "creating batch report", label = paste0("extracting ",basename(fileName[i_file])))
      tryCatch({
        if(is_obj) {
          i_obj = obj[[i_file]]
        } else {
          i_obj = readIFC(fileName=fileName[i_file], display_progress=FALSE,
                          extract_features=TRUE, extract_images=FALSE,
                          extract_offsets=FALSE, extract_stats=FALSE)
        }
        if(apply_gating) i_obj = applyGatingStrategy(obj=i_obj, gating=gating, display_progress=FALSE)
        plotGraph_args = c(plotGraph_args, list(obj=i_obj))
        foo = CreateGraphReport(i_obj, selection, onepage=TRUE,
                                color_mode=color_mode, add_key=add_key, precision=precision,
                                trunc_labels=trunc_labels, trans=trans, bin=bin, viewport=viewport,
                                display_progress=display_progress)
        if(create_csv) {
          write.table(x=basename(fileName[i_file]), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          lapply(seq_along(foo$graphs), FUN = function(i) {
            g = foo$graphs[[i]]
            if(length(g)!=0) {
              write.table(x=rbind(c(g$f1,"x")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
              if(g$type!="histogram") write.table(x=rbind(c(g$f2,"y")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
              write.table(x=foo$stats[[i]], file=export_to_csv, append=TRUE, sep=",", col.names = FALSE)
            }
          })
        }
        grobTree(
          arrangeGrob(grobs = foo$grobs[foo$layout$N], 
                       top=textGrob(paste0("\n\n\n",basename(fileName[i_file]),"\n"), 
                                    gp=gpar(fontsize=14, font=2, col="skyblue4", lineheight=0.5)),
                       newpage = FALSE,
                       layout_matrix = foo$matrix, as.table = FALSE),
          rectGrob(.5,.5,width=unit(.99,"npc"), height=unit(0.99,"npc"), gp=gpar(lwd=2, col="skyblue4", fill=NA)))
      }, error = function(e) {
        grobTree(arrangeGrob(grid.text(label = e$message, gp=gpar(col="red"), draw = FALSE),
                             top=textGrob(paste0("\n\n\n",basename(fileName[i_file]),"\n"), 
                                          gp=gpar(fontsize=14, font=2, col="skyblue4", lineheight=0.5))),
                 rectGrob(.5,.5,width=unit(.99,"npc"), height=unit(0.99,"npc"), gp=gpar(lwd=2, col="skyblue4", fill=NA)))
      })
    })
    if(create_pdf) {
      pdf(file=write_to, height=3*nrow(lay)*nrow(selection)*2.54, width=3*ncol(lay)*ncol(selection)*2.54,
          family = "sans",
          onefile=TRUE, pagecentre=TRUE, useDingbats=FALSE)
      if(!missing(main)) args=c(args, top=main)
      # no reason to have newpage = TRUE unless pdf is not open
      args = list(newpage=names(dev.cur()) != "pdf",
                  layout_matrix=lay, as.table=FALSE)
      tryCatch(do.call(what=grid.arrange, args=c(list(grobs=quote(grobs)), args)), 
               error = function(e) { stop(e$message, call.=FALSE) },
               finally = dev.off())
    }
    message(paste0("\n######################\n",ifelse(length(write_to)==2, "files have", "file has"), " been successfully ", ifelse(overwritten, "overwritten", "exported"),"\n"))
    return(invisible(write_to))
  }, error = function(e) {
    message(paste0(ifelse(length(write_to)==2, "files have", "file has"), " been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n"))
    stop(e$message, call. = FALSE)
  })
}
