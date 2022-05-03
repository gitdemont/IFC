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
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other parameters to be passed.
#' @return a list with onepage, layout, layout_matrix, graphs, grobs, and stats.
#' @keywords internal
CreateGraphReport <- function(obj, selection, onepage=TRUE,
                            color_mode=c("white","black")[1], add_key="panel", precision=c("light","full")[1],    # parameters to pass to plotGraph
                            trunc_labels=38, trans="asinh", bin, viewport="ideas",                                  # parameters to pass to plotGraph
                            display_progress=TRUE, ...) {
  dots = list(...)
  dv = dev.list()
  on.exit(while(!identical(dv, dev.list())) {
    dev.off(which = rev(dev.list())[1])
  })
  
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
  # check mandatory param
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
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
  plotGraph_args = c(dots, list(obj=obj, color_mode=color_mode, add_key=add_key,
                     precision=precision, trunc_labels=trunc_labels, viewport=viewport))
  if(length(bin) != 0) plotGraph_args = c(plotGraph_args, list(bin=bin))
  if(!missing(trans)) plotGraph_args = c(plotGraph_args, list(trans=trans))
  title_progress = basename(obj$fileName)
  
  # shortcuts
  G = obj$graphs
  if(length(G)==0) stop("there is no graph defined in 'obj'")
  
  # defines layout
  lay=lapply(G, FUN=function(g) with(g, c(x=xlocation,y=ylocation)))
  lay=as.data.frame(do.call(rbind, lay))
  row.names(lay)=seq_along(G)
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
    pb_gr = newPB(session = dots$session, min = 0, max = gl, initial = 0, style = 3)
    on.exit(endPB(pb_gr), add = TRUE)
  }
  default_stats_1D = matrix(NA, nrow = 1, ncol = 8)
  colnames(default_stats_1D) = c("count","perc",
                                 "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.")
  default_stats_2D = matrix(NA, nrow = 1, ncol = 14)
  colnames(default_stats_2D) = c("count","perc",
                                 "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.",
                                 "y-Min.","y-1st Qu.","y-Median","y-Mean","y-3rd Qu.","y-Max.")
  suppressWarnings({
    stList = list()
    grList = lapply(seq_along(G), FUN=function(i) {
      if(display_progress) setPB(pb = pb_gr, value = i, title = title_progress, label = "computing graphs and stats")
      stats = default_stats_1D
      rownames(stats) = paste0("Error: ", G[[i]]$title)
      if(length(G[[i]]) != 0) {
        g = try({
          p = do.call(what = plotGraph, args = c(plotGraph_args, list(graph = G[[i]])))
          p$plot = plot_lattice(p)
          p
        }, silent = TRUE)
        if(!do_draw || inherits(x = g, what = "try-error")) {
          foo = arrangeGrob(grid.text(label = paste0(ifelse(do_draw,"Error: ","'draw' is set to FALSE "), attr(x = g, which = "condition")$message), gp=gpar(col="red"), draw = FALSE),
                            top = textGrob(paste0("\n",G[[i]]$title), gp = gpar(fontsize = 8, font=2, lineheight=0.5)))
        } else {
          if(G[[i]]$type=="histogram") {
            stats = default_stats_1D
          } else {
            stats = default_stats_2D
          }
          rownames(stats) = paste0("Error: ", G[[i]]$title)
          foo = grob(p=g$plot, vp = viewport(x=0.5, y=unit(0.5,"npc")), cl = "lattice")
        }
      } else {
        do_stats = FALSE
        do_draw = FALSE
        foo = grob()
      }
      if(do_stats) try({stats = plot_stats(g)}, silent = TRUE)
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
        tab = tableGrob(format(stats, scientific=FALSE, digits=5), theme = ttheme_default(base_size=4, base_family = "serif"))
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
                          trunc_labels=38, trans="asinh", bin, viewport="ideas",                                  # parameters to pass to plotGraph
                          display_progress=TRUE, ...) {
  dots = list(...)
  dv = dev.list()
  
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
                            trunc_labels=trunc_labels, trans=trans, bin=bin, viewport=viewport,
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
        pb_pdf = newPB(session = dots$session, title = title_progress,
                       label = ifelse(foo$onepage, "writing to pdf (no update but file is being processed)","writing to pdf"),
                       min = 0, max = gl, initial = 0, style = 3)
        on.exit(endPB(pb_pdf), add = TRUE)
      }
      if(foo$onepage) {
        # TODO add a progress bar
        pdf(file=export_to_pdf,
            width=3*max(foo$layout$x)*2.54, height=3*max(foo$layout$y)*2.54, 
            family="serif", onefile=TRUE, pagecentre=TRUE, useDingbats=FALSE)
        # no reason to have newpage = TRUE unless pdf is not open
        grid.arrange(grobs = foo$grobs[foo$layout$N], top = title_progress, newpage = names(dev.cur()) != "pdf", layout_matrix = foo$matrix, as.table = FALSE)
      } else {
        pdf(file=export_to_pdf, paper = "a4",
            family="serif", onefile=TRUE, pagecentre=TRUE, useDingbats=FALSE)
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
  },
  finally = {
    while(!identical(dv, dev.list())) {
      dev.off(which = rev(dev.list())[1])
    }
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
  
  dots$obj <- obj
  dots$display_progress <- display_progress
  title_progress = basename(as.character(obj$fileName))
  
  foo = do.call(what = CreateGraphReport, args = dots)
  gl = length(foo$graphs)
  if(display_progress) {
    pb_pdf = newPB(session = dots$session, title = title_progress,
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
#' @param fileName path to file(s).
#' @param selection indices of desired graphs. It can be provided as an integer vector or as a matrix.\cr
#' In such case, the layout of the matrix will reflect the layout of the extracted graphs for each fileName.\cr
#' NA value will result in an empty place.\cr
#' Otherwise, when 'selection' is provided as a vector not identical to seq_along(obj$graphs), 'onepage' parameter will be set to FALSE.\cr
#' Note that indices are read from left to right, from top to bottom. Default is missing for extracting all graphs.
#' @param write_to pattern used to export file(s).
#' Placeholders, like c("\%d/\%s_fromR.pdf", "\%d/\%s_fromR.csv"), will be substituted:\cr
#' -\%d: with full path directory of first element of 'fileName'\cr
#' -\%p: with first parent directory of first element of 'fileName'\cr
#' -\%e: with extension of first element of 'fileName' (without leading .)\cr
#' -\%s: with shortname from first element of 'fileName' (i.e. basename without extension).\cr
#' Exported file(s) extension(s) will be deduced from this pattern. Note that has to be a .pdf and/or .csv.
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' Note that if TRUE, it will overwrite file. In addition a warning message will be sent.
#' @param gating an `IFC_gating` object as extracted by readGatingStrategy(). Default is missing.
#' If not missing, each file provided in 'fileName' will be extracted by readIFC() and 'gating' will be applied by applyGatingStrategy() to the returned object before creating the report
#' @param main the main title of the document. Default is missing.
#' @param byrow whether to add selected graphs for each file by row or not. Deafult is FALSE.
#' @param times number of files to add before starting a new row or column (depending on 'byrow').
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
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... other parameters to be passed.
#' @return It invisibly returns full path of exported .pdf and/or .csv file(s).
#' @keywords internal
BatchReport <- function(fileName, selection, write_to, overwrite=FALSE, 
                        gating, main, byrow = FALSE, times = 5, 
                        color_mode=c("white","black")[1], add_key="panel", precision=c("light","full")[1],    # parameters to pass to plotGraph
                        trunc_labels=38, trans="asinh", bin, viewport="ideas",                                  # parameters to pass to plotGraph
                        display_progress=TRUE, ...) {
  dots = list(...)
  create_pdf = FALSE
  create_csv = FALSE
  
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
  # check file(s) can be created
  can_write = tryReportFileCreation(fileName[1], write_to, overwrite)
  export_to_csv = can_write$csv
  export_to_pdf = can_write$pdf
  overwritten = can_write$overwritten
  create_csv = length(export_to_csv) != 0
  create_pdf = length(export_to_pdf) != 0
  write_to = c(export_to_pdf, export_to_csv)
  
  # checks
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
  byrow = as.logical(byrow); assert(byrow, len=1, alw=c(TRUE,FALSE))
  times = na.omit(as.integer(times)); assert(times, len = 1)
  if(missing(selection)) stop("'selection' can't be missing")
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

  n_tiles = length(fileName)
  s_tiles = seq(from = 1, to = n_tiles, by = 1)
  lay = split(seq(from = 1, to = n_tiles, by = 1), ceiling(seq(from = 1, to = n_tiles, by = 1) / times) )
  if(length(lay) > 1) lay[[length(lay)]] = c(lay[[length(lay)]], rep(NA, times - length(lay[[length(lay)]])))
  lay = do.call(what = ifelse(byrow, "cbind", "rbind"), args = lay)
  L = length(selection)
  if(!is.matrix(selection)) {
    selection = na.omit(selection)
    selection = matrix(selection, nrow=length(selection))
    if(!byrow) selection = t(selection)
  }
  apply_gating = FALSE
  if(!missing(gating)) apply_gating = TRUE

  # backup last state of graphic device
  dv = dev.list()
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = length(fileName), initial = 0, style = 3)
    on.exit(endPB(pb), add = TRUE)
  }
  tryCatch({
    grobs = lapply(seq_along(fileName), FUN = function(i_file) {
      if(display_progress) setPB(pb = pb, value = i_file, title = "creating batch report", label = paste0("extracting ",basename(fileName[i_file])))
      tryCatch({
        obj = readIFC(fileName=fileName[i_file], display_progress=FALSE, ...)
        if(apply_gating) obj = applyGatingStrategy(obj=obj, gating=gating, display_progress=FALSE)
        plotGraph_args = c(plotGraph_args, list(obj=obj))
        foo = CreateGraphReport(obj, selection, onepage=TRUE,
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
          grid.arrange(grobs = foo$grobs[foo$layout$N], 
                       top=textGrob(paste0("\n\n\n",basename(fileName[i_file]),"\n"), 
                                    gp=gpar(fontsize=14, font=2, col="skyblue4", lineheight=0.5)),
                       newpage = names(dev.cur()) != "pdf", 
                       layout_matrix = foo$matrix, as.table = FALSE),
          rectGrob(.5,.5,width=unit(.99,"npc"), height=unit(0.99,"npc"), gp=gpar(lwd=2, col="skyblue4", fill=NA)))
      }, error = function(e) {
        if(names(dev.cur()) != "pdf") grid.newpage()
        grobTree(arrangeGrob(grid.text(label = e$message, gp=gpar(col="red"), draw = FALSE),
                             top=textGrob(paste0("\n\n\n",basename(fileName[i_file]),"\n"), 
                                          gp=gpar(fontsize=14, font=2, col="skyblue4", lineheight=0.5))),
                 rectGrob(.5,.5,width=unit(.99,"npc"), height=unit(0.99,"npc"), gp=gpar(lwd=2, col="skyblue4", fill=NA)))
      })
    })
    if(create_pdf) {
      pdf(file=write_to, height=3*nrow(lay)*nrow(selection)*2.54, width=3*ncol(lay)*ncol(selection)*2.54,
          family = "serif",
          onefile=TRUE, pagecentre=TRUE, useDingbats=FALSE)
      if(!missing(main)) args=c(args, top=main)
      # no reason to have newpage = TRUE unless pdf is not open
      args = list(newpage=names(dev.cur()) != "pdf", layout_matrix=lay, as.table=FALSE)
      tryCatch(do.call(what=grid.arrange, args=c(list(grobs=grobs), args)), 
               error = function(e) { stop(e$message, call.=FALSE) },
               finally = dev.off())
    }
    message(paste0("\n######################\n",ifelse(length(write_to)==2, "files have", "file has"), " been successfully ", ifelse(overwritten, "overwritten", "exported"),"\n"))
    return(invisible(write_to))
  }, error = function(e) {
    message(paste0(ifelse(length(write_to)==2, "files have", "file has"), " been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n"))
    stop(e$message, call. = FALSE)
  },
  finally = {
    while(!identical(dv, dev.list())) {
      dev.off(which = rev(dev.list())[1])
    }
  })
}
