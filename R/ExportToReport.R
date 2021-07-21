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

#' @title Graphical and Statistic Report Generation
#' @description
#' Generates report from `IFC_data` object.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param selection when provided, indices of desired graphs.\cr
#' In such case onepage parameter is set to FALSE.\cr
#' Note that indices are read from left to right, from top to bottom. 
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
#' @param color_mode Whether to extract colors from obj in white or black mode. Default is 'white'.
#' @param add_key whether to draw a 'global' key under title or in the first 'panel' or 'both'. Default is 'panel'.\cr
#' Accepted values are either: FALSE, 'panel', 'global', 'both' or c('panel', 'global').\cr
#' Note that it only applies when display is seen as overlaying populations.
#' @param precision when graphs is a 2D scatter with population overlay, this argument controls amount of information displayed. Default is "light".\cr
#' -"light", the default, will only display points of same coordinates that are amoung the other layers.\cr
#' -"full" will display all the layers.
#' @param trunc_labels maximum number of characters to display for labels. Default is 38.
#' @param trans transformation function for density graphs. If missing the default, the BasePop[[1]]$densitytrans, if any, will be retrieved, otherwise asinh will be used.
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
                          trunc_labels=38, trans=asinh, bin, viewport="ideas",                                  # parameters to pass to plotGraph
                          display_progress=TRUE, ...) {
  dots = list(...)
  
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
  pdf_dev = TRUE
  create_pdf = FALSE
  create_csv = FALSE
  from = sys.nframe()
  if(from > 1) pdf_dev = (sub("^.*:(.*)$", "\\1", as.list(sys.call(-1))[[1]]) != "DisplayReport")
  
  # check mandatory param
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
  assert(onepage, len=1, alw=c(TRUE,FALSE))
  if(missing(bin)) {
    bin = NULL
  } else {
    bin=na.omit(as.integer(bin)); assert(bin, len=1, typ="integer")
  }
  plotGraph_args = list(obj=obj, draw=FALSE, color_mode=color_mode, add_key=add_key,
                        precision=precision, trunc_labels=trunc_labels, viewport=viewport)
  if(length(bin) != 0) plotGraph_args = c(plotGraph_args, list(bin=bin))
  if(!missing(trans)) plotGraph_args = c(plotGraph_args, list(trans=trans))
  
  title_progress = basename(obj$fileName)
  splitf_obj = splitf(obj$fileName)
  
  if(pdf_dev) {
    if(missing(write_to)) stop("'write_to' can't be missing")
    assert(write_to, typ = "character")
    file_extension = getFileExt(write_to)
    assert(file_extension, alw = c("pdf", "csv"))
    if(any(duplicated(file_extension))) stop("'write_to' has to be only one .pdf and / or only one .csv")
    export_to_pdf = NULL
    export_to_csv = NULL
    if(any(file_extension%in%"pdf")) {
      create_pdf = TRUE
      splitp_obj_pdf = splitp(write_to[file_extension=="pdf"])
      if(any(splitp_obj_pdf$channel > 0)) message("'write_to' (pdf part) has %c argument but channel information can't be retrieved with ExportToReport()")
      if(any(splitp_obj_pdf$object > 0)) message("'write_to' (pdf part) has %o argument but channel information can't be retrieved with ExportToReport()")
      export_to_pdf = formatn(splitp_obj_pdf, splitf_obj)
    }
    if(any(file_extension%in%"csv")) {
      create_csv = TRUE
      splitp_obj_csv = splitp(write_to[file_extension=="csv"])
      if(any(splitp_obj_csv$channel > 0)) message("'write_to', (csv part) has %c argument but channel information can't be retrieved with ExportToReport()")
      if(any(splitp_obj_csv$object > 0)) message("'write_to' (csv part) has %o argument but channel information can't be retrieved with ExportToReport()")
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
  }
  # shortcuts
  G = obj$graphs
  if(length(G)==0) stop("there is no graph defined in 'obj'")
  # defines layout
  lay=lapply(G, FUN=function(g) with(g, c(x=xlocation,y=ylocation)))
  lay=as.data.frame(do.call(rbind, lay))
  row.names(lay)=1:length(G)
  lay=by(lay, lay$y, FUN=function(d) {
    d=d[order(d$x), ]
    d$x=seq_along(d$x)
    cbind("N"=as.integer(row.names(d)), d, stringsAsFactors=FALSE)
  })
  lay=lapply(1:length(lay), FUN=function(i) {
    d=lay[[i]]
    d$y=i
    return(d)
  })
  lay=do.call("rbind", c(lay, make.row.names=FALSE))
  lay_mat=ftable(by(lay$N, lay[,c("y","x")], FUN=function(x) x))
  if(missing(selection)) {
    selection = 1:length(G)
  } else {
    if(!all(selection%in%(1:length(G)))) stop("'selection' refers to graph absent from 'obj'")
    onepage = FALSE
  }
  if(onepage) { 
    G=G[order(lay$N)]
  } else {
    G=G[lay$N[selection]]
  }
  gl = length(G)
  
  # backup last state of graphic device
  dv = dev.cur()
  if(display_progress) pb_gr = newPB(session = dots$session, min = 0, max = gl, initial = 0, style = 3)
  tryCatch({
    suppressWarnings({
      graphs = lapply(1:gl, FUN=function(i) {
        if(display_progress) setPB(pb = pb_gr, value = i, title = title_progress, label = paste0("computing ",ifelse(create_pdf,"graphs and ",""),"stats"))
        g = try(do.call(what = plotGraph, args = c(plotGraph_args, list(graph = G[[i]]))))
        if(inherits(x = g, what = "try-error")) {
          foo = arrangeGrob(grid.text(label = paste0("Error: ", attr(x = g, which = "condition")$message), gp=gpar(col="red"), draw = FALSE),
                            top = textGrob(paste0("\n",G[[i]]$title), gp = gpar(fontsize = 8, font=2, lineheight=0.5)))
          if(G[[i]]$type=="histogram") {
            stats = matrix(NA, nrow = 1, ncol = 8)
            colnames(stats) = c("count","perc",
                                "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.")
          } else {
            stats = matrix(NA, nrow = 1, ncol = 14)
            colnames(stats) = c("count","perc",
                                "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.",
                                "y-Min.","y-1st Qu.","y-Median","y-Mean","y-3rd Qu.","y-Max.")
          }
          rownames(stats) = paste0("Error: ", G[[i]]$title)
        } else {
          foo = g$plot
          stats = g$stats
        }
        if(create_csv) {
          write.table(x=rbind(c(G[[i]]$f1,"x")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          if(G[[i]]$type!="histogram") write.table(x=rbind(c(G[[i]]$f2,"y")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          write.table(x=stats, file=export_to_csv, append=TRUE, sep=",", col.names = FALSE)
        }
        stats = stats[,!grepl("Qu", colnames(stats)),drop=FALSE]
        foo_lay = matrix(c(rep(1,9), c(NA,2,NA)), nrow=4, ncol=3, byrow = TRUE)
        if(onepage) {
          if(nrow(stats)> 7) {
            stats = stats[1:8,]
            stats[8, ] = "..."
          }
        }
        foo$vp <- viewport(x=0.5, y=0.5)
        tab = tableGrob(format(stats, scientific=FALSE, digits=5), theme = ttheme_default(base_size=4, base_family = "serif"))
        tab$vp <- viewport(x=0.5, y=unit(1,"npc") - 0.5*sum(tab$heights))
        foo = arrangeGrob(foo, tab, layout_matrix = foo_lay, respect = TRUE)
        return(foo)
      })
    })
    if(create_pdf || !pdf_dev) {
      if(display_progress) {
        pb_pdf = newPB(session = dots$session, title = title_progress, label = "writing to pdf (no update but file is being processed)", min = 0, max = gl, initial = 0, style = 3)
        on.exit(endPB(pb_pdf), add = TRUE)
      }
      if(onepage) {
        if(pdf_dev) pdf(file=export_to_pdf, width=3*max(lay$x)*2.54, height=3*max(lay$y)*2.54, 
                        family="serif", onefile=TRUE, pagecentre=TRUE, useDingbats=FALSE)
        # TODO add a progress bar
        # for(i in 1:gl) {
        #   pos = matrix(NA, ncol = max(lay$x), nrow = max(lay$y))
        #   pos[lay_mat == lay$N[i]] <- i
        #   print(pos)
        #   grid.arrange(grobs = graphs[lay$N[1]], top = title_progress, newpage = TRUE, layout_matrix = lay_mat, as.table = FALSE)
        #   if(i == 1) {
        #     plot(graphs[lay$N[i]])
        #     grid.arrange(grobs = graphs[lay$N[i]], top = title_progress, newpage = TRUE, layout_matrix = pos, as.table = FALSE)
        #   } else {
        #     grid.arrange(grobs = graphs[lay$N[i]], newpage = FALSE, layout_matrix = pos, as.table = FALSE)
        #   }
        #   if(display_progress) {
        #     setPB(pb_pdf, value = i, title = title_progress, label = "writing to pdf")
        #   }
        # }
        grid.arrange(grobs = graphs[lay$N], top = title_progress, newpage = TRUE, layout_matrix = lay_mat, as.table = FALSE)
      } else {
        if(pdf_dev) pdf(file=export_to_pdf, paper = "a4", onefile = TRUE, pagecentre = TRUE, useDingbats = FALSE, family = "serif")
        for(i in 1:gl) {
          grid.arrange(graphs[[i]], top = title_progress, newpage = TRUE)
          if(display_progress) setPB(pb_pdf, value = i, title = title_progress, label = "writing to pdf")
        }
      }
    }
    if(pdf_dev) {
      message(paste0("\n######################\n",ifelse(length(write_to)==2, "files have", "file has"), " been successfully ", ifelse(overwritten, "overwritten", "exported"),"\n"))
      return(invisible(write_to))
    }
  }, error = function(e) {
    if(pdf_dev) message(paste0(ifelse(length(write_to)==2, "files have", "file has"), " been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n"))
    stop(e$message, call. = FALSE)
  },
  finally = {
    if(display_progress) endPB(pb_gr)
    if(pdf_dev) while(!all(dv == dev.cur())) {
      dev.off(which = rev(dev.cur())[1])
    }
  })
}

#' @title Graphical and Statistic Report Display
#' @description
#' Displays report from `IFC_data` object.
#' @param ... arguments to pass to \code{\link{ExportToReport}} where write_to is inoperative.
#' @return It invisibly returns NULL
#' @keywords internal
DisplayReport = function( ...) {
  ExportToReport(...)
}

#' @title Batch Generation of Graphical and Statistic Report
#' @description
#' Batch creates graphical an statistical report.
#' @param fileName path to file(s).
#' @param selection indices of desired graphs. It can be provided as an integer vector or as a matrix.
#' In such case, the layout of the matrix will reflect the layout of the extracted graphs for each single 'fileName'.
#' NA value will result in an empty place.
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
#' @param color_mode Whether to extract colors from obj in white or black mode. Default is 'white'.
#' @param add_key whether to draw a 'global' key under title or in the first 'panel' or 'both'. Default is 'panel'.\cr
#' Accepted values are either: FALSE, 'panel', 'global', 'both' or c('panel', 'global').\cr
#' Note that it only applies when display is seen as overlaying populations.
#' @param precision when graphs is a 2D scatter with population overlay, this argument controls amount of information displayed. Default is "light".\cr
#' -"light", the default, will only display points of same coordinates that are amoung the other layers.\cr
#' -"full" will display all the layers.
#' @param trunc_labels maximum number of characters to display for labels. Default is 38.
#' @param trans transformation function for density graphs. If missing the default, the BasePop[[1]]$densitytrans, if any, will be retrieved, otherwise asinh will be used.
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
                        trunc_labels=38, trans=asinh, bin, viewport="ideas",                                  # parameters to pass to plotGraph
                        display_progress=TRUE, ...) {
  dots = list(...)
  create_pdf = FALSE
  create_csv = FALSE
  
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
  # check
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
  byrow = as.logical(byrow); assert(byrow, len=1, alw=c(TRUE,FALSE))
  times = na.omit(as.integer(times)); assert(times, len = 1)
  if(missing(selection)) stop("'selection' can't be missing")
  if(missing(bin)) {
    bin = NULL
  } else {
    bin=na.omit(as.integer(bin)); assert(bin, len=1, typ="integer")
  }
  plotGraph_args=list(draw=FALSE, color_mode=color_mode, add_key=add_key,
                      precision=precision, trunc_labels=trunc_labels,
                      viewport=viewport)
  if(length(bin) != 0) plotGraph_args = c(plotGraph_args, list(bin=bin))
  if(!missing(trans)) plotGraph_args = c(plotGraph_args, list(trans=trans))
  
  if(missing(write_to)) stop("'write_to' can't be missing")
  assert(write_to, typ="character")
  file_extension = getFileExt(write_to)
  assert(file_extension, alw=c("pdf", "csv"))
  if(any(duplicated(file_extension))) stop("'write_to' has to be only one .pdf and / or only one .csv")
  export_to_pdf = NULL
  export_to_csv = NULL
  splitf_obj = splitf(fileName[1])
  if(any(file_extension%in%"pdf")) {
    create_pdf = TRUE
    splitp_obj_pdf = splitp(write_to[file_extension=="pdf"])
    if(any(splitp_obj_pdf$channel > 0)) message("'write_to' (pdf part) has %c argument but channel information can't be retrieved with ExportToReport()")
    if(any(splitp_obj_pdf$object > 0)) message("'write_to' (pdf part) has %o argument but channel information can't be retrieved with ExportToReport()")
    export_to_pdf = formatn(splitp_obj_pdf, splitf_obj)
  }
  if(any(file_extension%in%"csv")) {
    create_csv = TRUE
    splitp_obj_csv = splitp(write_to[file_extension=="csv"])
    if(any(splitp_obj_csv$channel > 0)) message("'write_to', (csv part) has %c argument but channel information can't be retrieved with ExportToReport()")
    if(any(splitp_obj_csv$object > 0)) message("'write_to' (csv part) has %o argument but channel information can't be retrieved with ExportToReport()")
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
  
  n_tiles = length(fileName)
  s_tiles = seq(from = 1, to = n_tiles, by = 1)
  lay = split(seq(from = 1, to = n_tiles, by = 1), ceiling(seq(from = 1, to = n_tiles, by = 1) / times) )
  if(length(lay) > 1) lay[[length(lay)]] = c(lay[[length(lay)]], rep(NA, times - length(lay[[length(lay)]])))
  lay = do.call(what = ifelse(byrow, "cbind", "rbind"), args = lay)
  L = length(selection)
  if(byrow) { LX = 1; LY = L } else { LX = L; LY = 1 }
  if(inherits(x = selection, what = "matrix")) {
    LX = ncol(selection)
    LY = ncol(selection)
  }
  apply_gating = FALSE
  if(!missing(gating)) apply_gating = TRUE

  # backup last state of graphic device
  dv = dev.cur()
  if(display_progress) pb = newPB(session = dots$session, min = 0, max = length(fileName), initial = 0, style = 3)
  tryCatch({
    grobs = lapply(1:length(fileName), FUN = function(i_file) {
      if(display_progress) setPB(pb = pb, value = i_file, title = "creating batch report", label = paste0("extracting ",basename(fileName[i_file])))
      obj = readIFC(fileName=fileName[i_file], display_progress=FALSE, ...)
      if(apply_gating) obj = applyGatingStrategy(obj=obj, gating=gating, display_progress=FALSE)
      G = obj$graphs[c(selection)]
      plotGraph_args = c(plotGraph_args, list(obj=obj))
      file_lay = matrix(1:L, ncol = LX, nrow = LY)
      file_grob = lapply(1:L, FUN = function(i) {
        if(length(G[[i]]) == 0) return(grob())
        g = try(do.call(what = plotGraph, args = c(plotGraph_args, list(graph = G[[i]]))))
        if(inherits(x = g, what = "try-error")) {
          foo = arrangeGrob(grid.text(label = paste0("Error: ", attr(x = g, which = "condition")$message), gp=gpar(col="red"), draw = FALSE),
                            top = textGrob(paste0("\n",G[[i]]$title), gp = gpar(fontsize = 8, font=2, lineheight=0.5)))
          if(G[[i]]$type=="histogram") {
            stats = matrix(NA, nrow = 1, ncol = 8)
            colnames(stats) = c("count","perc",
                                "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.")
          } else {
            stats = matrix(NA, nrow = 1, ncol = 14)
            colnames(stats) = c("count","perc",
                                "x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.",
                                "y-Min.","y-1st Qu.","y-Median","y-Mean","y-3rd Qu.","y-Max.")
          }
          rownames(stats) = paste0("Error: ", G[[i]]$title)
        } else {
          foo = g$plot
          stats = g$stats
        }
        if(create_csv) {
          write.table(x=basename(fileName[i_file]), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          write.table(x=rbind(c(G[[i]]$f1,"x")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          if(G[[i]]$type!="histogram") write.table(x=rbind(c(G[[i]]$f2,"y")), file=export_to_csv, append=TRUE, sep=",", col.names = FALSE, row.names = FALSE)
          write.table(x=stats, file=export_to_csv, append=TRUE, sep=",", col.names = FALSE)
        }
        stats = stats[,!grepl("Qu", colnames(stats)),drop=FALSE]
        foo_lay = matrix(c(rep(1,9), c(NA,2,NA)), nrow=4, ncol=3, byrow = TRUE)
        if(nrow(stats)> 7) {
          stats = stats[1:8,]
          stats[8, ] = "..."
        }
        foo$vp <- viewport(x=0.5, y=0.5)
        tab = tableGrob(format(stats, scientific=FALSE, digits=5), theme = ttheme_default(base_size=4, base_family = "serif"))
        tab$vp <- viewport(x=0.5, y=unit(1,"npc") - 0.5*sum(tab$heights))
        foo = arrangeGrob(foo, tab, layout_matrix=foo_lay, respect=TRUE)
        return(foo)
      })
      bar = grobTree(arrangeGrob(grobs=file_grob, ncol=LX, nrow=LY,
                                 top=textGrob(paste0("\n\n\n",basename(fileName[i_file]),"\n"), 
                                              gp=gpar(fontsize=14, font=2, col="skyblue4", lineheight=0.5))),
                     rectGrob(.5,.5,width=unit(.99,"npc"), height=unit(0.99,"npc"),
                              gp=gpar(lwd=2, col="skyblue4", fill=NA)))
    })
    if(create_pdf) {
      args = list(newpage=TRUE, layout_matrix=lay, as.table=FALSE)
      if(!missing(main)) args=c(args, top=main)
      pdf(file=write_to, height=3*nrow(lay)*LY*2.54, width=3*ncol(lay)*LX*2.54,
          family="serif", onefile=TRUE, pagecentre=TRUE, useDingbats=FALSE)
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
    if(display_progress) endPB(pb)
    while(!all(dv == dev.cur())) {
      dev.off(which = rev(dev.cur())[1])
    }
  })
}
