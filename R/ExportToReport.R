#' @title Graphical and Statistic Report Generation
#' @description
#' Generates report from IFC_data object.
#' @param obj an IFC_data object extracted with features extracted.
#' @param selection when provided, indices of desired graphs.
#' In such case onepage parameter is set to FALSE.
#' Note that indices are read from left to right, from top to bottom. 
#' @param export_to Pattern used to export file(s).
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
#' - 'light', the default, will only display points of same coordinates that are amoung the other layers.\cr
#' - 'full' will display all the layers.
#' @param trunc_labels maximum number of characters to display for labels. Default is 38.
#' @param trans transformation function for density graphs. Default is asinh.
#' @param bin default number of bin used for histogram. Default is missing.
#' @param viewport Either "ideas", "data" or "max" defining limits used for the graph. Default is "ideas".\cr
#' "ideas" will use same limits as the one defined in ideas.\cr
#' "data" will use data to define limits.\cr
#' "max" will use data and regions drawn to define limits.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param ... Other parameters to be passed.
#' @details depending on 'export_to', function will create .pdf and/or .csv file(s) report with according to graphs found in 'obj'.\cr
#' - csv file if created will contain "Min.","1st Qu.","Median","Mean","3rd Qu.","Max." for each graph found for x and y (if not histogram) for drawn populations and regions.\cr
#' - pdf file if created will contain graphs and to a certain extent some stats "Min.", "Median", "Mean", "Max." (no more than 7 rows).\cr
#' Note that only graphs will be exported (no images, features values, population stats, ...) in the same layout they were created and without sizing.
#' @examples
# #' \dontrun{
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
#'                    export_to = paste0(tmp, "\\test.pdf"), overwrite = TRUE)
#'   }
#' } else {
#'   message(sprintf('Please type `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
# #' }
#' @return It invisibly returns full path of exported .pdf and/or .csv file(s).
#' @export
ExportToReport = function(obj, selection, export_to, overwrite=FALSE, onepage=TRUE,
                          color_mode=c("white","black")[1], add_key="panel", precision=c("light","full")[1],    # parameters to pass to plotGraph
                          trunc_labels=38, trans=asinh, bin, viewport="ideas",                                  # parameters to pass to plotGraph
                          display_progress=TRUE, ...) {
  dots = list(...)
  # backup last state of device ask newpage and set to FALSE
  old_ask <- devAskNewPage(ask = FALSE)
  on.exit(devAskNewPage(ask = old_ask))
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  
  # check mandatory param
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class IFC_data")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  if(missing(export_to)) stop("'export_to' can't be missing")
  assert(export_to, typ = "character")
  file_extension = getFileExt(export_to)
  assert(file_extension, alw = c("pdf", "csv"))
  if(any(duplicated(file_extension))) stop("'export_to' has to be only one .pdf and / or only one .csv")

  splitf_obj = splitf(obj$fileName)
  create_pdf = FALSE
  create_csv = FALSE
  export_to_pdf = NULL
  export_to_csv = NULL
  if(any(file_extension%in%"pdf")) {
    create_pdf = TRUE
    splitp_obj_pdf = splitp(export_to[file_extension=="pdf"])
    if(any(splitp_obj_pdf$channel > 0)) message("'export_to' (pdf part) has %c argument but channel information can't be retrieved with ExportToReport()")
    if(any(splitp_obj_pdf$object > 0)) message("'export_to' (pdf part) has %o argument but channel information can't be retrieved with ExportToReport()")
    export_to_pdf = formatn(splitp_obj_pdf, splitf_obj)
  }
  if(any(file_extension%in%"csv")) {
    create_csv = TRUE
    splitp_obj_csv = splitp(export_to[file_extension=="csv"])
    if(any(splitp_obj_csv$channel > 0)) message("'export_to', (csv part) has %c argument but channel information can't be retrieved with ExportToReport()")
    if(any(splitp_obj_csv$object > 0)) message("'export_to' (csv part) has %o argument but channel information can't be retrieved with ExportToReport()")
    export_to_csv = formatn(splitp_obj_csv, splitf_obj)
  }
  export_to = c(export_to_pdf, export_to_csv)
  
  assert(color_mode, len=1, alw=c("white","black"))
  assert(precision, len=1, alw=c("light","full"))
  if(!all(add_key%in%c("panel","global","both",FALSE))) stop("Accepted values for add_key are either: FALSE, 'panel', 'global', 'both' or c('panel', 'global')")
  assert(onepage, len=1, alw=c(TRUE,FALSE))
  assert(overwrite, len=1, alw=c(TRUE,FALSE))
  trunc_labels=na.omit(as.integer(trunc_labels));trunc_labels=trunc_labels[trunc_labels>=0]
  assert(trunc_labels, len=1, typ="integer")
  if(missing(bin)) {
    bin = NULL
  } else {
    bin=na.omit(as.integer(bin)); assert(bin, len=1, typ="integer")
  }
  display_progress = as.logical(display_progress); assert(display_progress, len=1, alw=c(TRUE,FALSE))
  if(display_progress) {
    if(.Platform$OS.type == "windows") {
      pb_fun = setWinProgressBar
    } else {
      pb_fun = setTxtProgressBar
    }
  }
  overwritten = FALSE
  tmp = file.exists(export_to)
  if(any(tmp)) {
    if(!overwrite) stop(paste0("file ",paste0(export_to[tmp]," already exists"), collapse="\n"))
    if(create_pdf) if(file.exists(export_to_pdf)) export_to_pdf = normalizePath(export_to_pdf, winslash = "/")
    if(create_csv) if(file.exists(export_to_csv)) export_to_csv = normalizePath(export_to_csv, winslash = "/")
    overwritten = TRUE
  }
  if(create_pdf) {
    tryCatch({
      if(!dir.exists(dirname(export_to_pdf))) {
        if(!dir.create(dirname(export_to_pdf), showWarnings = FALSE, recursive = TRUE)) {
          stop(paste0(export_to_pdf,"\ncan't be created: check dirname"), call. = FALSE)
        }
      }
      pdf(file=export_to_pdf)
    }, error = function(e) stop(paste0(export_to_pdf,"\ncan't be created: check name / currently opened ?"), call. = FALSE))
    dev.off(dev.cur())
    export_to_pdf = normalizePath(export_to_pdf, winslash = "/")
  }
  if(create_csv) {
    tryCatch({
      if(!dir.exists(dirname(export_to_csv))) {
        if(!dir.create(dirname(export_to_csv), showWarnings = FALSE, recursive = TRUE)) {
          stop(paste0(export_to_csv,"\ncan't be created: check dirname"), call. = FALSE)
        }
      }
      write.table(x=rbind(c("pop","count","perc","x-Min.","x-1st Qu.","x-Median","x-Mean","x-3rd Qu.","x-Max.",
                                   "y-Min.","y-1st Qu.","y-Median","y-Mean","y-3rd Qu.","y-Max.")),
                         sep=",", row.names = FALSE, col.names = FALSE, file=export_to_csv)
      }, error = function(e) stop(paste0(export_to_csv,"\ncan't be created: check name / currently opened ?"), call. = FALSE))
    export_to_csv = normalizePath(export_to_csv, winslash = "/")
  }
  export_to = c(export_to_pdf,export_to_csv)
  message(paste0(ifelse(length(export_to)==2, "files", "file")," will be exported in :\n",paste0(normalizePath(dirname(export_to), winslash = "/"),collapse="\n")))
  tryCatch({
    # shortcuts
    G = obj$graphs
    P = obj$pops
    R = obj$regions
    # defines layout
    lay=lapply(G, FUN=function(g) with(g, c(x=xlocation,y=ylocation)))
    lay=as.data.frame(do.call(rbind, lay))
    lay=by(lay, lay$y, FUN=function(d) {
      d=d[order(d$x), ]
      d$x=seq_along(d$x)
      cbind("N"=as.numeric(row.names(d)), d, stringsAsFactors=FALSE)
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
    if(display_progress) {
      if(.Platform$OS.type == "windows") {
        pb_gr = winProgressBar(title = basename(obj$fileName), min = 0, max = 100, initial=1, label="done")
      } else {
        pb_gr = txtProgressBar(title = basename(obj$fileName), min = 0, max = 100, initial=1, label="done", style = 3)
      }
      on.exit(close(pb_gr), add = TRUE)
    }
    suppressWarnings({
      graphs = lapply(1:gl, FUN=function(i) {
        if(length(bin) == 0) {
          g = plotGraph(obj, G[[i]], draw=FALSE, color_mode=color_mode, add_key=add_key,
                        precision=precision, trunc_labels=trunc_labels, trans=trans, viewport=viewport)
        } else {
          g = plotGraph(obj, G[[i]], draw=FALSE, color_mode=color_mode, add_key=add_key,
                        precision=precision, trunc_labels=trunc_labels, trans=trans, viewport=viewport, bin=bin)
        }
        foo = g$plot
        stats = g$stats
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
        if(display_progress) {
          k=i/gl*100
          pb_fun(pb = pb_gr, value = k, label = paste0("Computing ",ifelse(create_pdf,"graphs and ",""),"stats ", sprintf("%.0f%%", k)))
        }
        return(foo)
      })
    })
    if(create_pdf) {
      if(display_progress) {
        if(.Platform$OS.type == "windows") {
          pb_pdf = winProgressBar(title = basename(obj$fileName), label = "Writing to pdf\nPlease wait...", min = 0, max = 100, initial=1)
        } else {
          pb_pdf = txtProgressBar(title = basename(obj$fileName), label = "Writing to pdf. Please wait...", min = 0, max = 100, initial=1, style = 3)
        }
        on.exit(close(pb_pdf), add = TRUE)
      }
      if(onepage) {
        pdf(file=export_to_pdf, width = 3*max(lay$x)*2.54, height = 3*max(lay$y)*2.54, 
            family = "serif", onefile = TRUE, pagecentre = TRUE, useDingbats = FALSE)
        on.exit(dev.off(which = dev.cur()), add = TRUE)
        # TODO add a progress bar
        # for(i in 1:gl) {
        #   pos = matrix(NA, ncol = max(lay$x), nrow = max(lay$y))
        #   pos[lay_mat == lay$N[i]] <- i
        #   print(pos)
        #   grid.arrange(grobs = graphs[lay$N[1]], top = basename(obj$fileName), newpage = TRUE, layout_matrix = lay_mat, as.table = FALSE)
        #   if(i == 1) {
        #     plot(graphs[lay$N[i]])
        #     grid.arrange(grobs = graphs[lay$N[i]], top = basename(obj$fileName), newpage = TRUE, layout_matrix = pos, as.table = FALSE)
        #   } else {
        #     grid.arrange(grobs = graphs[lay$N[i]], newpage = FALSE, layout_matrix = pos, as.table = FALSE)
        #   }
        #   if(display_progress) {
        #     k=i/gl*100
        #     pb_fun(pb_pdf, value = k, label = sprintf("Writing to pdf %.0f%%", k))
        #   }
        # }
        grid.arrange(grobs = graphs[lay$N], top = basename(obj$fileName), newpage = TRUE, layout_matrix = lay_mat, as.table = FALSE)
      } else {
        pdf(file=export_to_pdf, paper = "a4", onefile = TRUE, pagecentre = TRUE, useDingbats = FALSE, family = "serif")
        on.exit(dev.off(which = dev.cur()), add = TRUE)
        for(i in 1:gl) {
          grid.arrange(graphs[[i]], top=basename(obj$fileName), newpage = TRUE) #, respect = TRUE)
          if(display_progress) {
            k=i/gl*100
            pb_fun(pb_pdf, value = k, label = sprintf("Writing to pdf %.0f%%", k))
          }
        }
      }
    }
  }, error = function(e) {
    message(paste0(ifelse(length(export_to)==2, "files have", "file has"), " been incompletely ", ifelse(overwritten, "overwritten", "exported"), "\n"))
    stop(e$message, call. = FALSE)
  })
  message(paste0("\n######################\n",ifelse(length(export_to)==2, "files have", "file has"), " been successfully ", ifelse(overwritten, "overwritten", "exported"),"\n"))
  return(invisible(export_to))
}
