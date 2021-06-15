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

#' @title Gating Strategy File Reader
#' @description
#' Extracts Gating Strategy from files.
#' @param fileName path to file. It should be a .ast, .cif, .daf, .ist, .rif or .xml file.
#' @return A named list of class `IFC_gating`, whose members are:\cr
#' -graphs, a list of graphical elements found,\cr
#' -pops, a list describing populations found,\cr
#' -regions, a list describing how regions are defined.
#' @param ... other arguments to be passed.
#' @export
readGatingStrategy <- function(fileName, ...) {
  dots=list(...)
  if(missing(fileName)) stop("'fileName' can't be missing")
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = c("daf","cif","rif","xml","ast","ist"))
  if(!file.exists(fileName)) stop(paste("can't find",fileName,sep=" "))
  title_progress = basename(fileName)
  
  fileName = normalizePath(fileName, winslash = "/", mustWork = FALSE)
  if(file_extension == "xml") return(readGatingML(fileName, ...))
  if(file_extension %in% c("daf", "ast")) {
    assay = switch(file_extension,
                   daf = "/Assay",
                   ast = "/AssayTemplate")
    toskip=cpp_scanFirst(fname = fileName, target = paste0("<",assay,">"), start = 0, end = 0)
    if(toskip==0) stop(paste0(fileName, "\ndoes not seem to be well formatted: <",assay,"> not found"))
    toskip = toskip + nchar(paste0("<",assay,">")) - 1
    tmp=read_xml(readBin(con = fileName, what = "raw", n = toskip), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  } else {
    if(file_extension %in% c("cif", "rif")) {
      IFD = getIFD(fileName = fileName, offsets = "first", trunc_bytes = 8, force_trunc = FALSE, verbose = FALSE, verbosity = 0, bypass = TRUE)
      tmp=read_xml(as_list(xml_find_first(read_xml(getFullTag(IFD = IFD, which = 1, "33027"),
                                          options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN")), "//Imaging//DafFile"))[[1]],
                   options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
    } else { # for ist
      tmp=read_xml(as_list(xml_find_first(read_xml(fileName, 
                                                   options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN")), "//Imaging//DafFile"))[[1]],
                   options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
    }
    assay = "/AssayTemplate"
  }

  ##### extracts description
  # assay_attr = xml_attrs(xml_find_all(tmp, paste0("/",assay)))
  # description=list("Assay" = assay_attr,
  #                  "FCS"=xml_attrs(xml_find_all(tmp, "//FCS")),
  #                  "SOD"=xml_attrs(xml_find_all(tmp, "//SOD")))
  # description=lapply(description, FUN=function(x) {as.data.frame(do.call(what="rbind", x), stringsAsFactors=FALSE)})
  # if(length(description$FCS)==0) {
  #   description$ID = description$SOD
  #   is_fcs = FALSE
  # } else {
  #   description$ID = description$FCS
  #   is_fcs = TRUE
  # }
  # obj_count = as.integer(description$ID$objcount)
  
  ##### extracts graphs information
  plots=lapply(xml_attrs(xml_find_all(tmp, "//Graph")), FUN=function(x) as.list(x))
  if(length(plots)!=0) {
    # plots=mapply(plots, FUN=c, SIMPLIFY = FALSE)
    plots_tmp=lapply(plots, FUN=function(plot) {
      pat=paste0("//Graph[@xlocation='",plot$xlocation,"'][@ylocation='",plot$ylocation,"']")
      sapply(c("Legend","BasePop","GraphRegion","ShownPop"), FUN=function(i_subnode){
        lapply(xml_attrs(xml_find_all(tmp, paste(pat,i_subnode,sep="//"))), FUN=function(x) as.list(x))
      })
    })
    plots=mapply(plots, plots_tmp, FUN = append, SIMPLIFY = FALSE)
    plots_tmp=c("xlocation","ylocation","scaletype","xmin","xmax","ymin","ymax","axislabelsfontsize","axistickmarklabelsfontsize",
                "graphtitlefontsize","regionlabelsfontsize","bincount","histogramsmoothingfactor","xsize","ysize","splitterdistance")
    plots=lapply(plots, FUN=function(x) {replace(x, plots_tmp, lapply(x[plots_tmp], as.numeric))})
    plot_order=sapply(plots, FUN=function(i_plot) as.numeric(i_plot[c("xlocation", "ylocation")]))
    plots=plots[order(unlist(plot_order[1,]),unlist(plot_order[2,]))]
    plots=plots[order(unlist(plot_order[2,]))]
    rm(list=c("plots_tmp", "plot_order"))
  }
  class(plots) <- "IFC_graphs"
  
  ##### extracts regions information
  regions=lapply(xml_attrs(xml_find_all(tmp, "//Region")), FUN=function(x) as.list(x))
  if(length(regions) != 0) {
    names(regions)=lapply(regions, FUN=function(x) x$label)
    # regions=mapply(regions, FUN=c, SIMPLIFY = FALSE)
    regions=lapply(regions, FUN=function(x) {
      N = names(x)
      ##### changes unknown color names in regions
      if("color" %in% N) x[["color"]] <- map_color(x[["color"]])
      if("lightcolor" %in% N) x[["lightcolor"]] <- map_color(x[["lightcolor"]])
      ##### convert label position to numeric
      if("cx" %in% N) x[["cx"]] <- as.numeric(x["cx"])
      if("cy" %in% N) x[["cy"]] <- as.numeric(x["cy"])
      x
    })
    regions_tmp=lapply(regions, FUN=function(i_region) {
      pat=paste0("//Region[@label='",i_region$label,"']//axy")
      axy=do.call(cbind, args = xml_attrs(xml_find_all(tmp, pat)))
      list(x=as.numeric(axy["x",]), y=as.numeric(axy["y",]))
    })
    regions=mapply(FUN = append, regions, regions_tmp, SIMPLIFY = FALSE)
    rm(regions_tmp)
  }
  class(regions) <- "IFC_regions"
  
  ##### extracts populations information
  pops=lapply(xml_attrs(xml_find_all(tmp, "//Pop")), FUN=function(x) as.list(x))
  if(length(pops)>0) {
    names(pops)=lapply(pops, FUN=function(x) x$name)
    if((length(dots$display_progress) != 0) && dots$display_progress) {
      pb_pops = newPB(session = dots$session, min = 0, max = length(pops), initial = 0, style = 3)
      tryCatch({
        pops_=lapply(1:length(pops), FUN=function(i_pop) {
          setPB(pb_pops, value = i_pop, title = title_progress, label = "extracting tagged population objects")
          pat=paste0("//Pop[@name='",pops[[i_pop]]$name,"']//ob")
          list(obj=as.integer(unlist(xml_attrs(xml_find_all(tmp, pat)))))
        })
      }, error = function(e) {
        stop(e$message)
      }, finally = endPB(pb_pops))
    } else {
      pops_=lapply(1:length(pops), FUN=function(i_pop) {
        pat=paste0("//Pop[@name='",pops[[i_pop]]$name,"']//ob")
        list(obj=as.integer(unlist(xml_attrs(xml_find_all(tmp, pat)))))
      })
    }
    pops=mapply(FUN = append, pops, pops_, SIMPLIFY = FALSE)
    rm(pops_)
  }
  class(pops) <- "IFC_pops"
  
  ans = list("graphs"=plots, "pops"=pops, "regions"=regions)
  attr(ans, "class") <- c("IFC_gating")
  return(ans)
}
