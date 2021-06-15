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

#' @title GatingML Conversion
#' @description 
#' Helper to convert pops and regions to XML nodes in GatingML files.
#' @param pop a member of `IFC_pops` object.
#' @param reg a member of `IFC_regions` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a xml_node.
#' @keywords internal
toXML2_gating <- function(pop, reg, verbose = FALSE) {
  # region names conversion
  Nr = names(reg)
  Nr[Nr == "color"] <- "color1"
  Nr[Nr == "lightcolor"] <- "color2"
  names(reg) <- Nr
  reg = reg[!(Nr %in% c("ismarker", "doesnotoverride"))]
  
  # color conversion
  reg$color1 <- map_color(reg$color1, FALSE)
  reg$color2 <- map_color(reg$color2, FALSE)

  # nodes creation 
  # TODO add more gates ?
  ids_attrs = list("_ns_gating_ns_id" = pop$name, "_ns_gating_ns_parent_id" = pop$base)[1:ifelse(pop$base %in% c("","All"), 1, 2)]
  info_attrs1 = reg[!(names(reg) %in% c("x", "y"))]
  info_attrs2 = pop[c("name", "base", "color1", "color2", "style", "fx")]
  dim1_node = xml_new_node(name = "_ns_data-type_ns_fcs-dimension", attrs = list("_ns_data-type_ns_name"=pop$fx))
  if(reg$type == "line") {
    info_attrs1 = c(info_attrs1, y=reg$y[1])
  } else {
    dim2_node = xml_new_node("_ns_data-type_ns_fcs-dimension", attrs = list("_ns_data-type_ns_name"=pop$fy))
    info_attrs2 = c(info_attrs2, fy=pop$fy)
  }
  
  info_node = xml_new_node(name = "_ns_data-type_ns_custom_info",
               .children = list(xml_new_node(name = "region", attrs = info_attrs1),
                                xml_new_node(name = "population", attrs = info_attrs2)))
  switch(reg$type,
         "line" ={
           xml_new_node(name = "_ns_gating_ns_RectangleGate",
                        attrs = ids_attrs,
                        .children = c(list(info_node,
                                           xml_new_node(name = "_ns_gating_ns_dimension", attrs = list("_ns_gating_ns_min"=num_to_string(reg[["x"]][1]), "_ns_gating_ns_max"=num_to_string(reg[["x"]][2]), "_ns_gating_ns_compensation-ref"="uncompensated"),
                                                        .children = dim1_node)
                        )))
         },
         "rect" ={
           xml_new_node(name = "_ns_gating_ns_RectangleGate",
                        attrs = ids_attrs,
                        .children = c(list(info_node,
                                           xml_new_node(name = "_ns_gating_ns_dimension", attrs = list("_ns_gating_ns_min"=num_to_string(reg[["x"]][1]), "_ns_gating_ns_max"=num_to_string(reg[["x"]][2]), "_ns_gating_ns_compensation-ref"="uncompensated"),
                                                        .children = dim1_node),
                                           xml_new_node(name = "_ns_gating_ns_dimension", attrs = list("_ns_gating_ns_min"=num_to_string(reg[["y"]][1]), "_ns_gating_ns_max"=num_to_string(reg[["y"]][2]), "_ns_gating_ns_compensation-ref"="uncompensated"),
                                                        .children = dim2_node)
                        )))
         },
         "poly" ={
           xml_new_node(name = "_ns_gating_ns_PolygonGate",
                        attrs = ids_attrs,
                        .children = c(list(info_node,
                                           xml_new_node(name = "_ns_gating_ns_dimension", attrs = list("_ns_gating_ns_compensation-ref"="uncompensated"), .children = dim1_node),
                                           xml_new_node(name = "_ns_gating_ns_dimension", attrs = list("_ns_gating_ns_compensation-ref"="uncompensated"), .children = dim2_node)),
                                      lapply(1:length(reg[["x"]]), FUN = function(i_coord) {
                                        xml_new_node(name = "_ns_gating_ns_vertex",
                                                     .children = list(xml_new_node(name = "_ns_gating_ns_coordinate", attrs = list("_ns_data-type_ns_value" = num_to_string(reg[["x"]][i_coord]))),
                                                                      xml_new_node(name = "_ns_gating_ns_coordinate", attrs = list("_ns_data-type_ns_value" = num_to_string(reg[["x"]][i_coord])))))
                                      })
                        ))
         },
         "oval" ={
           center_x = mean(reg[["x"]])
           center_y = mean(reg[["y"]])
           xml_new_node(name = "_ns_gating_ns_EllipsoidGate",
                        attrs = ids_attrs,
                        .children = c(list(info_node,
                                           xml_new_node(name = "_ns_gating_ns_dimension", attrs = list("_ns_gating_ns_compensation-ref"="uncompensated"), .children = dim1_node),
                                           xml_new_node(name = "_ns_gating_ns_dimension", attrs = list("_ns_gating_ns_compensation-ref"="uncompensated"), .children = dim2_node),
                                           xml_new_node(name = "_ns_gating_ns_mean",
                                                        .children = list(xml_new_node(name = "_ns_gating_ns_coordinate", attrs = list("_ns_data-type_ns_value" = num_to_string(center_x))),
                                                                         xml_new_node(name = "_ns_gating_ns_coordinate", attrs = list("_ns_data-type_ns_value" = num_to_string(center_y))))),
                                           xml_new_node(name = "_ns_gating_ns_covarianceMatrix",
                                                        .children = list(xml_new_node(name = "_ns_gating_ns_row",
                                                                                      .children = list(xml_new_node(name = "_ns_gating_ns_entry", attrs = list("_ns_data-type_ns_value"=num_to_string((center_x-max(reg[["x"]]))^2))),  # remains to be checked
                                                                                                       xml_new_node(name = "_ns_gating_ns_entry", attrs = list("_ns_data-type_ns_value"=num_to_string((center_y-max(reg[["y"]]))^2))))),# remains to be checked
                                                                         xml_new_node(name = "_ns_gating_ns_row",
                                                                                      .children = list(xml_new_node(name = "_ns_gating_ns_entry", attrs = list("_ns_data-type_ns_value"=num_to_string((center_x-max(reg[["x"]]))^2))),  # remains to be checked
                                                                                                       xml_new_node(name = "_ns_gating_ns_entry", attrs = list("_ns_data-type_ns_value"=num_to_string((center_y-max(reg[["y"]]))^2))))) # remains to be checked
                                                        )),
                                           xml_new_node(name = "_ns_gating_ns_distanceSquare", attrs = list("_ns_data-type_ns_value"="1"))
                        )))
         }
  )
}

#' @title IFC_graphs GatingML Conversion
#' @description 
#' Helper to convert graphs (`IFC_graphs` object) to XML nodes in GatingML files.
#' @param graphs an `IFC_graphs` object.
#' @return a xml_node.
#' @keywords internal
toXML2_graphs_gs = function(graphs) {
  if(length(graphs)==0) return(list())
  graphs = lapply(graphs, FUN=function(i) {
    i$GraphRegion = lapply(i$GraphRegion, FUN = function(g) g[!grepl("def", names(g))]) # it is mandatory to remove def
    if(typeof(i) %in% c("integer","double")) {
      return(num_to_string(i))
    } else {
      return(i)
    }
  })
  lapply(graphs, FUN=function(i_graph) {
    region = paste0(unlist(i_graph[["GraphRegion"]]), collapse = "|")
    if(region == "") region = NULL
    overlay = paste0(unlist(i_graph[["ShownPop"]]), collapse = "|")
    if(overlay == "") overlay = NULL
    xml_new_node(name = "plot", attrs = i_graph[!grepl("Legend|BasePop|GraphRegion|ShownPop", names(i_graph)) & sapply(i_graph, FUN=function(ele) length(ele)!=0)],
                 .children = c(lapply(i_graph[["Legend"]], FUN=function(i) xml_new_node(name = "legend", attrs = i)),
                               lapply(i_graph[["BasePop"]], FUN=function(i) xml_new_node(name = "basepop", attrs = i)),
                               lapply(region, FUN=function(i) xml_new_node(name = "region", attrs = list(displayed = i))),
                               lapply(overlay, FUN=function(i) xml_new_node(name = "overlay", attrs = list(displayed = i)))))
  })
}

#' @title GatingML File Reader
#' @description
#' Extracts GatingML from file.
#' @param fileName path to file. It should be a .xml file.
#' @param ... other arguments to be passed.
#' @details reading GatingML files is in development and partly implemented.
#' For the moment only files generated with IFC package can be read.
#' @return A named list of class `IFC_gating`, whose members are:\cr
#' -graphs, a list of graphical elements found,\cr
#' -pops, a list describing populations found,\cr
#' -regions, a list describing how regions are defined.
#' @keywords internal
readGatingML <- function(fileName, ...) {
  dots=list(...)
  if(missing(fileName)) stop("'fileName' can't be missing")
  file_extension = getFileExt(fileName)
  assert(file_extension, len = 1, alw = "xml")
  if(!file.exists(fileName)) stop(paste("can't find",fileName,sep=" "))
  title_progress = basename(fileName)
  
  fileName = normalizePath(fileName, winslash = "/", mustWork = FALSE)
  tmp=read_xml(fileName, options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
  ns = xml_ns(tmp)
  gating = suppressWarnings(xml_find_first(tmp, "/gating:Gating-ML", ns = ns))
  if(length(gating) == 0) stop(paste0(fileName, "\ndoes not seem to be well formatted:  </gating:Gating-ML> not found"))
  
  general_custom = suppressWarnings(xml_find_first(gating, "data-type:custom_info", ns = ns))
  # check if xml has been generated by IFC package
  is_ifc = suppressWarnings(xml_find_first(general_custom, "info/IFC"))
  
  if(length(is_ifc) == 0) {
    #### extracts gating information
    gates = sapply(c("gating:RectangleGate", "gating:EllipsoidGate", "gating:BooleanGate","gating:PolygonGate", "gating:QuadrantGate"), simplify = FALSE, xml_find_all, x = gating)
    
    
  } else {
    ##### extracts regions information
    regions=lapply(xml_attrs(xml_find_all(general_custom, "regions/region")), FUN=function(x) as.list(x))
    if(length(regions) != 0) {
      names(regions)=lapply(regions, FUN=function(x) x$label)
      regions=lapply(regions, FUN=function(x) {
        ##### modify regions names
        N = names(x)
        N[N == "color1"] <- "color"
        N[N == "color2"] <- "lightcolor"
        names(x) <- N
        if("color" %in% N) x[["color"]] <- map_color(x[["color"]])
        if("lightcolor" %in% N) x[["lightcolor"]] <- map_color(x[["lightcolor"]])
        if("cx" %in% N) x[["cx"]] <- as.numeric(x["cx"])
        if("cy" %in% N) x[["cy"]] <- as.numeric(x["cy"])
        if("x" %in% N) x[["x"]] <- as.numeric(strsplit(x[["x"]], split="|", fixed=TRUE)[[1]])
        if("y" %in% N) x[["y"]] <- as.numeric(strsplit(x[["y"]], split="|", fixed=TRUE)[[1]])
        x
      })
    }
    
    ##### extracts populations information
    pops=lapply(xml_attrs(xml_find_all(general_custom, "pops//pop")), FUN=function(x) as.list(x))
    if(length(pops)>0) {
      names(pops)=lapply(pops, FUN=function(x) x$name)
      pops = lapply(pops, FUN=function(pop) {
        ##### modify pop names
        N = names(pop)
        N[N == "color1"] <- "color"
        N[N == "color2"] <- "lightModeColor"
        names(pop) <- N
        if(pop$type == "T") pop$obj = as.integer(strsplit(pop$obj, split="|",fixed=TRUE)[[1]])
        return(pop)
      })
    }
    
    ##### extracts graphs information
    plots=lapply(xml_attrs(xml_find_all(general_custom, "plots//plot")), FUN=function(x) as.list(x))
    if(length(plots)!=0) {
      # plots=mapply(plots, FUN=c, SIMPLIFY = FALSE)
      plots_tmp=lapply(plots, FUN=function(plot) {
        pat=paste0("//plot[@xlocation='",plot$xlocation,"'][@ylocation='",plot$ylocation,"']")
        sapply(c("legend","basepop","region","overlay"), simplify = FALSE, FUN=function(i_subnode){
          lapply(xml_attrs(xml_find_all(general_custom, paste(pat,i_subnode,sep="/"))), FUN=function(x) as.list(x))
        })
      })
      plots=mapply(plots, plots_tmp, FUN = append, SIMPLIFY = FALSE)
      plots_tmp=c("xlocation","ylocation","scaletype","xmin","xmax","ymin","ymax","axislabelsfontsize","axistickmarklabelsfontsize",
                  "graphtitlefontsize","regionlabelsfontsize","bincount","histogramsmoothingfactor","xsize","ysize","splitterdistance")
      plots=lapply(plots, FUN=function(x) {replace(x, plots_tmp, lapply(x[plots_tmp], as.numeric))})
      plots_tmp=c("region","overlay")
      plot_order=sapply(plots, FUN=function(i_plot) as.numeric(i_plot[c("xlocation", "ylocation")]))
      plots=plots[order(unlist(plot_order[1,]),unlist(plot_order[2,]))]
      plots=plots[order(unlist(plot_order[2,]))]
      rm(list=c("plots_tmp", "plot_order"))
      # now we need to modify plot with names of regions and pops
      plots = lapply(plots, FUN = function(g) {
        ##### modify pops and regions names
        N = names(g)
        N[N == "basepop"] <- "BasePop"
        N[N == "legend"] <- "Legend"
        N[N == "overlay"] <- "ShowPop"
        N[N == "region"] <- "GraphRegion"
        names(g) <- N
        if(length(g$ShowPop) != 0) g$ShowPop = lapply(splitn(definition = g$ShowPop[[1]]$displayed, all_names = names(pops)), FUN = function(x) list(name=x))
        if(length(g$GraphRegion) != 0) g$GraphRegion = lapply(splitn(definition = g$GraphRegion[[1]]$displayed, all_names = names(regions)), FUN = function(r) list(name = r))
        # {
        #   g$GraphRegion = lapply(splitn(definition = g$GraphRegion[[1]]$displayed, all_names = names(regions)), FUN = function(r) {
        #     foo = sapply(pops, FUN = function(p) {
        #       bar = (p$type == "G") &&
        #         (p$region == r) &&
        #         (p$base %in% unique(unlist(lapply(g$BasePop, FUN = function(b) b$name)))) &&
        #         (g$f1 == p$fx)
        #       if(regions[[r]]$type != "line") bar = bar && (g$f2 == p$fy)
        #       return(bar)
        #     })
        #     return(list(name = r)) #, def = names(which(foo))))
        #   })
        # }
        return(g)
      })
    }
  }
  
  # add classes
  class(plots) <- "IFC_graphs"
  class(regions) <- "IFC_regions"
  class(pops) <- "IFC_pops"
  # returned ans
  ans = list("graphs"=plots, "pops"=pops, "regions"=regions)
  attr(ans, "class") <- c("IFC_gating")
  return(ans)
}
