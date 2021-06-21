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

#' @title Graphical Population GatingML Conversion to XML2
#' @description 
#' Helper to convert pops and regions to XML nodes in GatingML files.
#' @param pop a member of `IFC_pops` object.
#' @param reg a member of `IFC_regions` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @return a xml_node.
#' @keywords internal
toXML2_graphpop_gs <- function(pop, reg, verbose = FALSE) {
  # color conversion
  reg$colors = paste0(map_color(c(reg$color,reg$lightcolor), FALSE),collapse="|")
  # only keep what is need
  reg = reg[!(names(reg) %in% c("ismarker", "doesnotoverride","color","lightcolor"))]

  # nodes creation 
  # TODO add more gates ?
  ids_attrs = list("_ns_gating_ns_id" = pop$name, "_ns_gating_ns_parent_id" = pop$base)[1:ifelse(pop$base %in% c("","All"), 1, 2)]
  
  # write custom_info
  info_attrs = list("rtext"=reg$label,
                     "xtext"=reg$cx,
                     "ytext"=reg$cy,
                     "lineat"=reg$y[1],
                     "rcolors"=reg$colors,
                     "colors"=pop$colors,
                     "pch"=map_style(pop$style, toR = FALSE),
                     "xlog"=reg$xlogrange,
                     "ylog"=reg$ylogrange)
  
  # only add xtrans/ytrans if they are not empty
  if(length(reg$xtrans)!=0) info_attrs = c(info_attrs, list("xtrans"=reg$xtrans))
  if(length(reg$ytrans)!=0) info_attrs = c(info_attrs, list("ytrans"=reg$ytrans))
  
  dim1_node = xml_new_node(name = "_ns_data-type_ns_fcs-dimension", attrs = list("_ns_data-type_ns_name"=pop$fx))
  if(reg$type != "line") dim2_node = xml_new_node("_ns_data-type_ns_fcs-dimension", attrs = list("_ns_data-type_ns_name"=pop$fy))
  
  info_node = xml_new_node(name = "_ns_data-type_ns_custom_info", attrs = info_attrs)
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
                                                                      xml_new_node(name = "_ns_gating_ns_coordinate", attrs = list("_ns_data-type_ns_value" = num_to_string(reg[["y"]][i_coord])))))
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
                                                                                      .children = list(xml_new_node(name = "_ns_gating_ns_entry", attrs = list("_ns_data-type_ns_value"=num_to_string((center_x-reg[["x"]][1])^2))),  # remains to be checked
                                                                                                       xml_new_node(name = "_ns_gating_ns_entry", attrs = list("_ns_data-type_ns_value"=0)))),# remains to be checked
                                                                         xml_new_node(name = "_ns_gating_ns_row",
                                                                                      .children = list(xml_new_node(name = "_ns_gating_ns_entry", attrs = list("_ns_data-type_ns_value"=0)),  # remains to be checked
                                                                                                       xml_new_node(name = "_ns_gating_ns_entry", attrs = list("_ns_data-type_ns_value"=num_to_string((center_y-reg[["y"]][1])^2))))) # remains to be checked
                                                        )),
                                           xml_new_node(name = "_ns_gating_ns_distanceSquare", attrs = list("_ns_data-type_ns_value"="1"))
                        )))
         }
  )
}

#' @title Boolean Population GatingML Conversion to XML2
#' @description 
#' Helper to convert boolean population to XML nodes in GatingML files.
#' @param obj an `IFC_data` object.
#' @param pop a member of `IFC_pops` object.
#' @return a list of xml_node.
#' @keywords internal
toXML2_boolpop_gs <- function(obj, pop) {
  # recover definition
  pop_def <- pop$split
  pop_def[pop_def == "And"]="&"
  pop_def[pop_def == "Or"] ="|"
  pop_def[pop_def == "Not"]="!"
  pop_names = !(pop_def %in% c("&", "|", "!", "(", ")"))
  pop_def[pop_names] = paste0("`",pop_def[pop_names],"`")
  pop_def = str2lang(paste0(pop_def,collapse = " "))
  
  # recover arithmetic tree
  expand_tree = function(x) {
    branch = lapply(x, unlist)
    if(length(branch) <= 1) return(branch)
    if(inherits(x, "("))  branch = branch[-1]
    lapply(branch, expand_tree)
  }
  
  tree = expand_tree(pop_def)
  # inialize values for recursive loop
  children = c()
  ids = c()
  # group by operation
  op = lapply(rev(unlist(tree)), FUN = function(x) {
    ele <- as.character(unlist(x))
    if(ele %in% c("|", "&", "!")) {
      # an operator is encountered
      # we generate a random id for the population resulting of the operation
      id = random_name(special = NULL, alpha = NULL, forbidden = c(ids, names(obj$pops)))
      names = children
      if(length(names) == 2) { # we have an operator and 2 names eg pop1 & pop2 , pop1 | pop2
        ids <<- c(id, ids) 
      } else {
        # we use children if any and add last defined population name
        names = c(names, ids[1])
        # we remove last population name from the stack
        ids <<- ids[-1]
        if(ele != "!" && length(names) == 1) {
          # we have only one name for a 2 side operation (& , |)
          # we add population name from the stack to have 2 names
          names = c(names, ids[1])
          # we remove last population name from the stack
          ids <<- ids[-1]
        }
        # finally new population name defined by the operation is added to the stack
        ids <<- c(id, ids)
      }
      # already encountered names are flushed
      children <<- c()
      return(list(bool = ele, def = names, id = id))
    } else {
      # no operator found, so this is a children name
      # children name is added
      children <<- c(children, ele)
    }
    return(NULL)
  })
  op = op[sapply(op, length) != 0]
  # if length(op) is 0 it means that population is an alias of another
  # so as a hack we declare it as being a `and` of this alias
  if(length(op) == 0) op = list(list(bool="&", def=c(pop$definition,pop$definition)))
  
  # id of the last operation is the population name
  op[[length(op)]]$id <- pop$name
  
  # create nodes
  lapply(op, FUN = function(x) {
    bool = switch(x$bool, "|" = "or", "&" = "and", "!" = "not")
    xml_new_node(name = "_ns_gating_ns_BooleanGate", attrs = list("_ns_gating_ns_id" = x$id),
                 .children = list(xml_new_node(name = "_ns_data-type_ns_custom_info", attrs = c(pop[!(names(pop) %in% "split")])),
                                  xml_new_node(name = paste0("_ns_gating_ns_", bool),
                                               .children = lapply(x$def, FUN = function(def) {
                                                 xml_new_node(name = "_ns_gating_ns_gateReference", attrs = list("_ns_gating_ns_ref" = def))
                                               }))))
  })
}

#' @title GatingML Conversion from XML2
#' @description 
#' Helper to convert boolean and graphical pops and corresponding regions from XML nodes in GatingML files.
#' @param xml_nodes a set of xml_nodes
#' @param type type of the gate to extract
#' @return a list of:\cr
#' -region, region information\cr
#' -pop, population information
#' @keywords internal
fromXML2_gating  <- function(xml_nodes, type = "rect") {
  lapply(xml_nodes, function(xml_node) {
    ids = unname(xml_attrs(xml_node))
    if(length(ids) == 1) ids = c(ids, "All")
    vals = to_list_node(xml_node)
    meta = vals[["custom_info"]]
    node = vals[names(vals) != "custom_info"]
    if(type == "bool") {
      op = c("And","Or","Not")[names(node) == c("and", "or", "not")]
      nam = unname(unlist(node))
      if(op=="Not") {
        if(length(nam) != 1) stop("[not] in BooleanGate can only have one operand") # FIXME check if Not can be chained
        def = paste0(op,"|",unname(unlist(node)))
      } else {
        def = paste0(unname(unlist(node)),collapse=paste0("|",op,"|"))
      }
      pop = list(name=ids[1], type="C", base="All", definition=def) # TODO find a way to remove intermediate pop so as to use original def = meta$definition
      reg = list()
    } else {
      dimension = do.call(what = cbind, args = sapply(node[names(node) == "dimension"], simplify = FALSE, unlist))
      names_loc = grep("name", rownames(dimension), value = TRUE)
      names_dim = unname(dimension[names_loc,])
      L = length(names_dim)
      if(L < 1) stop("non compatible dimensions: less than 1")
      if(L > 2) stop("non compatible dimensions: more than 2")
      if(L == 1) {
        if(type == "rect") {
          type = "line"
        } else {
          stop("incompatible type [",type,"] and number of dimension [",L,"]")
        }
      }
      switch(type,
             "line" = {
               x = as.numeric(dimension[c("min","max"), ]) # for Range gates only min or max is provided this will generate an error in buildRegion wich expects at least 2 finite values
               if(length(meta$lineat) == 0) {              # one trick would be to define Inf at 10^100 for example ?
                 y = c(0,0)
               } else {
                 y = rep(as.numeric(meta$lineat), length.out = 2)
               }
             },
             "rect" = {
               x = as.numeric(dimension[c("min","max"), 1]) # for Range gates only min or max is provided this will generate an error in buildRegion wich expects at least 2 finite values
               y = as.numeric(dimension[c("min","max"), 2]) # for Range gates only min or max is provided this will generate an error in buildRegion wich expects at least 2 finite values
             },
             "oval" = {
               mu = unlist(node[names(node) == "mean"])
               covarianceMatrix = matrix(sapply(unlist(node[names(node) == "covarianceMatrix"]), as.numeric), 2)
               eig = eigen(covarianceMatrix)
               eig_val = sqrt(eig$values)
               if(!all(abs(diag(eig$vectors))==1)) stop("can't deal with rotated ellipse")
               x = c(-1,1) * eig_val[1] + as.numeric(mu[1])
               y = c(-1,1) * eig_val[2] + as.numeric(mu[2])
             },
             "poly" = {
               vertices = matrix(unlist(node[names(node) == "vertex"]), ncol = L)
               x = vertices[, 1]
               y = vertices[, 2]
             })
      reg = list(label=ifelse(length(meta$rtext)==0,ids[1],meta$rtext), type=type, x=x, y=y, xlogrange="P", ylogrange="P")
      if(length(meta$rcolors) != 0) {
        meta$rcolors = strsplit(meta$rcolors, split="|", fixed=TRUE)[[1]]
        reg$color = meta$rcolors[1]
        reg$lightcolor = meta$rcolors[length(meta$rcolors)]
      }
      ref = do.call(what=buildRegion, args=reg)
      if(length(meta$xtext) != 0) reg$cx = as.numeric(meta$xtext)
      if(length(meta$ytext) != 0) reg$cy = as.numeric(meta$ytext)
      if(length(meta$xlog) != 0) reg$xlogrange = meta$xlog
      if(length(meta$ylog) != 0) reg$ylogrange = meta$ylog
      if(length(meta$xtrans) != 0) reg$xtrans = meta$xtrans
      if(length(meta$ytrans) != 0) reg$ytrans = meta$ytrans
      pop = list(name=ids[1], type="G", base=ids[2], region=ids[1], fx=names_dim[1])
      if(L == 2) pop = c(pop, list(fy = names_dim[2]))
    }
    if(length(meta$pch) != 0) pop$style = meta$pch
    if(length(meta$colors) != 0) {
      meta$colors = strsplit(meta$colors, split="|", fixed=TRUE)[[1]]
      pop$color = meta$colors[1]
      pop$lightModeColor = meta$colors[length(meta$colors)]
    }
    return(list(# meta = meta, no need to return meta 
      region = reg,
      pop = do.call(what=buildPopulation, pop)))
  })
}

#' @title GatingML File Reader
#' @description
#' Extracts GatingML from file.
#' @param fileName path to file. It should be a .xml file.
#' @param ... other arguments to be passed.
#' @details reading GatingML files is in development and partly implemented.
#' For the moment, only files generated with IFC package can be read.
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
  
  # extract gating from GatingML file
  gating = suppressWarnings(xml_find_first(tmp, "/gating:Gating-ML", ns = ns))
  if(length(gating) == 0) stop(paste0(fileName, "\ndoes not seem to be well formatted: </gating:Gating-ML> not found"))
  
  # initialize returned values
  regions = list()
  pops = list(buildPopulation(name="All",color="White"))
  plots = list()
  
  # types from GatingML specification
  gates_type = c("gating:RectangleGate", "gating:EllipsoidGate", "gating:PolygonGate", "gating:BooleanGate", "gating:QuadrantGate")
  gates = sapply(gates_type, simplify = FALSE, xml_find_all, x = gating)
  foo = sapply(gates_type[1:4], simplify = FALSE, FUN = function(type) {
    switch(type, 
           "gating:RectangleGate" = {
             fromXML2_gating(xml_nodes = gates[[type]], type = "rect")
           },
           "gating:EllipsoidGate" = {
             fromXML2_gating(xml_nodes = gates[[type]], type = "oval")
           },
           "gating:PolygonGate" = {
             fromXML2_gating(xml_nodes = gates[[type]], type = "poly")
           },
           "gating:BooleanGate" = {
             fromXML2_gating(xml_nodes = gates[[type]], type = "bool")
           },
           {
             stop("[",type,"] is not supported") # TODO, maybe we can add quadrant ?
           }
    )
  })
  
  # convert pops and regions to IFC compatible
  for(i in gates_type[1:4]) {
    regions = c(regions, lapply(foo[[i]], FUN = function(k) k$region))
    pops = c(pops, lapply(foo[[i]], FUN = function(k) k$pop))
  }
  regions = regions[sapply(regions, length) != 0]
  pops = pops[sapply(pops, length) != 0]
  names(regions) = sapply(regions, FUN = function(x) x$label)
  names(pops) = sapply(pops, FUN = function(x) x$name)
  
  # access general custom info
  general_custom = suppressWarnings(xml_find_first(gating, "data-type:custom_info", ns = ns))
  
  # check if xml has been generated by IFC package
  is_ifc = suppressWarnings(xml_find_first(general_custom, "info/IFC"))
  
  # if from ifc we can extract additional tagged pop and plots
  if(length(is_ifc) != 0) {
    ##### extracts tagged pop information
    tagged=lapply(xml_attrs(xml_find_all(general_custom, "tagged//pop")), FUN=function(x) as.list(x))
    if(length(tagged)!=0) {
      tagged=lapply(tagged, FUN = function(x) {
        # modify colors
        if(length(x$colors) != 0) {
          x$colors = strsplit(x$colors, split="|", fixed=TRUE)[[1]]
          x$color = x$colors[1]
          x$lightModeColor = x$colors[length(x$colors)]
          x = x[names(x) != "colors"]
        }
        # modify names
        N = names(x)
        if("pch" %in% N) names(x)[N == "pch"] <- "style"
        # modify obj
        if("obj" %in% N) x$obj = as.integer(strsplit(x$obj, split="|", fixed=TRUE)[[1]])
        return(x)
      })
      names(tagged) = sapply(tagged, FUN=function(x) x$name)
      pops = c(pops, tagged)
      All_ = strsplit(xml_attr(is_ifc, "All"), split="|", fixed=TRUE)[[1]]
      pops[["All"]]$style = All_[1]
      pops[["All"]]$color = All_[2]
      pops[["All"]]$lightModeColor = All_[3]
    }
    
    ##### extracts graphs information
    plots=lapply(xml_attrs(xml_find_all(general_custom, "plots//plot")), FUN=function(x) as.list(x))
    if(length(plots)!=0) {
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
      rm(list=c("plots_tmp", "plot_order"))
      # now we need to modify plot with names of regions and pops
      plots = lapply(plots, FUN = function(g) {
        ##### modify pops and regions names
        N = names(g)
        N[N == "basepop"] <- "BasePop"
        N[N == "legend"] <- "Legend"
        N[N == "overlay"] <- "ShownPop"
        N[N == "region"] <- "GraphRegion"
        names(g) <- N
        if(length(g$ShownPop) != 0) g$ShownPop = lapply(splitn(definition = g$ShownPop[[1]]$displayed, all_names = names(pops)), FUN = function(x) list(name=x))
        if(length(g$GraphRegion) != 0) g$GraphRegion = lapply(splitn(definition = g$GraphRegion[[1]]$displayed, all_names = names(regions)), FUN = function(r) list(name = r))
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
