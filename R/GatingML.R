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

#' @title Spillover GatingML Conversion
#' @description 
#' Helper to convert spillover to XML nodes in GatingML files.
#' @param spillover a spillover matrix. It has to have colnames and rownames.
#' @param name the name to use to identify the node.
#' @return a xml_node.
#' @keywords internal
toXML2_spillover_gs = function(spillover, name) {
  col_n = colnames(spillover)
  row_n = rownames(spillover)
  if(length(setdiff(col_n, "")) != ncol(spillover)) stop("'spillover' should have column names")
  if(length(setdiff(row_n, "")) != nrow(spillover)) stop("'spillover' should have row names")
  xml_new_node("_ns_transforms_ns_spectrumMatrix", attrs = list("_ns_transforms_ns_id" = name,
                                                                "_ns_transforms_ns_matrix-inverted-already"="false"),
               .children = c(list(xml_new_node("_ns_transforms_ns_fluorochromes",
                                               .children = lapply(colnames(spillover), FUN = function(x) {
                                                 xml_new_node("_ns_data-type_ns_fcs-dimension", attrs = list("_ns_data-type_ns_name"=x))
                                               })),
                                  xml_new_node("_ns_transforms_ns_detectors",
                                               .children = lapply(rownames(spillover), FUN = function(x) {
                                                 xml_new_node("_ns_data-type_ns_fcs-dimension", attrs = list("_ns_data-type_ns_name"=x))
                                               }))),
                             apply(spillover, 1, FUN = function(row) {
                               xml_new_node("_ns_transforms_ns_spectrum", .children = lapply(row, FUN = function(x){
                                 xml_new_node("_ns_transforms_ns_coefficient", attrs = list("_ns_transforms_ns_value"=num_to_string(x)))
                               }))
                             }))) 
}


#' @title IFC_graphs GatingML Conversion
#' @description 
#' Helper to convert graphs (`IFC_graphs` object) to XML nodes in GatingML files.
#' @param graphs an `IFC_graphs` object.
#' @return a xml_node.
#' @keywords internal
toXML2_graphs_gs = function(graphs) {
  if(length(graphs)==0) return(list())
  graphs = lapply(graphs, FUN=function(g) {
    do.call(what = buildGraph, args = g[!grepl("order", names(g))])
    g$GraphRegion = lapply(g$GraphRegion, FUN = function(r) r[!grepl("def", names(r))]) # it is mandatory to remove def
    foo = lapply(g, FUN = function(i) {
      if(typeof(i) %in% c("integer","double")) {
        return(num_to_string(i))
      } else {
        return(i)
      }
    })
  })
  lapply(graphs, FUN=function(g) {
    gg=list(type=g$type)
    gg$size=paste0(c(g$xsize, g$ysize), collapse="|")
    gg$xran=paste0(c(g$xmin, g$xmax), collapse="|")
    gg$yran=paste0(c(g$ymin, g$ymax), collapse="|")
    gg$xlog=g$xlogrange
    gg$ylog=g$ylogrange
    gg$at=paste0(c(g$xlocation,g$ylocation,g$splitterdistance), collapse="|")
    gg$x=g$f1
    gg$y=ifelse(length(g$f2)==0, g$freq, g$f2)
    gg$scale=g$scale
    gg$bin=g$bincount
    gg$smooth=g$histogramsmoothingfactor
    # additional parameters
    gg$maxpoints=g$maxpoints
    if((length(gg$maxpoints) != 0) && ((gg$maxpoints == +Inf) || (gg$maxpoints == 1))) gg$maxpoints=NULL
    gg$xtrans=g$xtrans
    gg$ytrans=g$ytrans
    labs = list(main= g$title, x=g$xlabel, y=g$ylabel,
                title=g$graphtitlefontsize,
                regions=g$regionlabelsfontsize,
                labs=g$axislabelsfontsize,
                ticks=g$axistickmarklabelsfontsize)
    stats = list(x=g$xstats, y=g$ystats, order=g$xstatsorder, show=g$stats)
    region = paste0(unlist(g[["GraphRegion"]]), collapse = "|")
    if(region == "") region = NULL
    overlay = paste0(unlist(g[["ShownPop"]]), collapse = "|")
    if(overlay == "") overlay = NULL
    density=list()
    if(gg$type == "density") {
      tmp = grepl("density", names(g$BasePop[[1]]))
      density = g$BasePop[[1]][tmp]
      if(length(density$densitytrans) != 0 && (density$densitytrans == "return" || density$densitytrans == "")) density = density[names(density) != "densitytrans"]
      names(density) = gsub("density", "", names(density))
      names(density)[names(density) == "bincount"] <- "bin"
      names(density)[names(density) == "colorslightmode"] <- "color1"
      names(density)[names(density) == "colorsdarkmode"] <- "color2"
      density = list(xml_new_node(name = "density", attrs = density))
      g$BasePop[[1]] = g$BasePop[[1]][!tmp]
    }
    xml_new_node(name = "plot", attrs = gg[sapply(gg, length) != 0],
                 .children = c(density,
                               list(xml_new_node(name = "labels", attrs = labs)),
                               # lapply(g[["Legend"]], FUN=function(i) xml_new_node(name = "legend", attrs = i)),
                               lapply(g[["BasePop"]], FUN=function(i) { 
                                 names(i)[names(i) == "linestyle"] <- "lty"
                                 xml_new_node(name = "base", attrs = i)
                               }),
                               lapply(region, FUN=function(i) xml_new_node(name = "region", attrs = list(displayed = i))),
                               lapply(overlay, FUN=function(i) xml_new_node(name = "overlay", attrs = list(displayed = i))),
                               list(xml_new_node(name = "stats", attrs = stats)),
                               list(xml_new_node(name = "order", text = g$order))))
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
  sync = attr(reg, "sync")
  # color conversion
  reg$colors = paste0(map_color(c(reg$color,reg$lightcolor), FALSE),collapse="|")
  # only keep what is needed
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
  if(length(sync) != 0) info_attrs = c(info_attrs, list("rsync"=sync))
  
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
#' @param already already attributed names. Default is character().
#' @return a list of xml_node.
#' @keywords internal
toXML2_boolpop_gs <- function(obj, pop, already = character()) {
  # recover definition
  pop_def <- pop$split
  pop_def[pop_def == "And"]="&"
  pop_def[pop_def == "Or"] ="|"
  pop_def[pop_def == "Not"]="!"
  pop_names = !(pop_def %in% c("&", "|", "!", "(", ")"))
  pop_def[pop_names] = paste0("`",pop_def[pop_names],"`")
  pop_def = str2lang(paste0(pop_def,collapse = " "))
  pch = map_style(pop$style, toR = FALSE)
  
  # define pseudo seed for id creation
  SEED = fetch_seed(list(seed = pseudo_seed(pop$definition), "Mersenne-Twister", "Inversion", "Rounding"))
  
  # recover arithmetic tree
  expand_tree = function(x) {
    branch = lapply(x, unlist)
    if(length(branch) <= 1) return(branch)
    if(inherits(x, "("))  branch = branch[-1]
    lapply(branch, expand_tree)
  }
  
  tree = expand_tree(pop_def)
  
  # inialize values for recursive loop
  op = list()
  ids = c()
  decomp_bool <- function(x) {
    xx = unlist(x, recursive = TRUE, use.names = FALSE)
    if(length(xx) > 3) return(sapply(rev(x), decomp_bool))
    ans = sapply(xx, as.character)
    ele = ans[1]
    if(ele %in% c("!","&", "|")) { # an operator is encountered
      # we generate a random id for the population resulting of the operation
      id = gen_altnames("foo", n = 8, forbidden = c(names(op), ids, already), random_seed = SEED)
      n = ans[-1]
      if(ele == "!") { # not operation takes only one operand
        if(length(ans) == 1) { n = tail(ids, 1); ids = ids[-length(ids)] } # we dequeue last population name from the stack
      }
      if(ele %in% c("&", "|")) {
        if(length(ans) == 1) { n = rev(tail(ids, 2)); ids = ids[-(length(ids) - 0:1)] } # we dequeue last 2 population names from the stack
        if(length(ans) == 2) { n = rev(c(tail(ids, 1), n)); ids = ids[-length(ids)] } # we dequeue last population name from the stack
      }
      op <<- c(op, structure(list(list(bool = ans[1], def = n, id = id)), names = id))
      ids <<- c(ids, id)
    } else {
      ids <<- c(ids, ans) # we add population name(s) to the stack
    }
  }
  decomp_bool(tree)
  
  op = op[sapply(op, length) != 0]
  # if length(op) is 0 it means that population is an alias of another
  # so as a hack we declare it as being a `and` of this alias
  if(length(op) == 0) op = list(list(bool="&", def=c(pop$definition,pop$definition)))
  
  # id of the last operation is the population name
  op[[length(op)]]$id <- pop$name

  # create nodes and keep track of already attributed names
  list(already = names(op), 
       xml = lapply(1:length(op), FUN = function(i) {
         x=op[[i]]
         bool = switch(x$bool, "|" = "or", "&" = "and", "!" = "not")
         xml_new_node(name = "_ns_gating_ns_BooleanGate", attrs = list("_ns_gating_ns_id" = x$id),
                      .children = list(xml_new_node(name = "_ns_data-type_ns_custom_info", attrs = list(op=pop$name, pch=pch, colors=pop$colors)),
                                       xml_new_node(name = paste0("_ns_gating_ns_", bool),
                                                    .children = lapply(x$def, FUN = function(def) {
                                                      xml_new_node(name = "_ns_gating_ns_gateReference", attrs = list("_ns_gating_ns_ref" = def))
                                                    }))))
       }))
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
      pop = list(name=ids[1], type="C", base="All", definition=def, op = meta$op)
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
               vertices = matrix(unlist(node[names(node) == "vertex"]), ncol = L, byrow = TRUE)
               x = as.numeric(vertices[, 1])
               y = as.numeric(vertices[, 2])
             })
      reg_label = ifelse(length(meta$rtext)==0,ids[1],meta$rtext)
      reg = list(label=reg_label, type=type, x=x, y=y, xlogrange="P", ylogrange="P")
      if(length(meta$rcolors) != 0) {
        meta$rcolors = map_color(strsplit(meta$rcolors, split="|", fixed=TRUE)[[1]], toR = TRUE)
        reg$color = meta$rcolors[1]
        reg$lightcolor = meta$rcolors[length(meta$rcolors)]
      }
      attr(reg, "sync") <- meta$rsync
      reg = keep_attributes(reg, what=buildRegion)
      if(length(meta$xtext) != 0) reg$cx = as.numeric(meta$xtext)
      if(length(meta$ytext) != 0) reg$cy = as.numeric(meta$ytext)
      if(length(meta$xlog) != 0) reg$xlogrange = meta$xlog
      if(length(meta$ylog) != 0) reg$ylogrange = meta$ylog
      if(length(meta$xtrans) != 0) reg$xtrans = meta$xtrans
      if(length(meta$ytrans) != 0) reg$ytrans = meta$ytrans
      pop = list(name=ids[1], type="G", base=ids[2], region=reg_label, fx=names_dim[1])
      if(L == 2) pop = c(pop, list(fy = names_dim[2]))
    }
    if(length(meta$pch) != 0) pop$style = meta$pch
    if(length(meta$colors) != 0) {
      meta$colors = map_color(strsplit(meta$colors, split="|", fixed=TRUE)[[1]], toR = TRUE)
      pop$color = meta$colors[1]
      pop$lightModeColor = meta$colors[length(meta$colors)]
    }
    return(list(# meta = meta, no need to return meta 
      region = reg,
      pop = keep_attributes(pop, what=buildPopulation)))
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
  spillover = list()
  
  # extract spillover
  compensation = suppressWarnings(xml_find_all(gating, "transforms:spectrumMatrix", ns = ns))
  if(length(compensation) != 0) {
    to_invert = lapply(compensation, FUN = xml_attr, attr = "transforms:matrix-inverted-already", ns = ns)
    spillover = lapply(1:length(compensation), FUN = function(i) {
      foo = to_list_node(compensation[[i]])
      col_n = unname(unlist(foo$detectors))
      row_n = unname(unlist(foo$fluorochromes))
      vals = suppressWarnings(as.numeric(unlist(foo[names(foo) == "spectrum"])))
      vals = matrix(vals, nrow = length(row_n), ncol = length(col_n), byrow = TRUE, dimnames = list(col_n, row_n))
      if("true" %in% to_invert[[i]]) vals = solve(t(vals))
      return(vals)
    }) 
    names(spillover) = sapply(compensation, FUN = xml_attr, attr = "transforms:id", ns = ns)
  }
  
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
  op = unique(na.omit(sapply(pops, FUN = function(p) ifelse(length(p$op) == 0, NA_character_, p$op))))
  class(pops) = "IFC_pops"
  pops = popsGetAffiliation(pops)
  
  recomp_bool <- function(pops, name) {
    p = pops[[name]]
    op_names = setdiff(names(which(sapply(pops, FUN = function(p) identical(p$op, name)))), name)
    pp = pops[op_names]
    while(length(pp) != 0) {
      escape = TRUE
      for(i in seq_along(p$split)) {
        n = p$split[i] 
        if(p$split[i] %in% names(pp)) {
          escape = FALSE
          beg = head(p$split, i - 1)
          end = tail(p$split, -i)
          p$split <- c(beg, "(", pp[[n]]$split, ")", end)
          pp <- pp[setdiff(names(pp), n)]
        }
      }
      if(escape) break;
    }
    if(length(pp) != 0) return(pops)
    p$definition = paste0(p$split,collapse="|")
    pops[[name]] <- p
    pops = pops[setdiff(names(pops), op_names)]
    pops
  }
  for(i in op) pops = recomp_bool(pops, i)
  pops = sapply(pops, FUN = function(p) p[setdiff(names(p), c("split","names","op"))])
  
  # check for duplicated regions
  dup = unique(names(regions)[duplicated(names(regions))])
  sapply(dup, FUN = function(x) {
    sapply(which(names(regions) == x)[-1], FUN = function(r) {
      if(!identical(regions[[x]], regions[[r]])) stop("non-identical regions with same name [",x,"]")
    })
  })
  regions = regions[!duplicated(names(regions))]
  
  # check for duplicated population
  dup = unique(names(pops)[duplicated(names(pops))])
  sapply(dup, FUN = function(x) {
    sapply(which(names(pops) == x)[-1], FUN = function(p) {
      if(!identical(pops[[x]], pops[[p]])) stop("non-identical pops with same name [",x,"]")
    })
  })
  pops = pops[!duplicated(names(pops))]
  
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
        x$type = "T"
        x$base = "All"
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
      plots = lapply(plots, FUN = function(g) {
        plot_node = xml_find_first(general_custom, paste0("//plot[@at='",g$at,"']"))
        foo = to_list_node(plot_node)
        foo$order = xml_text(plot_node)
        g[c("size","at","xran","yran")] = sapply(g[c("size","at","xran","yran")], strsplit, split="|", fixed=TRUE)
        ans = list(type = g$type, f1=g$x,
                   xlocation=g$at[1], ylocation=g$at[2], splitterdistance=g$at[3],
                   xsize=g$size[1], ysize=g$size[2], scaletype=g$scale,
                   xmin=g$xran[1], xmax=g$xran[2], ymin=g$yran[1], ymax=g$yran[2],
                   xlogrange=g$xlog,ylogrange=g$ylog,
                   axislabelsfontsize=foo$labels$labs, axistickmarklabelsfontsize=foo$labels$ticks,
                   graphtitlefontsize=foo$labels$title, regionlabelsfontsize=foo$labels$regions,
                   xlabel=foo$labels$x, ylabel=foo$labels$y, title=foo$labels$main,
                   stats=foo$stats$show,xstats=foo$stats$x,ystats=foo$stats$y,xstatsorder=foo$stats$order,
                   # Legend=list(foo$legend),
                   BasePop=unname(lapply(foo[names(foo)=="base"],
                                         FUN = function(i) { 
                                           names(i)[names(i)=="lty"] <- "linestyle"
                                           return(i)
                                         })), 
                   order=foo$order)
        if(g$type=="histogram") {
          ans = c(ans, list(freq=g$y, histogramsmoothingfactor=g$smooth, bincount=g$bin))
        } else {
          ans = c(ans, list(f2=g$y))
          if(g$type == "density") {
            names(foo$density)[names(foo$density) == "color2"] <- "colorsdarkmode"
            names(foo$density)[names(foo$density) == "color1"] <- "colorslightmode"
            names(foo$density)[names(foo$density) == "bin"] <- "bincount"
            names(foo$density) = paste0("density", names(foo$density))
            ans$BasePop[[1]] = c(ans$BasePop[[1]], foo$density)
          }
        }
        all_names = names(regions)
        alt_names = gen_altnames(all_names)
        if(length(foo$region) != 0) ans$GraphRegion=lapply(splitn(definition = foo$region$displayed, all_names = all_names, alt_names = alt_names),
                                                           FUN = function(x) list(name=x))
        all_names = names(pops)
        alt_names = gen_altnames(all_names)
        if(length(foo$overlay) != 0) ans$ShownPop=lapply(splitn(definition = foo$overlay$displayed, all_names = all_names, alt_names = alt_names),
                                                         FUN = function(x) list(name=x))
        if(length(g$maxpoints) != 0) ans$maxpoints = g$maxpoints
        if(length(g$xtrans) != 0) ans$xtrans = g$xtrans
        if(length(g$ytrans) != 0) ans$ytrans = g$ytrans
        return(ans)
      })
      plot_order=sapply(plots, FUN=function(i_plot) as.numeric(i_plot[c("xlocation", "ylocation")]))
      plots=plots[order(unlist(plot_order[1,]),unlist(plot_order[2,]))]
    }
  }
  
  # add classes
  class(plots) <- "IFC_graphs"
  class(regions) <- "IFC_regions"
  class(pops) <- "IFC_pops"
  # returned ans
  ans = list("spillover"=spillover, "graphs"=plots, "pops"=pops, "regions"=regions)
  attr(ans, "class") <- c("IFC_gating")
  return(ans)
}
