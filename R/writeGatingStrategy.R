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

#' @title Gating Strategy File Writer
#' @description
#' Writes GatingML from an `IFC_data` object to a xml file
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param write_to pattern used to export file.
#' Placeholders, like "\%d/\%s_fromR.\%e", will be substituted:\cr
#' -\%d: with full path directory of 'obj$fileName'\cr
#' -\%p: with first parent directory of 'obj$fileName'\cr
#' -\%e: with extension of 'obj$fileName' (without leading .)\cr
#' -\%s: with shortname from 'obj$fileName' (i.e. basename without extension).\cr
#' Exported file extension will be deduced from this pattern. Note that it has to be a .xml.
#' @param overwrite whether to overwrite file or not. Default is FALSE.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param verbose whether to display information (use for debugging purpose). Default is FALSE.
#' @param ... other arguments to be passed.
#' @details Partial implementation of ISAC's Gating-ML 2.0 data exchange standard for gating description.
#' See Josef Spidlen et al. Cytometry A 87 683-687 (2015). \url{https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22690}\cr
#' GatingML is partly implemented because:\cr
#' -Tagged population are not part of GatingML gates\cr
#' -IDEAS/INSPIRE regions are different from the collection of gates listed in GatingML. Notably,\cr
#' --only 1 or 2 dimensions gates will be used,
#' --range gates and quadrant gates are absent from IDEAS/INSPIRE
#' --ellipse gates exist in IDEAS/INSPIRE but are axis aligned and not rotated.
#' -Transformation applied in \pkg{IFC} is not part of GatingML.
#' Nonetheless, when possible additional information are provided in dedicated custom_info field.
#' @return It invisibly returns full path of exported file.
#' @export
writeGatingStrategy = function(obj, write_to, overwrite = FALSE,
                               display_progress = TRUE, verbose = FALSE, ...) {
  dots = list(...)
  # change locale
  locale_back = Sys.getlocale("LC_ALL")
  on.exit(suppressWarnings(Sys.setlocale("LC_ALL", locale = locale_back)), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_ALL", locale = "English"))
  now = format(Sys.time(), format = "%d-%b-%y %H:%M:%S")
  
  # check mandatory param
  assert(obj, cla = "IFC_data")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  is_tagged = sapply(obj$pops, FUN = function(p) p$type == "T")
  if(any(is_tagged)) {
    warning(paste0("some 'pops' in 'obj$pops' are dependent on tagged population(s) and can't be exported:\n", paste0(paste0("\t- ", names(obj$pops[is_tagged]), collapse = "\n"))), call. = FALSE, immediate. = TRUE)
    obj = data_rm_pops(obj = obj, pops = names(obj$pops[is_tagged]), list_only = FALSE, adjust_graph = FALSE)
  }
  
  if(missing(write_to)) stop("'write_to' can't be missing")
  assert(write_to, len = 1, typ = "character")
  assert(overwrite, len = 1, alw = c(TRUE, FALSE))
  assert(display_progress, c(TRUE, FALSE))
  assert(verbose, c(TRUE, FALSE))
  
  # tests if file can be written
  fileName = normalizePath(obj$fileName, winslash = "/", mustWork = FALSE)
  title_progress = basename(fileName)
  splitf_obj = splitf(fileName)
  splitp_obj = splitp(write_to)
  write_to = formatn(splitp_obj, splitf_obj)
  file_extension = getFileExt(write_to)
  assert(file_extension, len = 1, alw = "xml")
  if(any(splitp_obj$channel > 0)) message("'write_to' has %c argument but channel information can't be retrieved with writeGatingML()")
  if(any(splitp_obj$object > 0)) message("'write_to' has %o argument but channel information can't be retrieved with writeGatingML()")
  
  overwritten = FALSE
  if(file.exists(write_to)) {
    write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
    if(!overwrite) stop(paste0("file ",write_to," already exists"))
    if(tolower(fileName) == tolower(write_to)) stop("you are trying to overwrite source file which is not allowed")
    xmlEND_export = cpp_scanFirst(fname = write_to, target = "</gating:Gating-ML>", start = 0, end = 0)
    if(xmlEND_export > 0) {
      xml_export = read_xml(readBin(con = write_to, what = "raw", n = xmlEND_export + nchar("</gating:Gating-ML>") - 1), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
      tryCatch({
        is_fromR = as.character(na.omit(xml_attr(xml_find_first(xml_export, "//IFC"), attr = "IFC_version")))
      }, finally = rm(xml_export))
      if(length(is_fromR)==0) stop("you are trying to overwrite an original file which is not allowed")
    } else {
      stop(paste0(write_to, "\ndoes not seem to be well formatted: </gating:Gating-ML> not found")) 
    }
    tmp_file = normalizePath(tempfile(), winslash = "/", mustWork = FALSE)
    overwritten = TRUE
  }
  dir_name = dirname(write_to)
  if(!dir.exists(dir_name)) if(!dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)) stop(paste0("can't create\n", dir_name))
  file_w = ifelse(overwritten, tmp_file, write_to)
  tryCatch(suppressWarnings({
    towrite = file(description = file_w, open = "wb")
  }), error = function(e) {
    stop(paste0(ifelse(overwritten,"temp ","'write_to' "), "file: ", file_w, "\ncan't be created: check name ?"))
  })
  close(towrite)
  write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
  
  # defines some variables
  now = format(Sys.time(), format = "%d-%b-%y %H:%M:%S")
  pkg_ver = paste0(unlist(packageVersion("IFC")), collapse = ".")
  is_fcs = length(obj$description$FCS)!=0
  
  # create root with gratingML namespaces
  root <- read_xml('
   <gating:Gating-ML xmlns:gating="http://www.isac-net.org/std/Gating-ML/v2.0/gating" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:transforms="http://www.isac-net.org/std/Gating-ML/v2.0/transformations" xmlns:data-type="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes" xsi:schemaLocation="http://www.isac-net.org/std/Gating-ML/v2.0/gating http://www.isac-net.org/std/Gating-ML/v2.0/transformations http://www.isac-net.org/std/Gating-ML/v2.0/datatypes" >
   </gating:Gating-ML>
', encoding = "utf-8")
  
  # identify gatingML node
  gating = xml_find_first(root, "//gating:Gating-ML", ns = xml_ns(root))
  
  # add general custom info node
  gating %>% xml_add_child(xml_new_node(name = "_ns_data-type_ns_custom_info",
                                        .children = lapply(c("info","regions","pops","plots"), xml_new_node)))
  custom_info = xml_find_first(root, "_ns_data-type_ns_custom_info//info")
  custom_pops = xml_find_first(root, "_ns_data-type_ns_custom_info//pops")
  custom_regions = xml_find_first(root, "_ns_data-type_ns_custom_info//regions")
  custom_plots = xml_find_first(root, "_ns_data-type_ns_custom_info//plots")
  
  # fill general custom info node with info nodes
  custom_info %>% xml_add_child(xml_new_node(name = "IFC", attrs = list(IFC_version = pkg_ver, date = now, IDEAS_version = obj$description$Assay$IDEAS_version, objcount = obj$description$ID$objcount)))
  
  # fill general custom info node with region custom nodes
  lapply(obj$regions, FUN = function(reg) {
    # region names conversion
    Nr = names(reg)
    Nr[Nr == "color"] <- "color1"
    Nr[Nr == "lightcolor"] <- "color2"
    names(reg) <- Nr
    
    # color conversion
    reg$color1 <- map_color(reg$color1, FALSE)
    reg$color2 <- map_color(reg$color2, FALSE)
    # collapse regions vertices
    reg$x = paste0(reg$x, collapse = "|")
    reg$y = paste0(reg$y, collapse = "|")
    custom_regions %>% xml_add_child(xml_new_node(name = "region", attrs = reg))
  })
  
  # fill general custom info node with graph custom nodes
  graph_nodes = toXML2_graphs_gs(obj$graphs)
  lapply(graph_nodes, FUN =function(node) {
    custom_plots %>% xml_add_child(node)
  })
  
  #' @title Boolean Population GatingML Conversion
  #' @description 
  #' Helper to convert boolean population to XML nodes in GatingML files.
  #' @param pop a member of `IFC_pops` object.
  #' @return a list of xml_node.
  #' @keywords internal
  toXML2_boolpop_gs <- function(pop) {
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
  
  # create pop nodes either custom or gatingML ones
  pop_nodes = lapply(obj$pops, FUN = function(pop) {
    # pop names conversion
    Np = names(pop)
    Np[Np == "color"] <- "color1"
    Np[Np == "lightModeColor"] <- "color2"
    names(pop) <- Np
    
    # color conversion
    pop$color1 <- map_color(pop$color1, FALSE)
    pop$color2 <- map_color(pop$color2, FALSE)
    
    gatingML_node = switch(pop$type, 
                           "B" = {
                             pop = pop[!(Np %in% c("obj","names"))]
                             list()
                           },
                           "C" = {
                             pop = pop[!(Np %in% c("obj","names"))]
                             toXML2_boolpop_gs(pop)
                           },
                           "G" = {
                             pop = pop[!(Np %in% c("obj","names"))]
                             reg = obj$regions[[pop$region]]
                             toXML2_gating(pop, reg)
                           },
                           "T" = {
                             pop = pop[!(Np %in% c("names"))]
                             pop$obj = paste0(which(pop$obj) - 1L, collapse = "|")
                             list()
                           })
    # whatever what pop is we export its definition in custom info
    return(list(gatingML = gatingML_node, custom = xml_new_node(name = "pop", attrs = pop[!(names(pop) %in% "split")])))
  })
  
  # fill general custom info node with population custom nodes
  lapply(pop_nodes, FUN = function(node) {
    if(length(node$custom) == 0) return(NULL)
    custom_pops %>% xml_add_child(node$custom)
  })
  
  # fill gatingML node
  foo = lapply(pop_nodes, FUN = function(node) {
    if(length(node$gatingML) == 0) return(NULL)
    if(inherits(node$gatingML, "xml_node")) {
      gating %>% xml_add_child(node$gatingML) 
    } else {
      lapply(node$gatingML, FUN = function(node_ele) {
        gating %>% xml_add_child(node_ele)
      })
    }
  })
  
  tryCatch({
    # FIXME
    # here we use an awful trick to write the file
    # because it looks impossible to create a node 
    # in XML2 with a namespace e.g. gating:Gating-ML
    # so xml file:
    # -is 1st generated with _ns_gating_ns_ 
    # -then it is read with readLines and each _ns_gating_ns_
    # is substituted with gating: within each lines
    # (same process with data-type: and transformations:)
    # -then file is written with writeLines
    # -Finally, to get good indentation, file is
    # read with read_xml and write back with write_xml
    write_xml(root, file = file_w, encoding = "utf-8")
    lines <- readLines(file_w)
    lines <- lapply(lines, FUN = function(x) {
      x = gsub("_ns_gating_ns_", "gating:", x)
      x = gsub("_ns_data-type_ns_", "data-type:", x)
      x = gsub("_ns_transformation_ns_", "transformations:", x)
    })
    doc = read_xml(paste0(lines[seq_along(lines)[-1]], collapse = ""))
    write_xml(doc, file = file_w, encoding = "utf-8")
  }, error = function(e) {
    stop(paste0("Can't create 'write_to' file.\n", write_to,
                ifelse(overwritten,"\nFile was not modified.\n","\n"),
                "See pre-file @\n", normalizePath(file_w, winslash = "/"), "\n",
                e$message), call. = FALSE)
  })
  if(overwritten) {
    mess = paste0("\n######################\n", write_to, "\nhas been successfully overwritten\n")
    if(!suppressWarnings(file.rename(to = write_to, from = file_w))) { # try file renaming which is faster
      if(!file.copy(to = write_to, from = file_w, overwrite = TRUE)) { # try file.copy if renaming is not possible
        stop(paste0("Can't copy temp file@\n", normalizePath(file_w, winslash = "/"), "\n",
                    "Can't create 'write_to' file.\n", write_to,
                    "\nFile was not modified.\n"), call. = FALSE)
      } else {
        file.remove(file_w, showWarnings = FALSE)
      }
    }
  } else {
    mess = paste0("\n######################\n", write_to, "\nhas been successfully exported\n")
  }
  message(mess)
  return(invisible(write_to))
}
