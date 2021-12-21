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
#' @keywords internal
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
  
  # create root with gratingML namespaces
  root <- read_xml('
   <gating:Gating-ML xmlns:gating="http://www.isac-net.org/std/Gating-ML/v2.0/gating" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:transforms="http://www.isac-net.org/std/Gating-ML/v2.0/transformations" xmlns:data-type="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes" xsi:schemaLocation="http://www.isac-net.org/std/Gating-ML/v2.0/gating http://www.isac-net.org/std/Gating-ML/v2.0/transformations http://www.isac-net.org/std/Gating-ML/v2.0/datatypes" >
   </gating:Gating-ML>
', encoding = "utf-8")
  
  overwritten = FALSE
  if(file.exists(write_to)) {
    write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
    if(!overwrite) stop(paste0("file ",write_to," already exists"))
    if(tolower(fileName) == tolower(write_to)) stop("you are trying to overwrite source file which is not allowed")
    xmlEND_export = cpp_scanFirst(write_to, charToRaw('</gating:Gating-ML>'), start = 0, end = 0)
    if(xmlEND_export > 0) {
      xml_export = read_xml(readBin(con = write_to, what = "raw", n = xmlEND_export + nchar("</gating:Gating-ML>") - 1),
                            options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
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
  tryCatch({
    towrite = file(description = file_w, open = "wb")
  }, error = function(e) {
    stop(paste0(ifelse(overwritten,"temp ","'write_to' "), "file: ", file_w, "\ncan't be created: check name ?"))
  })
  close(towrite)
  write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
  
  # defines some variables
  now = format(Sys.time(), format = "%d-%b-%y %H:%M:%S")
  pkg_ver = paste0(unlist(packageVersion("IFC")), collapse = ".")
  is_fcs = length(obj$description$FCS)!=0
  
  # identify gatingML node
  gating = xml_find_first(root, "//gating:Gating-ML", ns = xml_ns(root))
  
  # check if there are some tagged population
  is_tagged = sapply(obj$pops, FUN = function(x) x$type == "T")
  
  # add general custom info node
  to_create = c("info","tagged","plots")[c(TRUE,any(is_tagged),TRUE)]
  gating %>% xml_add_child(xml_new_node(name = "_ns_data-type_ns_custom_info", .children = lapply(to_create, xml_new_node)))
  custom_info = xml_find_first(root, "_ns_data-type_ns_custom_info//info")
  custom_tagged = xml_find_first(root, "_ns_data-type_ns_custom_info//tagged")
  custom_plots = xml_find_first(root, "_ns_data-type_ns_custom_info//plots")
  
  # fill general custom info node with info nodes
  custom_info %>% xml_add_child(xml_new_node(name = "IFC", attrs = list(IFC_version = pkg_ver,
                                                                        date = now, 
                                                                        IDEAS_version = obj$description$Assay$IDEAS_version, 
                                                                        objcount = obj$description$ID$objcount, 
                                                                        All = paste0(c(map_style(obj$pops[["All"]]$style, toR=FALSE),
                                                                                       map_color(c(obj$pops[["All"]]$color, obj$pops[["All"]]$lightModeColor), toR = FALSE)), collapse="|"))))
  
  # fill general custom info node with tagged pop nodes
  tagged_nodes = lapply(obj$pops[is_tagged], FUN = function(pop) {
    custom_tagged %>% xml_add_child(xml_new_node(name = "pop", attrs = list(name=pop$name,
                                                                            colors=paste0(map_color(c(pop$color,pop$lightModeColor), FALSE),collapse="|"), 
                                                                            pch=map_style(pop$style, toR=FALSE), 
                                                                            obj=paste0(which(pop$obj) - 1L, collapse = "|"))))
  })
  
  # fill general custom info node with graph nodes
  graph_nodes = toXML2_graphs_gs(obj$graphs)
  lapply(graph_nodes, FUN =function(node) xml_add_child(custom_plots, node))
  
  # fill spillover if any
  gatingML_comp = list()
  if(length(obj$description$spillraw) != 0) gatingML_comp = c(gatingML_comp, list(toXML2_spillover_gs(obj$description$spillraw, name = "spillraw")))
  if(length(obj$description$spillover) != 0) gatingML_comp = c(gatingML_comp, list(toXML2_spillover_gs(obj$description$spillover, name = "spillover")))
  lapply(gatingML_comp, FUN = function(gatingML_node) xml_add_child(gating, gatingML_node))
  
  # Now, we fill GatingML with Boolean and Graphical pop + regions (C and G) in gating node
  pop_nodes = lapply(obj$pops, FUN = function(pop) {
    # color conversion
    pop$colors = paste0(map_color(c(pop$color,pop$lightModeColor), FALSE),collapse="|")
    # only keep what is need
    pop = pop[!(names(pop) %in% c("color","lightModeColor","names"))]
    
    gatingML_node = switch(pop$type, 
                           "C" = {
                             pop = pop[!(names(pop) %in% "obj")]
                             toXML2_boolpop_gs(obj, pop)
                           },
                           "G" = {
                             pop = pop[!(names(pop) %in% "obj")]
                             reg = obj$regions[[pop$region]]
                             toXML2_graphpop_gs(pop, reg)
                           })
    if(inherits(gatingML_node, "xml_node")) {
      gating %>% xml_add_child(gatingML_node) 
    } else {
      lapply(gatingML_node, FUN = function(node_ele) {
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
      x = gsub("_ns_transforms_ns_", "transforms:", x)
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
