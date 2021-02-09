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
#' Writes Gating Strategy from an `IFC_data` object to a xml file
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
    warning(paste0("some 'pops' in 'obj$pops' are tagged population and can't be exported:\n", paste0(paste0("\t- ", names(obj$pops[is_tagged]), collapse = "\n"))), call. = FALSE, immediate. = TRUE)
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
  if(any(splitp_obj$channel > 0)) message("'write_to' has %c argument but channel information can't be retrieved with writeNetwork()")
  if(any(splitp_obj$object > 0)) message("'write_to' has %o argument but channel information can't be retrieved with writeNetwork()")
  
  overwritten = FALSE
  if(file.exists(write_to)) {
    write_to = normalizePath(write_to, winslash = "/", mustWork = FALSE)
    if(!overwrite) stop(paste0("file ",write_to," already exists"))
    if(tolower(fileName) == tolower(write_to)) stop("you are trying to overwrite source file which is not allowed")
    xmlEND_export = cpp_scanFirst(fname = write_to, target = "</Network>", start = 0, end = 0)
    if(xmlEND_export > 0) {
      xml_export = read_xml(readBin(con = write_to, what = "raw", n = xmlEND_export + nchar("</Network>") - 1), options=c("HUGE","RECOVER","NOENT","NOBLANKS","NSCLEAN"))
      tryCatch({
        is_fromR = as.character(na.omit(xml_attr(xml_find_first(xml_export, "//Network"), attr = "IFC_version")))
      }, finally = rm(xml_export))
      if(length(is_fromR)==0) stop("you are trying to overwrite an original file which is not allowed")
    } else {
      stop(paste0(write_to, "\ndoes not seem to be well formatted: </Network> not found")) 
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
  
  now = format(Sys.time(), format = "%d-%b-%y %H:%M:%S")
  pkg_ver = paste0(unlist(packageVersion("IFC")), collapse = ".")
  root <- xml_new_root("Network")
  root %>% xml_set_attrs(value = c(IFC_version = pkg_ver, date = now, IDEAS_version = obj$description$Assay$IDEAS_version))
  xml_add_child(root, .value = xml_new_node(name = "SampleName", text = splitf_obj["short"]))
  # adds shared children subnodes
  sub_nodes = list(toXML2_regions(obj$regions, verbose = verbose),
                   toXML2_pops(obj$pops, verbose = verbose, display_progress = display_progress, title_progress = title_progress, ...))
  if(is_fcs) {
    xml_add_child(root, .value = xml_new_node(name = "FCS", attrs = obj$description$ID, .children = sub_nodes))
  } else {
    xml_add_child(root, .value = xml_new_node(name = "SOD", attrs = obj$description$ID, .children = sub_nodes))
  }
  xml_add_child(root, .value = toXML2_graphs(obj$graphs, verbose = verbose))
  tryCatch({
    write_xml(root, file = file_w, encoding = "utf-8")
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
