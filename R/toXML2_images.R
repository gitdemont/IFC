#' @title IFC_images XML Conversion
#' @description 
#' Helper to convert images (`IFC_images` object) to XML nodes.
#' @param images an `IFC_images` object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @return a xml_node.
#' @keywords internal
toXML2_images = function(images, verbose = FALSE, display_progress = TRUE, title_progress = "") {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating images node")
  assert(images, cla = "IFC_images")
  bgm = grep("^bgmean", names(images))
  bgs = grep("^bgstd", names(images))
  satc = grep("^satcount", names(images))
  satp = grep("^satpercent", names(images))
  lapply(1:nrow(images), FUN=function(i) {
    xml_new_node(name = "SO",
               attrs = c("id" = cpp_num_to_string(images[i, "id"]),
                         "imgIFD" = cpp_num_to_string(images[i, "imgIFD"]),
                         "mskIFD" = cpp_num_to_string(images[i, "mskIFD"]),
                         "spIFD" = cpp_num_to_string(images[i, "spIFD"]),
                         "w" = cpp_num_to_string(images[i, "w"]),
                         "l" = cpp_num_to_string(images[i, "l"]),
                         "fs" = cpp_num_to_string(images[i, "fs"]),
                         "cl" = cpp_num_to_string(images[i, "cl"]),
                         "ct" = cpp_num_to_string(images[i, "ct"]),
                         "objCenterX" = cpp_num_to_string(images[i, "objCenterX"]),
                         "objCenterY" = cpp_num_to_string(images[i, "objCenterY"]),
                         "bgmean" = paste0(cpp_num_to_string(unlist(images[i, bgm])), collapse = "|"),
                         "bgstd" = paste0(cpp_num_to_string(unlist(images[i, bgs])), collapse = "|"),
                         "satcount" = paste0(cpp_num_to_string(unlist(images[i, satc])), collapse = "|"),
                         "satpercent" = paste0(cpp_num_to_string(unlist(images[i, satp])), collapse = "|")))
    
  })
}
