#' @title Batch File Writer
#' @description
#' Writes an XML file to batch files
#' @param batch list of batch nodes as created by \code{\link{buildBatch}}.
#' @return It invisibly returns full path of xml batch file.
#' @export
ExportToBATCH = function(batch) {
  path = unique(unlist(lapply(batch, FUN=function(x) x$batch_dir)))
  if(length(path)!=1) stop("you ar trying to save batch file in multiple places which is not allowed")
  if(attr(x = batch[[1]]$batch_dir, which = "tempdir")) { # batch_dir is tempdir(), allowed by CRAN
    export_to = normalizePath(paste(path,"batch.xml",sep="/"), winslash = "/", mustWork = FALSE)
    write_xml(file = export_to, xml_new_node(name="batches",.children=lapply(batch, FUN=function(x) x$xml)), encoding = "utf-8")
    message(paste0("\n######################\n", export_to, "\nhas been successfully exported\n"))
  } else { # batch_dir is not tempdir(), so file will be saved to tempdir() and user will receive a message to MANUALLY rename the file to the desired location
    export_to = normalizePath(paste(tempdir(check = TRUE),"batch.xml",sep="/"), winslash = "/", mustWork = FALSE)
    write_xml(file = export_to, xml_new_node(name="batches",.children=lapply(batch, FUN=function(x) x$xml)), encoding = "utf-8")
    message(paste0("\n######################\n",
                   export_to,
                   "\nhas been successfully exported",
                   sprintf("\nTo complete the process, please run:\nfile.copy(from = '%s', to = '%s', overwrite = TRUE)\n", export_to, normalizePath(paste(path,"batch.xml",sep="/"), winslash = "/", mustWork = FALSE))))
  }
  return(invisible(export_to))
}
