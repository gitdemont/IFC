#' @title File Extension Retrieval
#' @description 
#' Determines extension from alpha numeric file name.
#' @param x path to file.
#' @source derived from file_ext() in \pkg{tools}, R Core Team, Kurt Hornik and Friedrich Leisch.
#' @return the file extension in lower case if found. Otherwise, "".
#' @export
getFileExt <- function (x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  tolower(ifelse(pos > -1L, substring(x, pos + 1L), ""))
}
