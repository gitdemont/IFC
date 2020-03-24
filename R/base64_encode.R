#' @title Raw Images to Base64 Encoding
#' @description
#' Encodes raw image vector to base64 string
#' @param x a raw vector.
#' @source \url{https://en.wikibooks.org/wiki/Algorithm_Implementation/Miscellaneous/Base64}
#' @keywords internal
base64_encode <- function(x) {
  if(is.null(x)) return("")
  return(cpp_base64_encode(x))
}
