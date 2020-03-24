#' @title Character Color Check and RGB Conversion
#' @description 
#' Checks that color is R compatible by converting it to RGB matrix.
#' @param color a character string.
#' @keywords internal
checkColor = function(color) {
  color=na.omit(as.character(color))
  assert(color, len = 1, typ = "character")
  RGB = col2rgb(color)
  return(RGB)
}
