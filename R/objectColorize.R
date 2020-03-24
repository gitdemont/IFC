#' @title IFC_object Colorizer
#' @description
#' Colorizes a [0,1] grayscale image.
#' @param mat a [0,1] numeric matrix.
#' @param color a color
#' @return a 3D array where 3rd dimension is rgb.
#' @export
objectColorize <- function(mat, color) {
  col = rgb2hsv(col2rgb(color)) # this converts named color to hsv
  return(cpp_M_HSV2RGB(mat, h = col[1], s = col[2]))
}
