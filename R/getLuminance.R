#' @title Luminance Computation
#' @description 
#' Determines color's luminance.
#' @param color a character string.
#' @keywords internal
getLuminance = function(color) {
  lum = checkColor(color)
  lum = sqrt(0.299*lum["red",1]^2 + 0.587*lum["green",1]^2 + 0.114*lum["blue",1]^2)
  return(unname(lum))
}
