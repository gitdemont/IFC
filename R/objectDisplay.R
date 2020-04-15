#' @title Object Display
#' @description
#' This function is intended to display object extracted by \code{\link{objectExtract}}.
#' @param image An object extracted by \code{\link{objectExtract}} of class `IFC_img` or `IFC_msk`.\cr
#' Note that a matrix with finite values can also be used.
#' @param input_range a finite numeric vector of 2 values, sets the range of the input intensity values.\cr
#' Values exceeding this range are clipped. Default is 'c(0, 4095)'.
#' @param full_range if 'full_range' is TRUE, then 'input_range' will be set to 'c(0, 4095)' and 'gamma' forced to 1. Default is FALSE.
#' @param force_range if 'force_range' is TRUE, then 'input_range' will be adjusted to object range in [-4095, +inf] and 'gamma' forced to 1. Default is FALSE.\cr
#' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
#' @param gamma gamma correction. Default is 1, for no correction.
#' @param color a color. Default is "Green".
#' @details If input 'image' is of class `IFC_img` or `IFC_msk`, then if 'input_range', 'full_range', 'force_range', 'gamma' and / or 'color' parameters is/are missing,
#' it/they will be extracted from 'image' attributes.\cr
#' If input 'image' is not of one of class `IFC_img` or `IFC_msk`, then force_range will be forced to TRUE.\cr
#' An error will be thrown if input image contains non finite values.
#' @param dpi display resolution. Default is 300.
#' @export
objectDisplay = function(image, input_range = c(0, 4095), full_range = FALSE, force_range = FALSE, gamma = 1, color = "Green", dpi = 300) {
  dpi = na.omit(as.integer(dpi)); dpi = dpi[dpi>0]; dpi = dpi[is.finite(dpi)]
  assert(dpi, len = 1, typ = "integer")
  d = dim(image); 
  switch(class(image), 
         "IFC_img" = { 
           if(missing(input_range)) input_range = attr(image, "input_range")
           if(missing(full_range)) force_range = attr(image, "full_range")
           if(missing(force_range)) force_range = attr(image, "force_range")
           if(missing(gamma)) gamma = attr(image, "gamma")
           if(missing(color)) color = attr(image, "color")
           checkColor(color)
           foo = objectColorize(objectNormalize(attr(image, "raw"), input_range = input_range, full_range = full_range, force_range = force_range, gamma = gamma), color)
         },
         "IFC_msk" = {
           if(missing(input_range)) input_range = attr(image, "input_range")
           if(missing(full_range)) force_range = attr(image, "full_range")
           if(missing(force_range)) force_range = attr(image, "force_range")
           if(missing(gamma)) gamma = attr(image, "gamma")
           if(missing(color)) color = attr(image, "color")
           if(missing(color)) color = attr(image, "color")
           checkColor(color)
           foo = objectColorize(objectNormalize(attr(image, "raw"), input_range = input_range, full_range = full_range, force_range = force_range, gamma = 1), color)
         },
         "matrix" = {
           checkColor(color)
           foo = objectColorize(objectNormalize(image, force_range = TRUE), color)
         },
         {
           stop("'image' is not compatible with objectDisplay")
         })
  grid.newpage()
  do.call(what = "grid.raster", args = list(image = foo,
                                            width = unit(dpi * d[2] / 96, "points"),
                                            height = unit(dpi * d[1] / 96, "points"),
                                            interpolate = FALSE))
}
