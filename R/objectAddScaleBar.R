#' @title Image Scale Bar Incrustation
#' @description
#' Adds scale bar to image
#' @param image a [0,1] image.
#' @param size positive integer. Scale's bar size in micro-meter. Default is '7'.\cr
#' This parameter can't be lesser than 6px and higher than image width + scale text.
#' @param style a character string. Scale's bar style, either 'dash' or 'line'. Default is 'dash'.
#' @param color a character string. color of the scale. Default is 'white'.
#' @param res positive integer. Resolution in ppi of 'image'. Default is '96'.
#' @param xoff positive integer. x offset in image to draw scale, starting from bottom left corner.
#' @param yoff positive integer. y offset in image to draw scale, starting from bottom left corner.
#' @return scale is added to the bottom left corner of the image.
objectAddScaleBar <- function(image, size, style=c("dash","line")[1], color="white", res=96, xoff = 0, yoff = 0) {
  # several checks
  size = na.omit(as.integer(size)); size = size[size>0]
  assert(size, len = 1, typ = "integer")
  size = as.character(size)
  style = na.omit(as.character(style))
  assert(style, len = 1, alw = c("dash","line"))
  color = na.omit(as.character(color))
  assert(color, len = 1, typ = "character")
  res = na.omit(as.integer(res)); res = res[res>=0]
  assert(res, len = 1, typ = "integer")
  xoff = na.omit(as.integer(xoff)); xoff = xoff[xoff>=0]
  assert(xoff, len = 1, typ = "integer")
  yoff = na.omit(as.integer(yoff)); yoff = yoff[yoff>=0]
  assert(yoff, len = 1, typ = "integer") 
  
  # add text
  lum = getLuminance(color)
  d = dim(image)
  bar_w = ceiling(as.numeric(size)*res/25.4)
  if(d[2] <= (bar_w+ 10 + xoff + 7*(nchar(size))+4)) stop("'scale' is outside of image width range")
  if(yoff > (d[1] - 12)) stop("'scale' is outside of image height range")
  if(style == "dash") bar_scheme = c(1,1,1,0,0,0)
  if(style == "line") bar_scheme = c(1,1,1,1,1,1)
  bar_msk = rep_len(bar_scheme,bar_w)
  bar_msk = t(sapply(1:4, FUN=function(x) bar_msk))
  bar_msk[,c(1,bar_w)] <- 1
  bar_msk[c(1,4),] <- 1
  bar_img = bar_msk
  bar_img = objectColorize(bar_img,color)
  ret = array(sapply(1:d[3], FUN=function(x) cpp_mark(A = image[,,x], B = bar_img[,,x], mask = bar_msk, xoff = 2 + xoff, yoff = d[1] - 12 + 4 - yoff, invert = ifelse(lum<128,TRUE,FALSE))),dim = d)
  ret = objectAddText(image = ret, text = size, color = color, xoff = bar_w + 4 + xoff,  yoff = d[1] - 12 - yoff)
  ret = objectAddText(image = ret, text = "|", color = color, xoff = bar_w+ 4 + xoff + 6*nchar(size),  yoff = d[1] - 10 - yoff)
  return(objectAddText(image = ret, text = "m", color = color, xoff = bar_w+ 10 + xoff + 6*nchar(size),  yoff = d[1] - 12 - yoff))
}
