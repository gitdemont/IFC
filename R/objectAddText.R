#' @title Image Text Incrustation
#' @description
#' Adds Text to image.
#' @param image a [0,1] image.
#' @param text a character string.
#' @param color a character string. color of the text.
#' @param xoff positive integer. x offset in image to start writing text.
#' @param yoff positive integer. y offset in image to start writing text.
#' @param corner a character string. where to position text in the image. Allowed are "TL", "TR", "BL", "BR", for top-left, top-right, bottom-left, bottom-right, respectively.
#' @details One-lined text will be added so has to be fully contained within image and anchored at desired corner plus x and y offset from it.
#' @return an image with text.
objectAddText <- function(image, text, color, xoff = 0, yoff = 0, corner = "TL") {
  # several checks
  text = na.omit(as.character(text))
  assert(text, len = 1, typ = "character")
  color = na.omit(as.character(color))
  assert(color, len = 1, typ = "character")
  xoff = na.omit(as.integer(xoff)); xoff = xoff[xoff>=0]
  assert(xoff, len = 1, typ = "integer")
  yoff = na.omit(as.integer(yoff)); yoff = yoff[yoff>=0]
  assert(yoff, len = 1, typ = "integer")
  corner = na.omit(as.character(corner))
  assert(corner, len = 1, alw = c("TL", "TR", "BL", "BR"))
  checkColor(color)
  
  # dim + check
  di = dim(image)
  txt_msk = texttomatrix(text)
  dt = dim(txt_msk)
  if((dt[1] + yoff > di[1]) | dt[2] + xoff > di[2]) {
    txt_msk = cpp_resize2(mat = txt_msk, new_height = di[2] - xoff, new_width = di[1] - yoff, add_noise = FALSE, bg = 0, sd= 0)
  }
  
  # modify xoff, yoff according to corner anchorage
  switch(corner,
         "TR" = {
           xoff = di[2] - dt[2] - xoff
         }, 
         "BL" = {
           yoff = di[1] - dt[1] - yoff
         }, 
         "BR" = {
           xoff = di[2] - dt[2] - xoff
           yoff = di[1] - dt[1] - yoff
         })
  
  # place text in image
  invert = FALSE
  if(is.na(di[3])) return(cpp_mark(A = image, B = txt_msk, mask = txt_msk, xoff = xoff, yoff = yoff, invert = invert))
  color = tolower(color)
  if(color=="black") invert = TRUE
  txt_img = objectColorize(txt_msk,color)
  return(array(sapply(1:di[3], FUN=function(x) cpp_mark(A = image[,,x], B = txt_img[,,x], mask = txt_msk, xoff = xoff, yoff = yoff, invert = invert)),dim = di))
}
