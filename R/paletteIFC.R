#' @title R/IDEAS Color Palette Mapping
#' @description
#' Maps colors between IDEAS and R.
#' @param x either "", "palette","palette_R", to_light, to_dark. Default is "".
#' @param col a compatible color to transform to color or lightModeColor. Default is "White".\cr
#' if 'x' == to_light, function will convert 'col' to lightModeColor.\cr
#' if 'x' == to_dark, function will convert 'col' to color.\cr
#' if 'col' is not found or 'x' is anything else then a data.frame of compatible colors is returned.
#' @export
paletteIFC <- function(x = c("","palette","palette_R","to_light","to_dark")[1], col = "White") {
  assert(x, len = 1, alw = c("","palette","palette_R","to_light","to_dark"))
  Software_Colors=data.frame(matrix(c("White","LightSkyBlue","CornflowerBlue","MediumSlateBlue","Blue","Aquamarine","MediumSpringGreen","Cyan","DarkTurquoise",   
                                      "Teal","Yellow","Gold","DarkKhaki","Lime","Green","Lime","Wheat","SandyBrown","Orange",     
                                      "Tomato","Red","Pink","HotPink","Plum","Magenta","DarkOrchid","LightCoral","IndianRed",       
                                      "LightGray","Gray","Black",
                                      "Black","CornflowerBlue","CornflowerBlue","MediumSlateBlue","Blue","Aquamarine","MediumSpringGreen","Teal","DarkTurquoise",    
                                      "Teal","Gold","Gold","IndianRed","Green","Green", "Lime","DarkOrange","Tomato","DarkOrange",       
                                      "Tomato","Red","DeepPink","DeepPink","DarkOrchid","Magenta","DarkOrchid","IndianRed","IndianRed",        
                                      "Black","Gray","White"), ncol=2, byrow = FALSE), stringsAsFactors = FALSE)
  names(Software_Colors) = c("color", "lightModeColor")
  # Software_Colors = Software_Colors[!(duplicated(Software_Colors[,1]) | duplicated(Software_Colors[,2])),]
  Software_Colors$color_R = gsub("^Teal","Cyan4", Software_Colors$color)
  Software_Colors$color_R = gsub("^Green","Green4", Software_Colors$color_R)
  # Software_Colors$color_R = gsub("^Gray","Darkgray", Software_Colors$color_R)
  Software_Colors$color_R = gsub("^Lime","Chartreuse", Software_Colors$color_R)
  Software_Colors$lightModeColor_R = gsub("^Teal","Cyan4", Software_Colors$lightModeColor)
  Software_Colors$lightModeColor_R = gsub("^Green","Green4", Software_Colors$lightModeColor_R)
  # Software_Colors$lightModeColor_R = gsub("^Gray","Darkgray", Software_Colors$lightModeColor_R)
  Software_Colors$lightModeColor_R = gsub("^Lime","Chartreuse", Software_Colors$lightModeColor_R)
  if(x == "palette") return(as.character(unique(c(Software_Colors[,1], Software_Colors[,2]))))
  if(x == "palette_R") return(tolower(as.character(unique(c(Software_Colors[,3], Software_Colors[,4])))))
  columns = c(0,0)
  M = FALSE
  if(x == "to_light") { columns = c(1,3,2,4) }
  if(x == "to_dark") { columns = c(2,4,1,3) }
  if(!identical(columns,c(0,0))) M = Software_Colors[,columns[1]]==col | Software_Colors[,columns[2]]==col
  if(any(M)) return(Software_Colors[which(M), columns[3:4]])
  return(Software_Colors)
}
