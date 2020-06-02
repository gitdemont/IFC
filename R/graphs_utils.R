#' @title LinLog Transformation of Ticks and Labels
#' @description Helper to rescale and label axes when linlog transformation is used
#' @keywords internal
scale_trans=function(hyper, b=10, lin_comp=log(b)){
  at=outer(-10:10,-10:10,FUN=function(m,p) {m*b^p})
  toKeep=abs(at)>=hyper
  at=c(at[toKeep],seq(-hyper,hyper,length.out = 9)[2:8])
  at=smoothLinLog(x=at, hyper=hyper, base=b, lin_comp=lin_comp)
  lab=outer(-10:10,-10:10,FUN=function(m,p) {paste0(m,"*",b,"^",p)})
  lab[grep(paste0(1,"\\*",b,"\\^"),lab,invert=TRUE)]=""
  lab=c(lab[toKeep],seq(-hyper,hyper,length.out = 9)[2:8])
  lab=gsub("\\*10\\^","e",lab)
  return(list(at=at,labels=lab))
}

#' @title LinLog Transformation for IFC Graphs Plotting Scales
#' @description Helper to rescale and label axes when linlog transformation is used
#' @keywords internal
myScales=function(x=list(), y=list()) {
  if(length(x$alternating)==0) x$alternating=1
  if(length(x$tck)==0) x$tck=c(TRUE,FALSE)
  if(length(x$hyper)==0) x$hyper="P"
  if(length(x$b)==0) x$b=10
  if(length(x$lin_comp)==0) x$lin_comp=log(x$b)
  if(length(x$rot)==0) x$rot=45
  
  if(length(y$alternating)==0) y$alternating=1
  if(length(y$tck)==0) y$tck=c(TRUE,FALSE)
  if(length(y$hyper)==0) y$hyper="P"
  if(length(y$b)==0) y$b=10
  if(length(y$lin_comp)==0) y$lin_comp=log(y$b)
  if(length(y$rot)==0) y$rot=0
  
  x_scale=list("alternating"=x$alternating,"tck"=x$tck,"rot"=x$rot)
  y_scale=list("alternating"=y$alternating,"tck"=y$tck,"rot"=y$rot)
  if(x$hyper!="P") x_scale=c(x_scale, do.call("scale_trans", list(hyper=x$hyper,b=x$b,lin_comp=x$lin_comp)))
  if(y$hyper!="P") y_scale=c(y_scale, do.call("scale_trans", list(hyper=y$hyper, b=y$b, lin_comp=y$lin_comp)))
  
  if(length(x$at)!=0) x_scale=c(x_scale,"at"=x$at)
  if(length(x$lab)!=0) x_scale=c(x_scale,"labels"=x$labels)
  if(length(y$at)!=0) y_scale=c(y_scale,"at"=y$at)
  if(length(y$lab)!=0) y_scale=c(y_scale,"labels"=y$labels)
  return(list("x"=x_scale,"y"=y_scale))
}

#' @title 2D Binned Kernel Density Estimation
#' @name calcDensity
#' @description Helper to compute density plot
#' @keywords internal
calcDensity <- getFromNamespace(x = ".smoothScatterCalcDensity", ns = "grDevices")

#' @title Colors for Smooth Density Plots
#' @description Helper to map density to colors
#' @source derived from \pkg{grDevices} R Core Team, Florian Hahne at FHCRC, originally
#' @keywords internal
densCols=function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("blue","green","red")), transformation=function(x) (x^0.125)) {
  xy <- xy.coords(x, y)
  select <- is.finite(xy$x) & is.finite(xy$y)
  x <- cbind(xy$x, xy$y)[select, ]
  map <- calcDensity(x, nbin, bandwidth)
  mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
  xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
  ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
  dens <- map$fhat[cbind(xbin, ybin)]
  dens <- transformation(dens) # slightly modifyed: a transformation function has been introduced 
  dens[is.na(dens)] <- 0
  colpal <- cut(dens, length(dens), labels = FALSE)
  cols <- rep(NA_character_, length(select))
  cols[select] <- colramp(length(dens))[colpal]
  cols
}

#' @title Integer to Hexadecimal Color Conversion
#' @description Helper to convert color to hex
#' @keywords internal
colConv=function(col){
  col=strsplit(col,split="|",fixed=TRUE)[[1]]
  col=suppressWarnings(sprintf("%02X",as.numeric(col)))
  return(sapply(col, FUN=function(x) {
    foo = paste0("#",substr(x,3,8), substr(x,1,2), collapse="")
    checkColor(foo)
    foo
  }))
}

#' @title Histogram Constructor
#' @name hist_constr
#' @description Helper to construct histogram 
#' @keywords internal
hist_constr <- getFromNamespace(x = "hist.constructor", ns = "lattice")

#' @title Histogram Type Constructor
#' @description Helper to construct histogram 
#' @source derived from \pkg{lattice} from Deepayan Sarkar
#' @keywords internal
type_constr=function(x, br, type, include.lowest=TRUE, right=TRUE, plot = FALSE) {
  h=hist_constr(x=x, breaks=br, include.lowest=include.lowest, right=right, plot=plot)
  val_constr(x, h, type)
}

#' @title Histogram Val Constructor
#' @description Helper to construct histogram 
#' @source derived from \pkg{lattice} from Deepayan Sarkar
#' @keywords internal
val_constr=function(x, h, type) {
  switch(type, "count"=h$counts, "percent"=100*h$counts/length(x), "density"=h$density, "mids"=h$mid, "breaks"=h$breaks)
}

#' @title Histogram Smooth Constructor
#' @description Helper to smooth histogram
#' @source derived from \pkg{lattice} from Deepayan Sarkar
#' @keywords internal
pan_smooth = function(x, type, br, normalize, fill, lwd, lty, col, alpha, ylim, bin, border, include.lowest=TRUE, right=TRUE, factor=0) {
  h=hist_constr(x, br, include.lowest=include.lowest, right=right, plot = FALSE)
  xx=val_constr(x, h, "mids")
  yy=density(x, n=bin, na.rm=TRUE, from=min(br), to=max(br))$y
  yy[xx<min(x, na.rm=TRUE)]=0
  yy[xx>max(x, na.rm=TRUE)]=0
  # if(normalize) {yy=yy/max(yy)*max(ylim)} else {yy=yy/max(yy)*max(val_constr(x, h, type))}
  if(normalize) {yy=yy/max(yy)} else {yy=yy/max(yy)*max(val_constr(x, h, type))}
  if(fill) {
    x1=1
    x2=length(xx)
    panel.polygon(x=c(xx[c(x1,x1:x2,x2)]), y= c(0, yy[x1:x2], 0), col=col, lwd=lwd, lty=lty, alpha=alpha)
  } else {
    panel.xyplot(x=br, y=c(yy,0), col=col, type="l", lwd=lwd, lty=lty, alpha=alpha)
  }
}

#' @title Lattice Histogram Panel Contructor
#' @description Helper to create histogram panel
#' @source derived from \pkg{lattice} from Deepayan Sarkar
#' @keywords internal
pan_hist=function(x, type, br, normalize, fill, lwd, lty, col, alpha, ylim, bin, border, include.lowest=TRUE, right=TRUE){
  h=hist_constr(x, br, include.lowest=include.lowest, right=right, plot=FALSE)
  yy=val_constr(x, h, type)
  # if(normalize) { yy=yy/max(yy)*max(ylim) }
  if(normalize) { yy=yy/max(yy) }
  if(fill) {
    panel.rect(x=br[-length(br)], ybottom=0, ytop=yy, width=diff(br), col=col, lwd=lwd, lty=lty, alpha=alpha, border=border)
  } else {
    panel.xyplot(x=br, y=c(yy,0), col=col, type="s", lwd=lwd, lty=lty, alpha=alpha)
  } 
}

#' @title Histogram y-Axis Limits Constructor
#' @description Helper to extract ylim
#' @source derived from \pkg{lattice} from Deepayan Sarkar
#' @keywords internal
get_ylim=function(x, type, br, include.lowest=TRUE, right=TRUE) {
  h=hist_constr(x, br, include.lowest=include.lowest, right=right, plot=FALSE)
  yy=val_constr(x, h, type)
  return(range(0,yy,finite=TRUE,na.rm=TRUE))
}

#' @title Lattice Key Panel Contructor
#' @description Helper to add key to panel
#' @source derived from \pkg{lattice} from Deepayan Sarkar
#' @keywords internal
pan_key=function (key, corner = c(0, 0.98), x = corner[1], y = corner[2]) {
  key.gf <- draw.key(key, draw = FALSE)
  vp <- viewport(x = unit(x, "npc") + unit(0.5 - corner[1], "grobwidth", list(key.gf)), y = unit(y, "npc") + unit(0.5 - corner[2], "grobheight", list(key.gf)))
  pushViewport(vp)
  grid.draw(key.gf)
  upViewport()
}

#' @title Ellipsoid Polygon Constructor
#' @description Helper to transform gate boundaries to ellipsoid polygon.
#' @param gate list containing x and y ellipse boundaries coordinates.
#' @param theta rotation angle of the ellipse. Default is 2*pi. It should not be modified since ellipse gate are axis-aligned.
#' @param npoints number of polygon vertices desired.
#' @keywords internal
toEllipse=function(gate, theta=2*pi, npoints=100) {
  ell = cpp_ell_coord(gate$x, gate$y)
  xcoord <- seq(-ell[1],ell[1],length=npoints)
  ycoord_neg <- sqrt(ell[2]^2*(1-(xcoord)^2/ell[1]^2))
  ycoord_pos <- -sqrt(ell[2]^2*(1-(xcoord)^2/ell[1]^2))
  xx <- c(xcoord,xcoord[npoints:1])
  yy <- c(ycoord_neg,ycoord_pos)
  return(list("x"=xx*cos(2*pi-theta)+yy*sin(2*pi-theta)+ell[3],"y"=yy*cos(2*pi-theta)-xx*sin(2*pi-theta)+ell[4]))
}
