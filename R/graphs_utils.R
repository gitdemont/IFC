################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2020 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2020                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,             #
#                     Jean-Pierre Marolleau, Loïc Garçon,                      #
#                     INSERM, UPD, CHU Amiens                                  #
#                                                                              #
# DISCLAIMER:                                                                  #
# -You are using this package on your own risk!                                #
# -We do not guarantee privacy nor confidentiality.                            #
# -This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or #
# contributors be liable for any direct, indirect, incidental, special,        #
# exemplary, or consequential damages (including, but not limited to,          #
# procurement of substitute goods or services; loss of use, data, or profits;  #
# or business interruption) however caused and on any theory of liability,     #
# whether in contract, strict liability, or tort (including negligence or      #
# otherwise) arising in any way out of the use of this software, even if       #
# advised of the possibility of such damage.                                   #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with IFC. If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################


#' @title Scale Constructor for IFC Graphs Plotting
#' @description Helper to rescale and label axes
#' @keywords internal
myScales=function(x=list(), y=list()) {
  if(length(x$alternating)==0) x$alternating=1
  if(length(x$tck)==0) x$tck=c(TRUE,FALSE)
  if(length(x$hyper)==0) x$hyper="P"
  if(length(x$rot)==0) x$rot=45
  
  if(length(y$alternating)==0) y$alternating=1
  if(length(y$tck)==0) y$tck=c(TRUE,FALSE)
  if(length(y$hyper)==0) y$hyper="P"
  if(length(y$rot)==0) y$rot=0
  
  x_scale=list("alternating"=x$alternating,"tck"=x$tck,"rot"=x$rot)
  y_scale=list("alternating"=y$alternating,"tck"=y$tck,"rot"=y$rot)
  x_scale=c(x_scale, base_axis_constr(x$lim, x$hyper))
  y_scale=c(y_scale, base_axis_constr(y$lim, y$hyper))
  
  return(list("x"=x_scale,"y"=y_scale))
}

#' @title 2D Binned Kernel Density Estimation
#' @name calcDensity
#' @description Helper to compute density plot
#' @keywords internal
calcDensity=getFromNamespace(x = ".smoothScatterCalcDensity", ns = "grDevices")

#' @title Colors for Smooth Density Plots
#' @description Helper to map density to colors
#' @source derived from \pkg{grDevices} R Core Team, Florian Hahne at FHCRC, originally
#' @keywords internal
densCols=function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("blue","green","red")), transformation=function(x) (asinh(x))) {
  x_features = attr(x, "features")
  xy <- xy.coords(x, y)
  select <- is.finite(xy$x) & is.finite(xy$y)
  check_fun <- try(suppressWarnings(formals(transformation)), silent = TRUE)
  if (!inherits(check_fun, what = "try-error")) {
    x <- cbind(xy$x, xy$y)[select, ]
    map <- calcDensity(x, nbin, bandwidth)
    mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
    xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
    ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
    if ("features" %in% names(check_fun)) {
      args = list(features = x_features)
    }
    else {
      args = list(x = map$fhat[cbind(xbin, ybin)])
    }
    dens <- do.call(what = transformation, args = args)
  }
  else {
    ran <- range(x_features, finite = TRUE)
    dens <- ((x_features - ran[1])/diff(ran))
    map = NULL
  }
  dens[is.na(dens)] <- 0
  colpal <- cut(dens, length(dens), labels = FALSE)
  cols <- rep(NA_character_, length(select))
  cols[select] <- colramp(length(dens))[colpal]
  attr(cols, "matrix") <- map
  attr(cols, "density") <- dens
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

#' @title Hexadecimal to Integer Color Conversion
#' @description Helper to convert color from hex
#' @keywords internal
inv_colConv=function(col){
  x = col2rgb(col,T)
  paste0(apply(x, 2, FUN = function(i) { sum(i * c(2^16,2^8,2^1,0)) - 2^24 - i[3] }),"|",collapse="")
}


#' @title Histogram Constructor
#' @name hist_constr
#' @description Helper to construct histogram 
#' @keywords internal
hist_constr=getFromNamespace(x = "hist.constructor", ns = "lattice")

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
pan_smooth=function(x, type, br, normalize, fill, lwd, lty, col, alpha, ylim, bin, border, include.lowest=TRUE, right=TRUE, factor=0) {
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

#' @title Axis Constructor
#' @name base_axis_constr
#' @description Helper to rescale and label axes when linlog / asinh transformation is used.
#' @param lim vector of length 2 of axis extents.
#' @param trans transformation applied. Defaut is "P".
#' @param nint positive integer value indicating (approximately) the desired number of intervals. Default is 10.
#' @keywords internal
base_axis_constr=function(lim, trans = "P", nint = 10) {
  nint = na.omit(as.integer(nint)); assert(nint, len = 1, typ = "integer")
  assert(trans, len = 1)
  trans_ = parseTrans(trans)
  if(trans_$what %in% c("smoothLinLog","smoothAsinh")) {
    hyper = formals(trans_$what)$hyper
    if(length(trans_$args$hyper)!=0) hyper = trans_$args$hyper
    base = formals(trans_$what)$base
    if(length(base) == 0) base = 10
    if(length(trans_$args$base)!=0) base = trans_$args$base
    n_ticks = 0
    p_ticks = 0
    neg_log_ticks = 0
    pos_log_ticks = 0
    # compute range without expansion in original scale
    ran = diff(lim) / 1.14 * c(0.07, -0.07) + lim
    ran = applyTrans(ran, trans_, inverse = TRUE)
    n_ticks = max(ran[1], -hyper)
    p_ticks = min(ran[2], hyper)
    
    # identify pos and neg base pow within displayed range
    ran_ = ran + c(hyper, -hyper)
    ran_ = sign(ran_) * log(abs(ran_) / hyper, base)
    if(ran_[1] < 0) neg_log_ticks = floor(ran_[1])
    if(ran_[2] > 0) pos_log_ticks = ceiling(ran_[2])
    tot = pos_log_ticks - neg_log_ticks + 
      ifelse((neg_log_ticks < 0) && (n_ticks >= -hyper), 0.5, 0) + 
      ifelse((pos_log_ticks > 0) && (p_ticks <= hyper), 0.5, 0)
    # create ticks and labels
    ticks_at = c()
    ticks_lab = c()
    if(neg_log_ticks != 0) {
      at = round(sort(unique(c(outer(1:base, (max(abs(neg_log_ticks),abs(pos_log_ticks))-1):0 + log(hyper, base),
                               FUN=function(m,p) {-m*base^p}))), decreasing = FALSE),10)
      at_scaled = applyTrans(at, trans_)
      at_lab = rep("", length(at_scaled))
      neg_nint = as.integer(-neg_log_ticks / tot * nint)
      if(neg_nint > 0) {
        if(length(at) < (1 * nint)) {
          at_lab = formatC(x = at, format = "g", width = -1, digits = 1, drop0trailing = TRUE) 
        } else {
          if(base == 10) {
            pretty_lab = round(-axisTicks(log10(-range(at)), log = TRUE, nint = neg_nint),10)
          } else {
            pretty_lab = round(-base^axisTicks(log(-range(at), base), log = FALSE, nint = neg_nint),10)
          }
          # log = log10(-at[at %in% pretty_lab])
          # tmp = (log - floor(log)) == 0
          # log[!tmp] <- floor(log[!tmp])
          # mul = at[at %in% pretty_lab] / 10^log
          # at_lab[at %in% pretty_lab] <- parse(text=sprintf("%i%%.%%10^~~%i", mul, log))
          at_lab[at %in% pretty_lab] = formatC(x = at[at %in% pretty_lab] , format = "g", width = -1, digits = 2, drop0trailing = TRUE) 
        }
      }
      ticks_at = c(ticks_at, at_scaled)
      ticks_lab = c(ticks_lab, at_lab)
    }
    if(n_ticks >= -hyper) {
      at = axisTicks(c(n_ticks,0), log = FALSE, nint = ceiling(3*nint*applyTrans(abs(n_ticks), trans_)/diff(lim)/tot))
      at = at[at > n_ticks]
      if(length(at) > 0) {
        at_scaled = applyTrans(at, trans_)
        at_lab = formatC(x = at, format = "g", width = -1, digits = 4, drop0trailing = TRUE) 
        ticks_at = c(ticks_at, at_scaled)
        ticks_lab = c(ticks_lab, at_lab)
      }
    }
    if(p_ticks <= hyper) {
      at = axisTicks(c(0,p_ticks), log = FALSE, nint = ceiling(3*nint*applyTrans(abs(p_ticks), trans_)/diff(lim)/tot))
      at = at[at < p_ticks]
      if(length(at) > 0) {
        at_scaled = applyTrans(at, trans_)
        at_lab = formatC(x = at, format = "g", width = -1, digits = 4, drop0trailing = TRUE) 
        ticks_at = c(ticks_at, at_scaled)
        ticks_lab = c(ticks_lab, at_lab)
      }
    }
    if(pos_log_ticks != 0) {
      at = round(sort(unique(c(outer(1:base, 0:(max(abs(neg_log_ticks),abs(pos_log_ticks))-1) + log(hyper, base),
                               FUN=function(m,p) {m*base^p}))), decreasing = FALSE),10)
      at_scaled = applyTrans(at, trans_)
      at_lab = rep("", length(at_scaled))
      pos_nint = as.integer(pos_log_ticks / tot * nint)
      if(pos_nint > 0) {
        if(length(at) < (1 * nint)) {
          at_lab = formatC(x = at, format = "g", width = -1, digits = 1, drop0trailing = TRUE) 
        } else {
          if(base == 10) {
            pretty_lab = round(axisTicks(log10(range(at)), log = TRUE, nint = pos_nint),10)
          } else {
            pretty_lab = round(base^axisTicks(log(range(at), base), log = FALSE, nint = pos_nint),10)
          }
          at_lab[at %in% pretty_lab] = formatC(x = at[at %in% pretty_lab] , format = "g", width = -1, digits = 2, drop0trailing = TRUE) 
        }
      }
      ticks_at = c(ticks_at, at_scaled)
      ticks_lab = c(ticks_lab, at_lab)
    }
    keep = (ticks_at >= lim[1]) & (ticks_at <= lim[2])
    dup = duplicated(ticks_at)
    ticks_at = ticks_at[!dup & keep]
    ticks_lab = ticks_lab[!dup & keep]
    if(length(ticks_at) < 2) {
      at = signif(seq(from = lim[1], to = lim[2], length.out = nint), 3)
      ticks_at = applyTrans(at, trans_)
      ticks_lab = at
    }
    ticks_lab = ticks_lab[is.finite(ticks_at)]
    ticks_at = ticks_at[is.finite(ticks_at)]
    # trick to fix 
    # grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : 
    # polygon edge not found (zero-width or zero-height?)
    # Calls: <Anonymous> -> heightDetails -> heightDetails.text -> grid.Call
    ticks_lab[ticks_lab == ""] <- " "
    ord = order(ticks_at)
    return(list("at" = ticks_at[ord], "labels" = ticks_lab[ord]))
  } else {
    at = axisTicks(lim, log = FALSE, nint = nint)
    at = at[is.finite(at)]
    return(list("at" = at, "labels" = formatC(x = at, format = "g", width = -1, digits = 4, drop0trailing = TRUE)))
  }
}

#' @title Histogram Constructor for 'base' Plot
#' @name base_hist_constr
#' @description Helper to create histogram.
#' @param x a vector of values for which the histogram is desired.
#' @param type histogram type. Default is missing. Allowed are "count" and "percent".
#' @param br breakpoints given an interval and the number of pieces to break it into.
#' @param normalize whether to normalize. Default is missing.
#' @param fill whether to fill. Default is missing.
#' @param smooth whether to smooth. Default is missing.
#' @param lwd,lty,col,alpha,border graphical parameters. See par() from package 'graphics'.
#' @param include.lowest logical; if TRUE, an x[i] equal to the breaks value will be included in the first (or last, for right = FALSE) bar. This will be ignored (with a warning) unless breaks is a vector.
#' @param right logical; if TRUE, the histogram cells are right-closed (left open) intervals.
#' @keywords internal
base_hist_constr=function(x, type, br, normalize, fill, smooth, lwd, lty, col, alpha, border, include.lowest=TRUE, right=TRUE){
  assert(type, len = 1, alw = c("count", "percent"))
  normalize = as.logical(normalize); assert(normalize, len = 1, alw = c(TRUE, FALSE))
  fill = as.logical(fill); assert(fill, len = 1, alw = c(TRUE, FALSE))
  smooth = as.logical(smooth); assert(smooth, len = 1, alw = c(TRUE, FALSE))
  include.lowest = as.logical(include.lowest); assert(include.lowest, len = 1, alw = c(TRUE, FALSE))
  right = as.logical(right); assert(right, len = 1, alw = c(TRUE, FALSE))
  h = hist_constr(x, br, include.lowest=include.lowest, right=right, plot=FALSE)
  k = col2rgb(col, alpha = TRUE)
  k["alpha",] <- alpha * k["alpha",]
  b = col2rgb(border, alpha = TRUE)
  b["alpha",] <- alpha * b["alpha",]
  if(smooth) {
    xx=val_constr(x, h, "mids")
    yy=density(x, n=length(br)-1, na.rm=TRUE, from=min(br), to=max(br))$y
    yy[xx<min(x, na.rm=TRUE)]=0
    yy[xx>max(x, na.rm=TRUE)]=0
    if(normalize) {yy=yy/max(yy)} else {yy=yy/max(yy)*max(val_constr(x, h, type))}
    if(fill) {
      x1=1
      x2=length(xx)
      polygon(x=c(xx[c(x1,x1:x2,x2)]), y= c(0, yy[x1:x2], 0), col=col, lwd=lwd, lty=lty)
    } else {
      lines(x=br, y=c(yy,0), col=col, type="l", lwd=lwd, lty=lty)
    }
  } else {
    yy = val_constr(x, h, type)
    if(normalize) { yy=yy/max(yy) }
    if(fill) {
      rect(xleft = br[-length(br)], ybottom=0, ytop=yy, xright = br[-1], col=k, lwd=lwd, lty=lty, border=border)
    } else {
      lines(x=br, y=c(yy,0), col=col, type="s", lwd=lwd, lty=lty)
    }
  }
}

#' @title Device Raw Coordinates
#' @name get_coordmap_raw
#' @description Helper to extract current device plotting region
#' @source computes drawing region in a similar way as shiny:::getPrevPlotCoordmap()
#' @keywords internal
get_coordmap_raw=function() {
  usr <- graphics::par("usr")
  list(domain = list(left = usr[1], 
                     right = usr[2],
                     bottom = usr[3],
                     top = usr[4]), 
       range = list(left = graphics::grconvertX(usr[1], "user", "ndc"),
                    right = graphics::grconvertX(usr[2], "user", "ndc"),
                    bottom = graphics::grconvertY(usr[3], "user", "ndc"),
                    top = graphics::grconvertY(usr[4], "user", "ndc")))
}

#' @title Device Adjusted Coordinates
#' @name get_coordmap_adjusted
#' @description Helper to extract current device plotting region adjusted to device size
#' @param coordmap current device plotting region. Default is missing.
#' @param width current device height in pixel. Default is grDevices::dev.size("px")[1].
#' @param height current device width in pixel. Default is grDevices::dev.size("px")[2].
#' @param ratio current device ratio. Default is graphics::par('din') / graphics::par('pin').
#' @source computes drawing region in a similar way as shiny:::getPrevPlotCoordmap()
#' @keywords internal
get_coordmap_adjusted=function(coordmap,
                               width = grDevices::dev.size("px")[1],
                               height = grDevices::dev.size("px")[2],
                               ratio = graphics::par('din') / graphics::par('pin')) {
  if(missing(coordmap)) coordmap = get_coordmap_raw()
  range = list(left = coordmap$range$left * (width - 1) * ratio[1],
               right = coordmap$range$right * (width - 1) * ratio[1],
               bottom = (1 - coordmap$range$bottom) * (height - 1) * ratio[2],
               top = (1 - coordmap$range$top) * (height - 1) * ratio[2])
  return(list(domain=coordmap$domain, range=range, width = width, height = height, ratio=list(x=ratio[1],y=ratio[2])))
}

#' @title User's Coordinates to Pixels Conversion
#' @name coord_to_px
#' @description Helper map user's coordinates to pixels
#' @param coord coordinates in user system. A matrix where rows are points and with at least 2 columns named "x" and "y" for x and y coordinates, respectively.
#' @param coordmap current device adjusted coordinates. Default is missing.
#' @param clipedge whether to clip points outside of plotting region to the edge. Default is FALSE.
#' @return a 2-columns matrix with "x" and "y" coordinates.
#' @keywords internal
coord_to_px=function (coord, coordmap, clipedge = FALSE) {
  if(missing(coordmap)) coordmap = get_coordmap_adjusted()
  ran_x = range(coordmap$domain$left, coordmap$domain$right)
  dx = coordmap$domain$right - coordmap$domain$left
  ran_y = range(coordmap$domain$bottom, coordmap$domain$top)
  dy = coordmap$domain$top - coordmap$domain$bottom
  ran_img_width = range(coordmap$range$right, coordmap$range$left)
  width = diff(ran_img_width)
  ran_img_height = range(coordmap$range$bottom, coordmap$range$top)
  height = diff(ran_img_height)
  cpp_coord_to_px(x = coord$x, y = coord$y, 
                  param = c(xmin = ran_x[1], xmax = ran_x[2],
                            ymin = ran_y[1], ymax = ran_y[2], 
                            width / dx, height / dy,
                            domain_l = coordmap$domain$left,
                            domain_b = coordmap$domain$bottom,
                            img_w_1 = ran_img_width[1],
                            img_h_2 = ran_img_height[2],
                            ratio_x = coordmap$ratio$x,
                            ratio_y = coordmap$ratio$y,
                            edge = clipedge))
  # return(cbind(x = ((coord$x - coordmap$domain$left)/dx * width + ran_img_width[1])     /coordmap$ratio$x,
  #              y = (ran_img_height[2] - (coord$y - coordmap$domain$bottom)/dy * height) /coordmap$ratio$y, 
  #              draw = (coord$x <= ran_x[2]) & (coord$x >= ran_x[1]) & 
  #                (coord$y <= ran_y[2]) & (coord$y >= ran_y[1])))
}

#' @title `IFC_plot` Conversion to 'base' Plot
#' @name plot_base
#' @description Helper to convert `IFC_plot` to 'base' plot.
#' @param obj an object of class `IFC_plot` as created by \code{\link{plotGraph}}.
#' @keywords internal
plot_base=function(obj) {
  old_mar = par("mar")
  on.exit(par("mar" = old_mar), add = TRUE)
  # variables for future use
  n_ticks = 10
  # check obj is `IFC_plot`
  assert(obj, cla = "IFC_plot")
  
  # short names
  D = obj$input$data
  displayed = obj$input$displayed
  basepop = obj$input$base
  Xlim = obj$input$xlim
  Ylim = obj$input$ylim
  disp_n = names(displayed)
  main = obj$input$title
  subtitle = FALSE
  if(any(obj$input$add_key %in% c("global", "both"))) {
    par("mar" = c(old_mar[1:2],old_mar[3]+length(displayed) * obj$plot$par.settings$add.text$cex * 0.5 - 1,old_mar[4]))
    main = " "
  }
  
  # common plot args
  args_plot = list(xlim = obj$input$xlim, ylim = obj$input$ylim, 
                   main = trunc_string(main, obj$input$trunc_labels), 
                   xlab = trunc_string(obj$input$xlab, obj$input$trunc_labels),
                   ylab = trunc_string(obj$input$ylab, obj$input$trunc_labels),
                   cex.lab = obj$plot$par.settings$par.xlab.text$cex,
                   cex.main = obj$plot$par.settings$par.main.text$cex,
                   cex.axis = obj$plot$par.settings$axis.text$cex,
                   axes = FALSE)
  if(args_plot$main == "") args_plot$main = " "
  if(args_plot$xlab == "") args_plot$xlab = " "
  if(args_plot$ylab == "") args_plot$ylab = " "
  
  # draw plot
  if(obj$input$type %in% c("percent", "count")) {
    # 1D
    br = do.breaks(Xlim, obj$input$bin)
    do.call(args = c(list(x = D[, "x2"],
                          border = "transparent",
                          freq = obj$input$type == "count",
                          breaks = br,
                          col = "transparent"),
                     args_plot),
            what = hist)
    if(length(displayed) > 0) {
      for(disp in disp_n) {
        if(any(D[,disp]))
          base_hist_constr(x = D[D[,disp], "x2"],
                           br = br,
                           type = obj$input$type, 
                           normalize = obj$input$normalize, 
                           smooth = obj$input$histogramsmoothingfactor,
                           fill = basepop[[obj$input$order[disp]]]$fill=="true",
                           alpha = 0.8, lwd=1,
                           col = displayed[[disp]][c("color","lightModeColor")][[obj$input$mode]],
                           border = displayed[[disp]][c("color","lightModeColor")][[obj$input$mode]],
                           lty = c(1,2,3,4,6)[match(basepop[[obj$input$order[disp]]]$linestyle,c("Solid","Dash","Dot","DashDot","DashDotDot"))])
        
      }
    }
  } else {
    # 2D
    hasdata=FALSE
    pch=16
    if((nrow(obj$input$data) > 0) && any(obj$input$subset) && (nrow(obj$input$data[obj$input$subset,]) > 0)) hasdata <- TRUE
    if(obj$input$type == "density") {
      col=c("black","white")[obj$input$mode]
      colramp=colorRampPalette(colConv(basepop[[1]][c("densitycolorsdarkmode","densitycolorslightmode")][[obj$input$mode]]))
      args_level = basepop[[1]][["densitylevel"]]
      if((length(args_level) != 0) && (args_level != "") && hasdata) {
        col = densCols(x = structure(obj$input$data$x2[obj$input$subset], features=attr(obj$input$data,"features")),
                       y = obj$input$data$y2[obj$input$subset],
                       colramp=colramp,
                       nbin=obj$input$bin,
                       transformation=function(x) {x})
        args_level = strsplit(args_level,split="|",fixed=TRUE)[[1]]
        fill = args_level[1] == "true"
        dolines = args_level[2] == "true"
        nlevels = as.integer(args_level[3])
        lowest = as.numeric(args_level[4])
        z = attr(col, "matrix")
        dd = attr(col, "density")
        zz = as.vector(z$fhat)
        at_probs = seq(from=lowest, to=1, length.out=1+nlevels+(lowest!=0))
        at = quantile(dd, at_probs, na.rm = TRUE, names = FALSE)
        do.call(args = c(list(x = obj$input$data$x2[obj$input$subset][1],
                              y = obj$input$data$y2[obj$input$subset][1],
                              pch = pch,
                              col = c("black","white")[obj$input$mode]),
                         args_plot),
                what = plot)
        low=0
        if(lowest != 0) {
          low = at[1]; at = at[-1]
          points(x=obj$input$data$x2[obj$input$subset][dd<low], y=obj$input$data$y2[obj$input$subset][dd<low], pch=".", col=c("white","black")[obj$input$mode]) #col[dd<low])
        }
        if(fill) {
          .filled.contour(x=z$x1, y=z$x2, z$fhat, levels=c(low, at), col=c(NA,colramp(length(at)-1)))
        } 
        if(dolines) {
          if(fill) {
            contour_cols = c("white","black")[obj$input$mode]
          } else {
            contour_cols = colramp(length(at))
          }
          contour(x=z$x1, y=z$x2, z$fhat, levels=at, col=contour_cols, drawlabels = FALSE, add = TRUE)
        }
      } else {
        if(hasdata) col = densCols(x = structure(obj$input$data$x2[obj$input$subset], features=attr(obj$input$data,"features")),
                                   y = obj$input$data$y2[obj$input$subset],
                                   colramp = colramp,
                                   nbin = obj$input$bin,
                                   transformation = obj$input$trans)
        do.call(args = c(list(x = obj$input$data$x2[obj$input$subset],
                              y = obj$input$data$y2[obj$input$subset],
                              pch = pch,
                              col = col),
                         args_plot),
                what = plot)
      }
      if((length(args_level) == 0) || (args_level == ""))
        if(!(inherits(obj$input$trans, what="function") || 
           !inherits(try(suppressWarnings(formals(obj$input$trans)), silent = TRUE), what="try-error"))) subtitle = TRUE
    } else {
      if(obj$input$precision == "full") {
        disp = disp_n[length(displayed)]
        do.call(args = c(list(x = obj$input$data[obj$input$data[,disp], "x2"][obj$input$subset],
                              y = obj$input$data[obj$input$data[,disp], "y2"][obj$input$subset],
                              pch = displayed[[disp]]$style, 
                              col = displayed[[disp]][c("color","lightModeColor")][[obj$input$mode]]),
                         args_plot),
                what = plot)
        if(length(displayed) > 1) {
          for(disp in rev(disp_n)[-1]) {
            points(x = obj$input$data[obj$input$data[,disp], "x2"][obj$input$subset],
                   y = obj$input$data[obj$input$data[,disp], "y2"][obj$input$subset],
                   pch = displayed[[disp]]$style, 
                   col = displayed[[disp]][c("color","lightModeColor")][[obj$input$mode]])
          }
        }
      } else {
        if(hasdata) {
          if(length(disp_n) > 1) {
            groups = apply(as.data.frame(D[obj$input$subset,disp_n]), 1, FUN=function(x) {
              tmp = which(x)[1]
              if(is.na(tmp)) return(NA)
              return(disp_n[tmp])
            })
          } else {
            groups = disp_n
          }
          pch = sapply(groups, FUN = function(disp) displayed[[disp]]$style)
          col = sapply(groups, FUN = function(disp) displayed[[disp]][c("color","lightModeColor")][[obj$input$mode]])
        } else {
          col = NA 
        }
        do.call(args = c(list(x = obj$input$data$x2[obj$input$subset],
                              y = obj$input$data$y2[obj$input$subset],
                              pch = pch,
                              col = col),
                         args_plot),
                what = plot)
      }
    }
  }
  
  # axis
  x_ticks = base_axis_constr(lim = Xlim, trans = obj$input$trans_x, nint = n_ticks)
  y_ticks = base_axis_constr(lim = Ylim, trans = obj$input$trans_y, nint = n_ticks)
  x_axis = axis(side = 1, at = x_ticks$at, labels = FALSE)
  text(x = x_axis, y = Ylim[1] - diff(Ylim) * 0.07, labels = x_ticks$labels, xpd=TRUE, adj = c(1, 1),
       cex = obj$plot$par.settings$axis.text$cex, cex.axis = obj$plot$par.settings$axis.text$cex, srt=45)
  y_axis = axis(side = 2, at = y_ticks$at, labels = FALSE)
  text(y = y_axis, x = Xlim[1] - diff(Xlim) * 0.07, labels = y_ticks$labels, xpd=TRUE, adj = c(1, 0.5),
       cex = obj$plot$par.settings$axis.text$cex, cex.axis = obj$plot$par.settings$axis.text$cex, las=1)
  box()
  
  # regions
  for(reg in obj$input$regions) {
    k = reg[c("color","lightcolor")][[obj$input$mode]]
    coords = reg[c("x","y")]
    trans_x = parseTrans(obj$input$trans_x)
    coords$x = applyTrans(coords$x, trans_x)
    reg$cx = applyTrans(reg$cx, trans_x)
    lab =  trunc_string(reg$label, obj$input$trunc_labels)
    if(reg$type=="line") {
      if(reg$cy == 0) reg$cy = diff(Ylim)*0.6 # allow to show label when it is on the axe
      if(coords$y[1] == 0) coords$y = rep(diff(Ylim)*.5, length.out=2) # allow to show line when on the axe
      text(x=reg$cx, y=reg$cy*diff(Ylim), col=k, labels=lab, pos=4, cex=obj$plot$par.settings$add.text$cex)
      polygon(x=coords$x, y=coords$y*diff(Ylim), col = k, border = k)
    } else {
      trans_y = parseTrans(obj$input$trans_y)
      coords$y = applyTrans(coords$y, trans_y)
      reg$cy = applyTrans(reg$cy, trans_y)
      if(reg$type=="rect") {
        coords$x=c(coords$x[1],coords$x[1],coords$x[2],coords$x[2])
        coords$y=c(coords$y[1],coords$y[2],coords$y[2],coords$y[1])
      }
      if(reg$type=="oval") {
        coords = toEllipse(coords)
      }
      text(x=reg$cx, y=reg$cy, col=k, labels=lab, pos=4, cex=obj$plot$par.settings$add.text$cex) 
      polygon(x=coords$x, y=coords$y, border=k, col="transparent", lwd=1, lty=1)
    }
  }
  
  # key / subtitle / title
  sub_lab = obj$input$trans
  if(sub_lab == "") sub_lab = " "
  args_sub = list(text = sub_lab, side = 3, line = 0.2, adj = 0.5, font = 3, cex = obj$plot$par.settings$par.main.text$cex * 0.8)
  if(any(obj$input$add_key %in% c("panel","global","both"))) {
    args_key = list(x="topleft",inset=0.025,
                    col=sapply(displayed, FUN=function(p) p[c("color","lightModeColor")][[obj$input$mode]]),
                    legend=disp_n,cex=obj$plot$par.settings$add.text$cex * 0.5,bg="#ADADAD99",pt.cex=1,bty="o",box.lty=0)
    if(obj$input$add_key %in% c("panel","both")) do.call(args=c(list(pch=sapply(displayed, FUN=function(p) p$style)), args_key), what=legend)
    if(obj$input$add_key %in% c("global","both")) {
      args_key$x = "top"
      args_key$inset = -graphics::grconvertY(1+ length(displayed), "chars", "nfc")
      args_key$xpd = TRUE
      do.call(args=c(list(pch=sapply(displayed, FUN=function(p) p$style)), args_key), what=legend)
      mtext(side = 3, line = 2 + length(displayed) * obj$plot$par.settings$add.text$cex * 0.5, adj = 0.5, args_plot$main, font = 2)
      args_sub$line = 1.1 + length(displayed) * obj$plot$par.settings$par.main.text$cex * 0.8
    } 
  }
  # subtitle
  if(subtitle) do.call(args = args_sub, what = mtext)
}

#' @title `IFC_plot` Conversion to 'raster' Plot
#' @name plot_raster
#' @description Helper to convert `IFC_plot` to 'raster' plot.
#' @param obj an object of class `IFC_plot` as created by \code{\link{plotGraph}}.
#' @keywords internal
plot_raster=function(obj) {
  if(obj$input$type %in% c("count", "percent")) return(plot_base(obj))
  if(obj$input$type == "density") {
    basepop = obj$input$base
    args_level = basepop[[1]][["densitylevel"]]
    if((length(args_level) != 0) && (args_level != "")) return(plot_base(obj))
  }
  # determines population order
  basepop = obj$input$base
  displayed = obj$input$displayed
  disp_n = names(displayed)
  
  old_mar = par("mar")
  on.exit(par("mar" = old_mar), add = TRUE)
  
  # copy obj and empty data
  subtitle = FALSE
  graph = obj
  graph$input$data <- graph$input$data[rep(FALSE, nrow(graph$input$data)),,drop = FALSE]
  graph$input$precision <- "full"
  graph$input$regions <- list()
  graph$input$add_key <- FALSE
  if(any(obj$input$add_key %in% c("global", "both"))) {
    par("mar" = c(old_mar[1:2],old_mar[3]+length(displayed) * obj$plot$par.settings$add.text$cex * 0.5 - 1,old_mar[4]))
    graph$input$title = ""
    graph$input$trans = ""
  }
  if(!(inherits(obj$input$trans, what="function") || 
     !inherits(try(suppressWarnings(formals(obj$input$trans)), silent = TRUE), what="try-error"))) subtitle = TRUE

  # create empty plot
  plot_base(graph)
  
  # determines current device plotting region
  coordmap = get_coordmap_adjusted()
  
  # create data specific list for raster plot
  if((obj$input$precision == "light") && (length(disp_n) > 1)) {
    set = disp_n[apply(obj$input$data[,disp_n, drop = FALSE], 1, FUN = function(x) {
      foo = which(x)[1]
    })]
  }
  data = sapply(rev(disp_n), simplify = FALSE, USE.NAMES = TRUE, FUN = function(p) {
    if((obj$input$precision == "light") && (length(disp_n) > 1)) {
      sub_ = obj$input$data[, p] & obj$input$subset & (set == p)
    } else {
      sub_ = obj$input$data[, p] & obj$input$subset
    } 
    if(sum(sub_) == 0) return(NULL)
    coords = obj$input$data[, c("x2","y2")]
    colnames(coords) = c("x","y")
    if(obj$input$type == "scatter") {
      size = 4
      col = map_color(obj$input$displayed[[p]]$lightModeColor)
    } else {
      size = 6
      colramp = colorRampPalette(colConv(obj$input$base[[1]][c("densitycolorsdarkmode", "densitycolorslightmode")][[obj$input$mode]]))
      if((sum(sub_) < 20000) || inherits(try(suppressWarnings(formals(obj$input$trans)), silent = TRUE), "try-error")) {
        col = densCols(x = structure(coords[sub_,"x"], features=attr(obj$input$data,"features")),
                       y = coords[sub_,"y"],
                       colramp = colramp,
                       nbin = obj$input$bin,
                       transformation = obj$input$trans)
      } else {
        col = colramp(255)
      }
    }
    list(size = size,
         pch = obj$input$displayed[[p]]$style,
         col = rbind(col2rgb(col), 255),
         coords = coord_to_px(coord=coords[sub_,,drop=FALSE], coordmap=coordmap),
         blur_size = 9,
         blur_sd = 3)
  })
  data = data[sapply(data, length) != 0]
  if(length(data) != 0) {
    zoom = 1
    # image is zoom fold more the size of the device and then raster size is zoom fold less, this should allow anti-aliasing
    # however this could lead to much more longer computation times
    # call c part to produce image raster
    img = cpp_raster(width = zoom * grDevices::dev.size("px")[1], height = zoom * grDevices::dev.size("px")[2], data)
    usr = par("usr")
    # subset img to drawing region
    lims = round(c(graphics::grconvertX(usr[1], "user", "ndc") * (grDevices::dev.size("px")[1]),
                   graphics::grconvertX(usr[2], "user", "ndc") * (grDevices::dev.size("px")[1]),
                   (1 - graphics::grconvertY(usr[3], "user", "ndc")) * (grDevices::dev.size("px")[2]),
                   (1 - graphics::grconvertY(usr[4], "user", "ndc")) * (grDevices::dev.size("px")[2])))
    bg = img[lims[3]:lims[4], lims[1]:lims[2],]
    # add image to plot
    rasterImage(bg / 255, xleft = usr[1], xright = usr[2], ybottom = usr[4], ytop = usr[3], interpolate = FALSE)
    # rasterImage is faster than grid.raster and allows to fit bg when it is resized
  }
  # redraw regions
  for(reg in obj$input$regions) {
    k = reg[c("color","lightcolor")][[obj$input$mode]]
    coords = reg[c("x","y")]
    trans_x = parseTrans(obj$input$trans_x)
    coords$x = applyTrans(coords$x, trans_x)
    reg$cx = applyTrans(reg$cx, trans_x)
    lab =  trunc_string(reg$label, obj$input$trunc_labels)
    if(reg$type=="line") {
      Ylim = obj$input$ylim
      if(reg$cy == 0) reg$cy = diff(Ylim)*0.6 # allow to show label when it is on the axe
      if(coords$y[1] == 0) coords$y = rep(diff(Ylim)*.5, length.out=2) # allow to show line when on the axe
      text(x=reg$cx, y=reg$cy*diff(Ylim), col=k, labels=lab, pos=4, cex=obj$plot$par.settings$add.text$cex)
      polygon(x=coords$x, y=coords$y*diff(Ylim), col = k, border = k)
    } else {
      trans_y = parseTrans(obj$input$trans_y)
      coords$y = applyTrans(coords$y, trans_y)
      reg$cy = applyTrans(reg$cy, trans_y)
      if(reg$type=="rect") {
        coords$x=c(coords$x[1],coords$x[1],coords$x[2],coords$x[2])
        coords$y=c(coords$y[1],coords$y[2],coords$y[2],coords$y[1])
      }
      if(reg$type=="oval") {
        coords = toEllipse(coords)
      }
      text(x=reg$cx, y=reg$cy, col=k, labels=lab, pos=4, cex=obj$plot$par.settings$add.text$cex) 
      polygon(x=coords$x, y=coords$y, border=k, col="transparent", lwd=1, lty=1)
    }
  }
  
  # redraw key / subtitle / title
  main = trunc_string(obj$input$title, obj$input$trunc_labels)
  if(main == "") main = " "
  sub_lab = obj$input$trans
  if(sub_lab == "") sub_lab = " "
  args_sub = list(text = sub_lab, side = 3, line = 0.2, adj = 0.5, font = 3, cex = obj$plot$par.settings$par.main.text$cex * 0.8)
  if(any(obj$input$add_key %in% c("panel","global","both"))) {
    args_key = list(x="topleft",inset=0.025,
                    col=sapply(displayed, FUN=function(p) p[c("color","lightModeColor")][[obj$input$mode]]),
                    legend=disp_n,cex=obj$plot$par.settings$add.text$cex * 0.5,bg="#ADADAD99",pt.cex=1,bty="o",box.lty=0)
    if(obj$input$add_key %in% c("panel","both")) do.call(args=c(list(pch=sapply(displayed, FUN=function(p) p$style)), args_key), what=legend)
    if(obj$input$add_key %in% c("global","both")) {
      args_key$x = "top"
      args_key$inset = -graphics::grconvertY(1+ length(displayed), "chars", "nfc")
      args_key$xpd = TRUE
      do.call(args=c(list(pch=sapply(displayed, FUN=function(p) p$style)), args_key), what=legend)
      mtext(side = 3, line = 2 + length(displayed) * obj$plot$par.settings$add.text$cex * 0.5, adj = 0.5, text = main, font = 2)
      args_sub$line = 1.1 + length(displayed) * obj$plot$par.settings$par.main.text$cex * 0.8
    } 
  }
  # subtitle
  if(subtitle) do.call(args = args_sub, what = mtext)
  box()
  return(invisible(img))
}

#' @title `IFC_plot` Conversion to 'lattice' Plot
#' @name plot_lattice
#' @description Helper to convert `IFC_plot` to 'lattice' plot.
#' @param obj an object of class `IFC_plot` as created by \code{\link{plotGraph}}.
#' @keywords internal
plot_lattice=function(obj) {
  # check obj is `IFC_plot`
  assert(obj, cla = "IFC_plot")
  
  # short names
  xy_subset = obj$input$subset
  P = obj$input$displayed
  R = obj$input$regions
  D = obj$input$data[xy_subset,,drop=FALSE]
  #browser()
  nbin = obj$input$bin
  basepop = obj$input$base
  Xlim = obj$input$xlim
  Ylim = obj$input$ylim
  Xtrans = obj$input$trans_x
  Ytrans = obj$input$trans_y
  trans_x = parseTrans(Xtrans)
  trans_y = parseTrans(Ytrans)
  reg_n = names(R)
  displayed = P
  displayed_n = names(P)
  displayed_o = obj$input$order
  type = obj$input$type
  normalize = obj$input$normalize
  color_mode = obj$input$mode
  trans = obj$input$trans
  precision = obj$input$precision
  add_key = obj$input$add_key
  lt = obj$plot$par.settings
  histogramsmoothingfactor = obj$input$histogramsmoothingfactor
  
  trunc_labels = obj$input$trunc_labels
  main = trunc_string(obj$input$title, trunc_labels) 
  if(main == "") main = " "
  xlab = trunc_string(obj$input$xlab, trunc_labels)
  if(xlab == "") xlab = " "
  ylab = trunc_string(obj$input$ylab, trunc_labels)
  if(ylab == "") ylab = " "
  if(type %in% c("percent", "count")) {
    # define legend
    KEY = list("text"=list(displayed_n),
               "cex"=lt$add.text$cex * 0.5,
               "lines"=list(col = sapply(displayed_n, FUN=function(p) P[[p]][c("color","lightModeColor")][[color_mode]]),
                            lty = sapply(displayed_n, FUN=function(r) c(1,2,3,4,6)[match(basepop[[displayed_o[r]]]$linestyle,c("Solid","Dash","Dot","DashDot","DashDotDot"))])))
    # make histogram
    if(nrow(D) > 0) {
      br = do.breaks(Xlim, nbin)
      foo = histogram(~ D[,"x2"], auto.key=FALSE,
                      xlim = Xlim, ylim = Ylim, main = main, xlab = xlab,
                      ylab = obj$input$ylab,
                      scales =  myScales(x=list(lim = Xlim, "hyper"=Xtrans), y=list(lim = Ylim, "hyper"=Ytrans)),
                      border = "transparent", nint = nbin, type = type, breaks = br, normalize = normalize,
                      panel = function(x, ...) { })
      for(l in length(displayed_n):1) {
        disp = displayed_n[l]
        if(any(D[,disp])) { # adds layer only if there is at least one point
          tmp = histogram(~ D[,"x2"], auto.key=FALSE, subset = D[,disp], alpha = 0.8,
                          col = P[[disp]][c("color","lightModeColor")][[color_mode]], border="transparent",
                          fill = as.logical(basepop[[displayed_o[disp]]]$fill=="true"),
                          lty = c(1,2,3,4,6)[match(basepop[[displayed_o[disp]]]$linestyle,c("Solid","Dash","Dot","DashDot","DashDotDot"))],
                          nint = nbin, type = type, breaks = br, normalize = normalize, Ylim = Ylim,
                          panel = function(x, type, breaks, normalize, fill, nint, border, col, alpha, lty, Ylim = Ylim, ...) {
                            if(histogramsmoothingfactor > 0) {
                              pan_smooth(x=x, type=type, br=breaks, normalize=normalize, fill=fill, lwd=1, lty=lty, col=col, alpha=alpha, ylim=Ylim, bin=nint, border=border, factor=histogramsmoothingfactor)
                            } else {
                              pan_hist(x=x, type=type, br=breaks, normalize=normalize, fill=fill, lwd=1, lty=1, col=col, alpha=alpha, ylim=Ylim, bin=nint, border=border)
                            }
                            if(l == 1) {
                              if(any(c("panel","both")%in%add_key)) pan_key(key=c(KEY,"background"="lightgrey","alpha.background"=0.8), x = 0.02)
                              lapply(reg_n, FUN=function(r) {
                                reg = R[[r]] 
                                col = reg[c("color","lightcolor")][[color_mode]]
                                coords = reg[c("x","y")]
                                coords$x = applyTrans(coords$x, trans_x)
                                reg$cx = applyTrans(reg$cx, trans_x)
                                lab =  trunc_string(reg$label, trunc_labels)
                                if(reg$cy == 0) reg$cy = diff(Ylim)*0.6 # allow to show label when it is on the axe
                                if(coords$y[1] == 0) coords$y = rep(diff(Ylim)*.5, length.out=2) # allow to show line when on the axe
                                panel.text(x=reg$cx, y=reg$cy*diff(Ylim), col=col, labels=lab, pos=4)
                                panel.lines(x=coords$x, y=coords$y*diff(Ylim),col=col)
                              })
                            }
                          })
          foo = foo + as.layer(tmp, opposite = FALSE, axes = NULL)
        }
      }
    } else {
      foo = histogram(0 ~ 0, auto.key=FALSE,
                      xlim = Xlim, ylim = Ylim, main = main, xlab = xlab,
                      ylab = obj$input$ylab,
                      scales =  myScales(x=list(lim = Xlim, "hyper"=Xtrans), y=list(lim = Ylim, "hyper"=Ytrans)), border = "transparent",
                      nint = nbin, type = type, normalize = normalize, Ylim = Ylim,
                      panel = function(x, Ylim = Ylim, ...) {
                        if(any(c("panel","both")%in%add_key)) pan_key(key=c(KEY,"background"="lightgrey","alpha.background"=0.8), x = 0.02)
                        lapply(reg_n, FUN=function(r) {
                          reg = R[[r]] 
                          col = reg[c("color","lightcolor")][[color_mode]]
                          coords = reg[c("x","y")]
                          coords$x = applyTrans(coords$x, trans_x)
                          reg$cx = applyTrans(reg$cx, trans_x)
                          lab = trunc_string(reg$label, trunc_labels)
                          if(reg$cy == 0) reg$cy = diff(Ylim)*0.6 # allow to show label when it is on the axe
                          if(coords$y[1] == 0) coords$y = rep(diff(Ylim)*.5, length.out=2) # allow to show line when on the axe
                          panel.text(x=reg$cx, y=reg$cy*diff(Ylim), col=col, labels=lab, pos=4)
                          panel.lines(x=coords$x, y=coords$y*diff(Ylim), col=col)
                        })
                      })
    }
  } else {
    # define legend
    KEY = list("text"=list(rev(displayed_n)),
               "cex"=lt$add.text$cex * 0.5,
               "points"=list(col = sapply(P[rev(displayed_n)], FUN=function(p) p[c("color","lightModeColor")][[color_mode]]),
                             pch = sapply(P[rev(displayed_n)], FUN=function(p) p$style)))
    # identify groups
    groups = NULL
    if(nrow(D) > 0) if(type == "scatter") if(precision=="light") {
      if(length(displayed_n) > 1) {
        groups=apply(as.data.frame(D[,displayed_n]), 1, FUN=function(x) {
          tmp = which(x)[1]
          if(is.na(tmp)) return(NA)
          return(displayed_n[which(x)[1]])
        })
      } else {
        groups = rep(displayed_n, nrow(D))
      }
    }
    # make xyplot
    xtop = NULL
    if(type == "density") {
      args_level = basepop[[1]][["densitylevel"]]
      if((length(args_level) == 0) || (args_level == ""))
        if(!(inherits(trans, what="function") || 
             !inherits(try(suppressWarnings(formals(trans)), silent = TRUE), what="try-error"))) xtop = trans
    }
    foo = xyplot(D[,"y2"] ~ D[,"x2"], auto.key=FALSE,
                 xlim = Xlim, ylim = Ylim,
                 main = main, xlab = xlab, ylab = ylab,
                 xlab.top = xtop,
                 groups=groups,
                 scales =  myScales(x=list(lim = Xlim, "hyper"=Xtrans), y=list(lim = Ylim, "hyper"=Ytrans)),
                 panel = function(x, y, groups=NULL, subscripts, ...) {
                   if(any(c("panel","both")%in%add_key)) if(type=="scatter") pan_key(key=c(KEY,"background"="lightgrey","alpha.background"=0.8), x = 0.02)
                   if(type == "density") {
                     colramp=colorRampPalette(colConv(basepop[[1]][c("densitycolorsdarkmode","densitycolorslightmode")][[color_mode]]))
                     args_level = basepop[[1]][["densitylevel"]]
                     if((length(args_level) != 0) && (args_level != "")) {
                       col = densCols(x=structure(x, features=attr(obj$input$data,"features")), y=y, colramp=colramp, nbin=nbin, transformation=function(x) {x})
                       args_level=strsplit(args_level,split="|",fixed=TRUE)[[1]]
                       fill = args_level[1] == "true"
                       dolines = args_level[2] == "true"
                       nlevels = as.integer(args_level[3])
                       lowest = as.numeric(args_level[4])
                       z  = attr(col, "matrix")
                       zz = as.vector(z$fhat)
                       dd = attr(col, "density")
                       at_probs = seq(from=lowest, to=1, length.out=1+nlevels+(lowest!=0))
                       at = quantile(dd, probs = at_probs, na.rm = TRUE, names = FALSE)
                       low=0
                       if(lowest != 0) {
                         low = at[1]; at = at[-1]
                         panel.xyplot(x=x[dd<low], y=y[dd<low], pch=".", col=c("white","black")[color_mode]) #col[dd<low])
                       }
                       contour_cols = colramp(length(at))
                       if(fill) {
                         # FIXME when filled the graph is pixelated
                         panel.levelplot(x=rep(z$x1, times=nbin), y=rep(z$x2, each=nbin), z=zz,
                                         at=c(low, at), col.regions=c(NA,colramp(length(at)-1)),
                                         subscripts=seq_along(as.vector(zz)),
                                         border = "transparent", border.lty = 1, border.lwd = 0,
                                         region=TRUE, contour=FALSE)
                       }
                       if(dolines) {
                         lines = contourLines(x=z$x1, y=z$x2, z=z$fhat, nlevels=length(at), levels=at)
                         if(fill) contour_cols = rep(c("white","black")[color_mode], length(lines))
                         lapply(1:length(lines), FUN = function(i_l) do.call(panel.lines, args = c(lines[i_l], list(col=contour_cols[lines[[i_l]]$level == at]))))
                       }
                     } else {
                       col = densCols(x=structure(x, features=attr(obj$input$data,"features")), y=y, colramp=colramp, nbin=nbin, transformation=trans)
                       panel.xyplot(x=x,y=y,pch=".", col=col)
                     }
                   }
                   if(type == "scatter") {
                     if(is.null(groups[subscripts])) {
                       panel.xyplot(x=x[1], y=y[1], pch="", alpha=0)
                     } else {
                       by(data.frame("x"=x,"y"=y,"g"=groups[subscripts], stringsAsFactors=FALSE), groups[subscripts], FUN=function(d) {
                         disp = d$g[1]
                         panel.xyplot(x=d$x, y=d$y, pch=P[[disp]]$style, col = P[[disp]][c("color","lightModeColor")][[color_mode]])
                       })
                     }
                   }
                   lapply(reg_n, FUN=function(r) {
                     reg = R[[r]]
                     k = reg[c("color","lightcolor")][[color_mode]]
                     coords = reg[c("x","y")]
                     coords$x = applyTrans(coords$x, trans_x)
                     reg$cx = applyTrans(reg$cx, trans_x)
                     coords$y = applyTrans(coords$y, trans_y)
                     reg$cy = applyTrans(reg$cy, trans_y)
                     if(reg$type=="rect") {
                       coords$x=c(coords$x[1],coords$x[1],coords$x[2],coords$x[2])
                       coords$y=c(coords$y[1],coords$y[2],coords$y[2],coords$y[1])
                     }
                     if(reg$type=="oval") {
                       coords = toEllipse(coords)
                     }
                     lab =  trunc_string(reg$label, trunc_labels)
                     panel.text(x=reg$cx, y=reg$cy, col=k, labels=lab, pos=4)
                     panel.polygon(x=coords$x, y=coords$y, border=k, col="transparent", lwd=1, lty=1)
                   })
                 })
    if(nrow(D) > 0) if(precision=="full") if(type == "scatter") for(l in length(displayed_n):1) {
      disp = displayed_n[l]
      if(any(D[,disp])) { # adds layer only if there is at least one point
        tmp = xyplot(D[,"y2"] ~ D[,"x2"], pch = P[[disp]]$style, col = P[[disp]][c("color","lightModeColor")][[color_mode]], subset = D[,disp])
        foo = foo + as.layer(tmp)
      }
    }
  }
  if(any(c("global","both")%in%add_key)) foo = update(foo, key=KEY)
  foo = update(foo, par.settings = lt)
}

#' @title `IFC_plot` Statistics Extraction
#' @name plot_stats
#' @description Helper to extract `IFC_plot` statistics.
#' @param obj an object of class `IFC_plot` as created by \code{\link{plotGraph}}.
#' @keywords internal
plot_stats=function(obj) {
  # check obj is `IFC_plot`
  assert(obj, cla = "IFC_plot")
  xy_subset = obj$input$subset
  R = obj$input$regions
  D = obj$input$data[xy_subset,,drop=FALSE]
  basepop = obj$input$base
  Xlim = obj$input$xlim
  Ylim = obj$input$ylim
  Xtrans = obj$input$trans_x
  Ytrans = obj$input$trans_y
  trans_x = parseTrans(Xtrans)
  trans_y = parseTrans(Ytrans)
  reg_n = names(R)
  displayed = obj$input$displayed
  displayed_n = names(displayed)
  displayed_o = obj$input$order
  type = obj$input$type
  
  base_n = unlist(lapply(basepop, FUN=function(x) x$name))
  base_o = sapply(base_n, FUN=function(x) which(displayed_n%in%x))
  base_n = base_n[order(base_o)]
  graph_n = obj$input$graphical
  shown_n = setdiff(rev(displayed_n), c(base_n, graph_n))
  
  stats = NULL
  if(type %in% c("count","percent")) {
    coln_stats = c("count","perc","Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
    stats = structure(matrix(numeric(), ncol = length(coln_stats), nrow = 0), dimnames = list(character(), coln_stats))
    base_s = lapply(base_n, FUN=function(d) {
      np = sum(D[,d])
      if(np == 0) return(structure(rep(NA, length(coln_stats)), names = coln_stats))
      p = c("count"=np, "perc"=100, summary(na.omit(D[D[,d],"x1"])))
    })
    kids_s = lapply(shown_n, FUN=function(s) {
      do.call(what = "rbind", args = lapply(base_n, FUN=function(d) {
        np = sum(D[,d])
        if(np == 0) return(structure(rep(NA, length(coln_stats)), names = coln_stats))
        isin = D[,d] & D[,s]
        n = sum(isin)
        c("count"=n, "perc"=n/np*100, summary(na.omit(D[isin,"x1"])))
      }))
    })
    kids_r = lapply(reg_n, FUN=function(r) {
      do.call(what = "rbind", args = lapply(base_n, FUN=function(d) {
        alg = 3
        reg = R[[r]]
        coords = reg["x"]
        coords$x = applyTrans(coords$x, trans_x)
        np = sum(D[,d])
        if(np == 0) return(structure(rep(NA, length(coln_stats)), names = coln_stats))
        isin = D[D[,d],"x2"]
        isin = (isin >= min(coords$x)) & (isin <= max(coords$x))
        n = sum(isin)
        c("count"=n, "perc"=n/np*100, summary(na.omit(D[isin,"x1"])))
      }))
    })
    stats = do.call(what=rbind, args=c(base_s, kids_s, kids_r))
    rnames = base_n
    if(length(reg_n) > 0) rnames = c(rnames, unlist(t(sapply(base_n, FUN = function(b) {if(b == "All") {graph_n} else {paste(reg_n, b, sep = " & ") }}))))
    rownames(stats) = rnames
    colnames(stats) = c(coln_stats[1:2], paste0("x-",coln_stats[3:8]))
  } else {
    coln_stats = c("count","perc","Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
    stats = structure(matrix(numeric(), ncol = length(coln_stats), nrow = 0), dimnames = list(character(), coln_stats))
    base_s = lapply(base_n, FUN=function(d) {
      np = sum(D[,d])
      if(np == 0) return(structure(rep(NA, length(coln_stats)), names = coln_stats))
      p = c("count"=np, "perc"=100, summary(na.omit(D[D[,d],"x1"])), summary(na.omit(D[D[,d],"y1"])))
    })
    kids_s = lapply(shown_n, FUN=function(s) {
      do.call(what = "rbind", args = lapply(base_n, FUN=function(d) {
        np = sum(D[,d])
        if(np == 0) return(structure(rep(NA, length(coln_stats)), names = coln_stats))
        isin = D[,d] & D[,s]
        n = sum(isin)
        c("count"=n, "perc"=n/np*100, summary(na.omit(D[isin,"x1"])), summary(na.omit(D[isin,"y1"])))
      }))
    })
    kids_r = lapply(reg_n, FUN=function(r) {
      alg = 1
      reg = R[[r]]
      coords = reg[c("x","y")]
      coords$x = applyTrans(coords$x, trans_x)
      coords$y = applyTrans(coords$y, trans_y)
      if(reg$type=="oval") alg = 3
      if(reg$type=="rect") alg = 2
      do.call(what = "rbind", args = lapply(base_n, FUN=function(d) {
        np = sum(D[,d])
        if(np == 0) return(structure(rep(NA, length(coln_stats)), names = coln_stats))
        isin = cpp_pnt_in_gate(pnts = cbind(D[D[,d],"x2"],D[D[,d],"y2"]), gate = cbind(coords$x,coords$y), algorithm = alg)
        n = sum(isin)
        c("count"=n, "perc"=n/np*100, summary(na.omit(D[isin,"x1"])), summary(na.omit(D[isin,"y1"])))
      }))
    })
    stats = do.call(what=rbind, args=c(base_s, kids_r, kids_s))
    rnames = base_n
    if(length(reg_n) > 0) rnames = c(rnames, unlist(t(sapply(base_n, FUN = function(b) {if(b == "All") {reg_n} else {paste(reg_n, b, sep = " & ") }}))))
    if(length(shown_n) > 0) rnames = c(rnames, unlist(sapply(shown_n, FUN = function(s) paste(base_n, s, sep = " & "))))
    rownames(stats) = rnames
    colnames(stats) = c(coln_stats[1:2], paste0("x-",coln_stats[3:8]), paste0("y-",coln_stats[9:14]))
  }
  as.table(stats)
}

#' @title IFC Graph Adjustment
#' @name adjustGraph
#' @description Helper to readjust `IFC_data` graphs in case of missing feature, region, population.
#' @param obj an object of class `IFC_data` extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param selection when provided, indices of desired graphs.\cr
#' Note that indices are read from left to right, from top to bottom. 
#' @param adjust_graph whether to try to adjust graph when possible. Default is TRUE.
#' @param ... other arguments to be passed.
#' @keywords internal
adjustGraph=function(obj, selection, adjust_graph = TRUE, ...) {
  dots = list(...)
  assert(obj, cla = "IFC_data")
  G = obj$graphs
  N = names(G)
  if(length(G) > 0) {
    if(missing(selection)) {
      selection = 1:length(G)
    } else {
      if(!all(selection%in%(1:length(G)))) stop("'selection' refers to graph absent from 'obj'")
    }
    foo = lapply(1:length(G), FUN = function(i_graph) {
      g = G[[i_graph]]
      if(i_graph %in% selection) {
        # check if x axis is present in obj
        if(!(g$f1 %in% names(obj$features))) return(NULL)
        # check if y axis is present in obj
        if(g$type != "histogram") if(!(g$f2 %in% names(obj$features))) return(NULL)
        # check that at least one base pop will be plot
        tmp = sapply(g$BasePop, FUN = function(p) p$name %in% names(obj$pops))
        if(!adjust_graph) if(!all(tmp)) return(NULL)
        if(!any(tmp)) return(NULL)
        # remove BasePop not present in obj
        g$BasePop = g$BasePop[tmp]
        # remove GraphRegion not found in obj
        if(length(g$GraphRegion) !=0 && length(g$GraphRegion[[1]]) != 0) {
          tmp = sapply(g$GraphRegion, FUN = function(r) r$name %in% names(obj$regions))
          if(!adjust_graph) if(!all(tmp)) return(NULL)
          g$GraphRegion = g$GraphRegion[tmp]
        }
        # remove ShownPop not found in obj
        if(length(g$ShownPop) != 0 && length(g$ShownPop[[1]]) != 0) {
          tmp = sapply(g$ShownPop, FUN = function(p) p$name %in% names(obj$pops))
          if(!adjust_graph) if(!all(tmp)) return(NULL)
          g$ShownPop = g$ShownPop[sapply(g$ShownPop, FUN = function(p) p$name %in% names(obj$pops))]
        }
        # rebuild Graph, mainly to recompute order
        g = do.call(what = buildGraph, args = g[!grepl("order", names(g))])
        if(inherits(x = g, what = "try-error")) return(NULL)
        # try to draw the graph
        drawable = plotGraph(obj = obj, graph = g, draw = FALSE, stats_print = FALSE)
        if(inherits(x = drawable, what = "try-error")) return(NULL)
      }
      return(g)
    })
    names(foo) = N
    bar = foo[sapply(foo, FUN = function(x) length(x) > 0)]
    class(bar) = class(obj$graphs)
    obj$graphs = bar
    return(obj)
  } else {
    if(!missing(selection) && (length(selection) != 0)) stop("'selection' refers to graph absent from 'obj'")
    return(obj)
  }
}
