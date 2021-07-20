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
calcDensity <- getFromNamespace(x = ".smoothScatterCalcDensity", ns = "grDevices")

#' @title Colors for Smooth Density Plots
#' @description Helper to map density to colors
#' @source derived from \pkg{grDevices} R Core Team, Florian Hahne at FHCRC, originally
#' @keywords internal
densCols=function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("blue","green","red")), transformation=function(x) (asinh(x))) {
  x_features = attr(x, "features")
  xy <- xy.coords(x, y)
  select <- is.finite(xy$x) & is.finite(xy$y)
  x <- cbind(xy$x, xy$y)[select, ]
  map <- calcDensity(x, nbin, bandwidth)
  mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
  xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
  ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
  dens <- map$fhat[cbind(xbin, ybin)]
  check_fun <- try(suppressWarnings(formals(transformation)), silent = TRUE)
  if(!inherits(check_fun, what="try-error")) {
    if("features" %in% names(check_fun)) {
      args = list(features = x_features)
    } else {
      args = list(x = dens)
    }
    dens <- do.call(what = transformation, args = args)  
  } else {
    dens = x_features
  }
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

#' @title Axis Constructor
#' @name base_axis_constr
#' @description Helper to rescale and label axes when linlog / asinh transformation is used.
#' @param lim vector of length 2 of axis extents.
#' @param trans transformation applied. Defaut is "P".
#' @param nint positive integer value indicating (approximately) the desired number of intervals. Default is 10.
#' @keywords internal
base_axis_constr = function(lim, trans = "P", nint = 10) {
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
          # log = log10(at[at %in% pretty_lab])
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
    keep = (ticks_at >= lim[1]) & (ticks_at <= lim[2])
    dup = duplicated(ticks_at)
    ticks_at = ticks_at[!dup & keep]
    ticks_lab = ticks_lab[!dup & keep]
    if(length(ticks_at) < 2) {
      at = signif(seq(from = lim[1], to = lim[2], length.out = nint), 3)
      ticks_at = applyTrans(at, trans_)
      ticks_lab = at
    }
    ord = order(ticks_at)
    return(list("at" = ticks_at[ord], "labels" = ticks_lab[ord]))
  } else {
    at = axisTicks(lim, log = FALSE, nint = nint)
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
base_hist_constr = function(x, type, br, normalize, fill, smooth, lwd, lty, col, alpha, border, include.lowest=TRUE, right=TRUE){
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

#' @title IFC Graph Conversion to 'base' Plot
#' @name convert_to_baseplot
#' @description Helper to convert `IFC_plot` to 'base' plot.
#' @param obj an object of class `IFC_plot` as created by \code{\link{plotGraph}}.
#' @keywords internal
convert_to_baseplot = function(obj) {
  # variables for future use
  n_ticks = 10
  pkg = "base"
  
  # check obj is `IFC_plot`
  assert(obj, cla = "IFC_plot")
  
  # short names
  D = obj$input$data
  displayed = obj$input$displayed
  basepop = obj$input$base
  Xlim = obj$input$xlim
  Ylim = obj$input$ylim
  disp_n = names(displayed)
  # draw plot
  if(obj$input$type %in% c("percent", "count")) {
    # 1D
    br = do.breaks(Xlim, obj$input$bin)
    hist(D[, "x2"], xlim = Xlim, ylim = Ylim, 
         main = trunc_string(obj$input$title, obj$input$trunc_labels), 
         xlab = trunc_string(obj$input$xlab, obj$input$trunc_labels),
         ylab = trunc_string(obj$input$ylab, obj$input$trunc_labels),
         cex.lab = obj$plot$par.settings$par.xlab.text$cex,
         cex.main = obj$plot$par.settings$par.main.text$cex,
         cex.axis = obj$plot$par.settings$axis.text$cex,
         col = "transparent", border = "transparent", freq = obj$input$type == "count",
         breaks = br, axes = FALSE)
    if(length(displayed) > 0) {
      for(disp in disp_n) {
        if(any(D[,disp]))
          base_hist_constr( D[D[,disp], "x2"], br = br, type = obj$input$type, 
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
    if(obj$input$type == "density") {
      pch=16
      col = "white"
      if((nrow(obj$input$data) > 0) && any(obj$input$subset)) if(nrow(obj$input$data[obj$input$subset,]) > 0) 
        col=densCols(x = structure(obj$input$data$x2[obj$input$subset], features=attr(obj$input$data,"features")),
                     y = obj$input$data$y2[obj$input$subset],
                     colramp=colorRampPalette(colConv(basepop[[1]][c("densitycolorsdarkmode","densitycolorslightmode")][[obj$input$mode]])),
                     nbin=obj$input$bin,
                     transformation=obj$input$trans)
      plot(x = obj$input$data$x2[obj$input$subset],
           y = obj$input$data$y2[obj$input$subset],
           xlim = Xlim , ylim = Ylim ,
           main = obj$input$title,
           xlab = trunc_string(obj$input$xlab, obj$input$trunc_labels),
           ylab = trunc_string(obj$input$ylab, obj$input$trunc_labels),
           cex.lab = obj$plot$par.settings$par.xlab.text$cex,
           cex.main = obj$plot$par.settings$par.main.text$cex,
           cex.axis = obj$plot$par.settings$axis.text$cex,
           pch = pch, col = col, 
           axes = FALSE)
    } else {
      if(obj$input$precision == "full") {
        disp = disp_n[length(displayed)]
        plot(x = obj$input$data[obj$input$data[,disp], "x2"][obj$input$subset],
             y = obj$input$data[obj$input$data[,disp], "y2"][obj$input$subset],
             xlim = Xlim , ylim = Ylim ,
             main = obj$input$title,
             xlab = trunc_string(obj$input$xlab, obj$input$trunc_labels),
             ylab = trunc_string(obj$input$ylab, obj$input$trunc_labels),
             cex.lab = obj$plot$par.settings$par.xlab.text$cex,
             cex.main = obj$plot$par.settings$par.main.text$cex,
             cex.axis = obj$plot$par.settings$axis.text$cex,
             pch = displayed[[disp]]$style, 
             col = displayed[[disp]][c("color","lightModeColor")][[obj$input$mode]],
             axes = FALSE)
        if(length(displayed) > 1) {
          for(disp in rev(disp_n)[-1]) {
            points(x = obj$input$data[obj$input$data[,disp], "x2"][obj$input$subset],
                   y = obj$input$data[obj$input$data[,disp], "y2"][obj$input$subset],
                   pch = displayed[[disp]]$style, 
                   col = displayed[[disp]][c("color","lightModeColor")][[obj$input$mode]])
          }
        }
      } else {
        groups = apply(as.data.frame(D[obj$input$subset,disp_n]), 1, FUN=function(x) {
          tmp = which(x)[1]
          if(is.na(tmp)) return(NA)
          return(disp_n[tmp])
        })
        pch = sapply(groups, FUN = function(disp) displayed[[disp]]$style)
        col = sapply(groups, FUN = function(disp) displayed[[disp]][c("color","lightModeColor")][[obj$input$mode]])
        plot(x = obj$input$data$x2[obj$input$subset],
             y = obj$input$data$y2[obj$input$subset],
             xlim = Xlim, ylim = Ylim,
             main = obj$input$title,
             xlab = trunc_string(obj$input$xlab, obj$input$trunc_labels),
             ylab = trunc_string(obj$input$ylab, obj$input$trunc_labels),
             cex.lab = obj$plot$par.settings$par.xlab.text$cex,
             cex.main = obj$plot$par.settings$par.main.text$cex,
             cex.axis = obj$plot$par.settings$axis.text$cex,
             pch = pch, col = col,
             axes = FALSE)
      }
    }
  }
  
  # axis
  if(pkg == "base") {
    x_ticks = base_axis_constr(lim = Xlim, trans = obj$input$trans_x, nint = n_ticks)
    y_ticks = base_axis_constr(lim = Ylim, trans = obj$input$trans_y, nint = n_ticks)
    x_axis = axis(side = 1, at = x_ticks$at, labels = FALSE)
    text(x = x_axis, y = Ylim[1] - diff(Ylim) * 0.07, labels = x_ticks$labels, xpd=TRUE, adj = c(1, 1),
         cex = obj$plot$par.settings$axis.text$cex, cex.axis = obj$plot$par.settings$axis.text$cex, srt=45)
    y_axis = axis(side = 2, at = y_ticks$at, labels = FALSE)
    text(y = y_axis, x = Xlim[1] - diff(Xlim) * 0.07, labels = y_ticks$labels, xpd=TRUE, adj = c(1, 0.5),
         cex = obj$plot$par.settings$axis.text$cex, cex.axis = obj$plot$par.settings$axis.text$cex, las=1)
    box()
  }
  
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
      switch(pkg,
             "lattice" = {
               foo = foo +
                 layer(panel.text(x=reg$cx, y=reg$cy*diff(Ylim), col=k, labels=lab, pos=4, cex=obj$plot$par.settings$add.text$cex)) +
                 layer(panel.lines(x=coords$x, y=coords$y*diff(Ylim),col=k))
             },
             "base" = {
               text(x=reg$cx, y=reg$cy*diff(Ylim), col=k, labels=lab, pos=4, cex=obj$plot$par.settings$add.text$cex)
               polygon(x=coords$x, y=coords$y*diff(Ylim), col = k, border = k)
             })
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
      switch(pkg,
             "lattice" = {
               foo = foo +
                 layer(panel.text(x=reg$cx, y=reg$cy, col=k, labels=lab, pos=4, cex=obj$plot$par.settings$add.text$cex)) +
                 layer(panel.polygon(x=coords$x, y=coords$y, border=k, col="transparent", lwd=1, lty=1))
             },
             "base" = {
               text(x=reg$cx, y=reg$cy, col=k, labels=lab, pos=4, cex=obj$plot$par.settings$add.text$cex) 
               polygon(x=coords$x, y=coords$y, border=k, col="transparent", lwd=1, lty=1)
             })
    }
  }
  if(pkg != "base") plot(foo)
  
  # key
  if(obj$input$type %in% c("percent", "count")) {
    legend("topleft", inset = 0.025, 
           lty = sapply(disp_n, FUN=function(p) c(1,2,3,4,6)[match(basepop[[obj$input$order[p]]]$linestyle,c("Solid","Dash","Dot","DashDot","DashDotDot"))]),
           col = sapply(displayed, FUN=function(p) p[c("color","lightModeColor")][[obj$input$mode]]),
           legend = disp_n, cex = obj$plot$par.settings$add.text$cex * 0.5, bg = "#ADADAD99", pt.cex = 1, bty ="o", box.lty = 0)
  } else {
    legend("topleft", inset = 0.025, 
           pch = sapply(displayed, FUN=function(p) p$style),
           col = sapply(displayed, FUN=function(p) p[c("color","lightModeColor")][[obj$input$mode]]),
           legend = disp_n, cex = obj$plot$par.settings$add.text$cex * 0.5, bg = "#ADADAD99", pt.cex = 1, bty ="o", box.lty = 0)
  }
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
adjustGraph = function(obj, selection, adjust_graph = TRUE, ...) {
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
    if(length(selection) != 0) stop("'selection' refers to graph absent from 'obj'")
    return(obj)
  }
}
