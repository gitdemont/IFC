################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2022 Yohann Demont                                             #
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

################################################################################
#             functions described hereunder are experimental                   #
#              inputs and outputs may change in the future                     #
################################################################################

#' @title Fast 2D plot
#' @description
#' Creates fast 2D plots with Rcpp
#' @param x,y the x and y coordinates for the plot. If 'y' is NULL, it will be same as 'x'.
#' x and y should have same length.
#' @param pch a vector of plotting symbols. Default is "." (resulting in a 1-pixel dot). Allowed are 0 to 20 and ".". It will be repeated along 'x'. 
#' With the exception of ".", NA resulting its coercion to integer(s) from the conversion will be omitted
#' (i.e. points won't be displayed, but their x-y coordinates will account for xlim, ylim range computation when not provided through ...).
#' Everything else (coercible to integer) will result in a dot (a 1-pixel pixel).
#' @param size an integer vector giving the size(s) of the 'pch'. Default is 7. It will be repeated along 'x'. 
#' @param alpha a [0,255] integer. Default is 255.
#' @param col a vector of desired colors of the symbols that will be passed by grDevices::col2rgb('x', alpha = TRUE). Default is "black".
#' If number of colors equals number of points every point will be assigned this color.
#' Otherwise, if color is of length 1 for a single combination of size / pch, all points with this combination will be assigned this color.
#' Finally, if there is only one combination of size / pch and color not equals 1 nor the total number of points,
#' then, colors will be used as a gradient for density (in such case 'blur_size' and 'blur_sd' will be taken in consideration)
#' This only applies when 'force' is FALSE.
#' @param rgba a 4 rows color matrix, with rows being Red, Green, Blue and Alpha and number of columns identical to number of points.\cr
#' /!\ When provided this argument will take precedence over 'col' and 'alpha'.
#' @param force whether to force scatter instead of density when multiple 'col' are provided
#' @param draw whether to draw to plot (when TRUE), or to image only (when FALSE). Default is TRUE.
#' @param new whether a new plot should be created, only applies when 'draw' is TRUE. Default is is.null(bg_).
#' If FALSE, the current plot will be used to draw points.
#' @param interpolate whether to use linear interpolation, only applies when 'draw' is TRUE. Default is FALSE.
#' @param width the desired width of the raster Default is 512. It only applies when draw is FALSE.
#' @param height the desired height of the raster Default is 512. It only applies when draw is FALSE.
#' @param pntsonedge whether points outside of plotting region should be bounded on the edge. Default is FALSE to clip points.
#' @param blur_size (for density) an integer controlling the size of the blurring gaussian kernel. Default is 9.
#' @param blur_sd (for density) a double controlling the sd of the blurring gaussian kernel. Default is 3.
#' @param bg_ an `rasterplot` object as returned by rasterplot() that will be used to add points to. Default is NULL.
#' If provided it will have to be compatible with current drawing size or 'width' and 'height' when 'draw' is FALSE.
#' @param bg_map whether to use 'bg_' when provided to compute points coordinates. Default is TRUE.
#' This allows to get same "user" to "pixel" coordinates conversion as the one used to create 'bg_'.
#' @param ... other arguments to pass to graphics::plot().
#' For example, providing xlim and/or ylim will controls if point will be shown or not
#' @details some examples:\cr
#' set.seed(2)\cr
#' n_points = 1e7; n_clusters = 5\cr
#' x = c(t(sapply(1:5, FUN = function(i) rnorm(n_points / n_clusters, mean = sample(-2:2, size = 1), sd = 1/sample(1:10, 1)))))\cr
#' y = c(t(sapply(1:5, FUN = function(i) rnorm(n_points / n_clusters, mean = sample(-2:2, size = 1), sd = 1/sample(1:10, 1)))))\cr
#' # plot points\cr
#' rasterplot(x = x, y = y, col = "black")\cr
#' # generate img\cr
#' rasterplot(x = x, y = y, col = "black", draw = FALSE)\cr
#' # plot multiple shapes\cr
#' rasterplot(x = x, y = y, pch = c(3,5,9,10,2))\cr
#' # plot multiple shapes + colors\cr
#' bg_ = rasterplot(x = x, y = y, pch = c(3,5,9,10,2), col = c("plum", "green", "indianred", "blue", "black"))\cr
#' # addition of new points to an already drawn background, it a kind of points(...)\cr
#' rasterplot(x = x[1:1e5], y = y[1:1e5], col = "black", bg_ = bg_, bg_map = TRUE)\cr
#' # plot 1 shape  + multiple colors\cr
#' rasterplot(x = x, y = y, pch = ".", col = c("plum", "green", "indianred", "blue", "black"), force = TRUE)\cr
#' # density\cr
#' rasterplot(x = x, y = y, pch = 20, size = 7, draw = TRUE, col = colorRampPalette(c("blue", "green", "red"))(100))\cr
#' # density with limits\cr
#' rasterplot(x = x, y = y, draw = TRUE, xlim = c(0, 1.5), pntsonedge = FALSE, col = colorRampPalette(c("blue", "green", "red"))(100))\cr
#' # density with limits + computation on drawn points only\cr
#' rasterplot(x = x, y = y, draw = TRUE, xlim = c(0, 1.5), pntsonedge = TRUE, col = colorRampPalette(c("blue", "green", "red"))(100))\cr
#' # using rgba\cr
#' col = c("plum", "green", "indianred", "blue", "black")\cr
#' rgba = col2rgb(col, alpha = TRUE)\cr
#' rgba = t(apply(rgba, 1, FUN = function(x) rep(x, length.out = n_points)))\cr
#' rasterplot(x = x, y = y, pch = ".", rgba = rgba, draw = TRUE)
#' @return an [0, 255] integer array of (height, width, 4) of class `rasterplot`
#' @keywords internal
rasterplot = function(x, y = NULL, 
                      pch = ".", size = 7,
                      alpha = 255, col = "black", rgba = NULL, force = FALSE, # parameters for colors
                      draw = TRUE,
                      new = is.null(bg_), interpolate = FALSE,                # only when draw == TRUE
                      width = 512, height = 512,                              # only when draw == FALSE
                      pntsonedge = FALSE, 
                      blur_size = 9, blur_sd = 3,                             # only for density
                      bg_ = NULL, bg_map = TRUE, ...) {
  dots = list(...)
  if(length(y) == 0) y = x
  xlim = dots$xlim
  if(length(xlim) == 0) xlim = cpp_fast_range(x)
  ylim = dots$ylim
  if(length(ylim) == 0) ylim = cpp_fast_range(y)
  xlab = dots$xlab
  if(length(xlab) == 0) xlab = "x"
  ylab = dots$ylab
  if(length(ylab) == 0) ylab = "y"
  main = dots$main
  if(length(main) == 0) main = "Raster Plot"
  if(missing(rgba)) {
    if(length(nrow(rgba) != 4) && (ncol(rgba) != nrow(d))) stop("when provided 'rgba' should be a 4 rows matrix of number of columns identical to x length")
  } else {
    if(length(col) == 0) stop("bad 'col' specification")
  }
  alpha = as.integer(alpha)
  if(length(alpha) != 1) stop("'alpha' should be of length 1")
  if((alpha < 0) || (alpha > 255)) stop("'alpha' should be a [0,255] integer")
  
  dots = dots[setdiff(names(dots), c("xlim", "ylim", "xlab", "ylab", "main"))]
  pch[pch == "."] <- 27
  pch = suppressWarnings(as.integer(pch))
  has_bg = !missing(bg_) && (length(bg_) != 0)
  if(has_bg && !inherits(bg_, "rasterplot")) stop("when provided 'bg_' should be of class `rasterplot`")
  cex = par("cex"); if(length(dots$cex) != 0) cex = dots$cex
  if(cex < 0) stop("'cex' should be positive numeric")
  lwd = par("lwd"); if(length(dots$lwd) != 0) lwd = dots$lwd
  if(lwd < 0) stop("'lwd' should be positive integer")
  size = size * cex
  size[size < 1] <- 1
  
  # create empty plot
  if(draw) {
    if(new) do.call(what = plot, 
                    args = c(list(x = quote(x[1]), y = quote(y[1]), col = "transparent",
                                  xlim = xlim, ylim = ylim,
                                  xlab = xlab, ylab = ylab,
                                  main = main),
                             dots))
    if(has_bg && bg_map) {
      coordmap = attr(bg_, "coordmap")
    } else {
      coordmap = get_coordmap_adjusted()
    }
  } else {
    if(has_bg && bg_map) {
      coordmap = attr(bg_, "coordmap")
    } else {
      coordmap = list(domain = c(list(bottom = ylim[1], top = ylim[2]),
                               list(left = xlim[1], right = xlim[2])),
                    range = c(list(right = 0, left = width - 1),
                              list(top = 0, bottom = height - 1)),
                    width = width,
                    height = height,
                    ratio = list(x = 1, y = 1))
    }
  }
  
  # setup data
  one_group = (length(pch) == 1) && (length(size) == 1)
  if(one_group) one_group = !force
  if(one_group) { # only one combination of size / pch
    data = list(list(size = size,
                     pch = pch,
                     lwd = lwd, 
                     coords = coord_to_px(coord=data.frame(x = x, y = y), coordmap = coordmap, pntsonedge = pntsonedge),
                     blur_size = blur_size,
                     blur_sd = blur_sd))
    if(missing(rgba)) {
      data[[1]]$col = col2rgb(col, alpha = TRUE)
      data[[1]]$col[4,] <- alpha
    } else {
      data[[1]]$col = rgba
    }
  } else {
    if(missing(rgba)) { # we draw every combinations with its color argument (one only !)
      d = data.frame(x = x, y = y, pch = pch, size = size, col = col) 
      g = group(d[, 3:5], keepNAlevels = FALSE)
      data = lapply(seq_along(g), FUN = function(i) {
        col_ = col2rgb(d$col[g[[i]][1]], alpha = TRUE)
        col_[4, ] <- alpha
        list(size = d$size[g[[i]][1]],
             pch = d$pch[g[[i]][1]],
             lwd = lwd,
             col = col_,
             coords = coord_to_px(coord=data.frame(x = d$x[g[[i]]], y = d$y[g[[i]]]), coordmap = coordmap, pntsonedge = pntsonedge),
             blur_size = blur_size,
             blur_sd = blur_sd)
      })
    } else { # with rgba we have a color for each single point
      d = data.frame(x = x, y = y, pch = pch, size = size)
      if((length(pch) == 1) && (length(size) == 1)) {
        g = list(seq_along(x))
      } else {
        g = group(d[, 3:4], keepNAlevels = FALSE) 
      }
      data = lapply(seq_along(g), FUN = function(i) {
        list(size = d$size[g[[i]][1]],
             pch = d$pch[g[[i]][1]],
             lwd = lwd,
             col = rgba[, g[[i]], drop = FALSE],
             coords = coord_to_px(coord=data.frame(x = d$x[g[[i]]], y = d$y[g[[i]]]), coordmap = coordmap, pntsonedge = pntsonedge),
             blur_size = blur_size,
             blur_sd = blur_sd)
      })
    }
  }
  
  # plot it / return it
  if(draw) {
    # create raster
    width = coordmap$width
    height = coordmap$height
    img = cpp_raster(width = width, height = height, data, bg_)
    
    # subset img to drawing region
    usr = unlist(recursive = FALSE, use.names = FALSE, coordmap$domain)
    lims = round(c(coord_to_px(coord = data.frame(x = usr[1:2], y = usr[3:4]),
                               coordmap = coordmap,
                               pntsonedge = F)) + c(1,1,0,0))
    
    # add image to plot background
    if(!identical(get_coordmap_adjusted(), coordmap)) {
      text(x = graphics::grconvertX(0.5, "npc", "user"),
           y = graphics::grconvertY(0.5, "npc", "user"),
           labels = "you should not modify graphic device while drawing raster",
           col = "red")
      img = NULL
    } else {
      rasterImage(image = cpp_as_nativeRaster(img[lims[3]:lims[4], lims[1]:lims[2],]),
                xleft = usr[1], xright = usr[2], ybottom = usr[4], ytop = usr[3],
                interpolate = interpolate)
    }
    return(invisible(structure(img, "coordmap" = coordmap, class = "rasterplot")))
  } else {
    return(invisible(structure(cpp_raster(width = width, height = height, data, bg_), "coordmap" = coordmap, class = "rasterplot")))
  }
}
