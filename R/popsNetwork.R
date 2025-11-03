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

#' @title IFC_pops Network Display
#' @description
#' Builds and displays populations network.
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param hierarchical whether to display network using a hierarchical layout or not. Default is \code{TRUE}, will be forced to \code{TRUE} when \code{'flatten'} is \code{TRUE}.
#' @param color_mode Whether to extract colors from \code{'obj'} in \code{"white"} or \code{"black"} mode. Default is \code{"white"}.
#' @param highlight population to permanently highlight. If found in \code{'obj'}, this population will be displayed with its color. Default is \code{NULL}.
#' @param seed If you provide a seed manually, the layout will be the same every time. Default is \code{NULL}, only apply when \code{'flatten'} is \code{FALSE}.
#' @param direction The direction of the hierarchical layout. Default is \code{'LR'}, will be forced to \code{'UD'} when \code{'flatten'} is \code{TRUE}.\cr
#' The available options are: \code{'UD'}, \code{'DU'}, \code{'LR'}, \code{'RL'}. To simplify: up-down, down-up, left-right, right-left.
#' @param weighted whether to scale population's node size according to count. Default is \code{TRUE}, only apply when \code{'flatten'} is \code{FALSE}.
#' @param flatten whether to return a flatten data.frame of class `IFC_network` containing stats, colors, styles and whose row names represent the hierarchy tree. Default is \code{FALSE}.
#' @param ... other argument to be passed.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#'   popsNetwork(obj = daf)
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return a \pkg{visNetwork} object or a data.frame of class `IFC_network` if \code{'flatten'} is \code{TRUE}.
#' @export
popsNetwork = function(obj, hierarchical=TRUE, color_mode="white", highlight=NULL, seed=NULL, direction="LR", weighted=TRUE, flatten=FALSE, ...) {
  dots = list(...)
  if(missing(obj)) stop("'obj' can't be missing")
  if(!("IFC_data"%in%class(obj))) stop("'obj' is not of class `IFC_data`")
  if(length(obj$pops)==0) stop("please use argument 'extract_features' = TRUE with ExtractFromDAF() or ExtractFromXIF() and ensure that features were correctly extracted")
  flatten = as.logical(flatten); assert(flatten, len=1, alw=c(TRUE,FALSE))
  assert(color_mode, len=1, alw=c("white","black"))
  if(flatten) {
    if(length(obj$stats) == 0) obj$stats = get_pops_stats(obj$pops, as.integer(obj$description$ID$objcount))
    tree = popsTree(obj$pops)
    name = names(tree)
    styl = popsStyle(obj)
    d = cbind(as.data.frame(obj$stats[name, 3:5]),
              style = sapply(name, FUN = function(n) map_style(obj$pops[[n]]$style, toR = FALSE)),
              color = sapply(name, FUN = function(n) obj$pops[[n]]$color),
              lightModeColor = sapply(name, FUN = function(n) obj$pops[[n]]$lightModeColor))
    rownames(d) = sapply(seq_along(tree), simplify = TRUE, FUN=function(i) {
      paste0(collapse="",
             paste0(tree[[i]],collapse=""),
             paste0(ifelse(i==1,"","\u202f"),styl[[name[i]]],"\u00a0",name[i]))
    })
    attr(d, "highlight") = highlight
    attr(d, "color_mode") = color_mode
    class(d) = unique(c("IFC_network",class(d)))
    return(d)
  }
  if(!is.null(seed)) {seed=na.omit(as.integer(seed)); assert(seed, len=1, typ="integer")}
  hierarchical = as.logical(hierarchical); assert(hierarchical, len=1, alw=c(TRUE,FALSE))
  weighted = as.logical(weighted); assert(weighted, len=1, alw=c(TRUE,FALSE))
  color_mode=which(c("black","white")%in%color_mode)
  assert(direction, len=1, alw=c("UD", "DU", "LR", "RL"))

  # determine levels for hierarchical layout
  lev = popsGetLevels(obj$pops)
  size = table(lev)
  size = c(length(size), max(size))
  width = NULL
  height = NULL
  session = NULL
  if(with_seed(requireNamespace("shiny", quietly = TRUE), NULL)) session = shiny::getDefaultReactiveDomain()
  if(hierarchical && (length(session) != 0)) {
    width = size[1] * 75
    height = size[2] * 75
  }
  # create nodes
  nodes = as.data.frame(t(sapply(obj$pops, FUN=function(p) {
    typ = c("Base","Boolean","Graphical","Tagged")[c("B","C","G","T")%in%p$type]
    col = rgb(t(checkColor(p[c("color","lightModeColor")][[color_mode]])), maxColorValue = 255)
    val = ifelse(weighted, sum(p$obj), 1)
    tit = paste0("<em>", p$name, "</em><br>", typ, ", n: <b>", val, "</b>, %All: <b>", round(100*sum(p$obj)/length(p$obj),2), "</b>")
    if(typ=="Boolean") tit = paste0(tit, "<br>def: ", paste0(p$split,collapse=" "))
    if(typ=="Graphical") {
      tit = paste0(tit, ", %", p$base, ": <b>", round(100*sum(p$obj)/sum(obj$pops[[p$base]]$obj),2), "</b>")
      sync_typ = sync_type(obj$regions[[p$region]])
      sync_num = sync_part(obj$regions[[p$region]], sync_typ)
      type = obj$regions[[p$region]]$type
      if(sync_num != "") type = paste0(type, " ", sync_typ, sync_num)
      tit = paste0(tit, "<br><b>", type, "</b><br><b>x</b>: ", p$fx)
      if(!is.null(p$fy)) tit = paste(tit, "<br><b>y</b>: ", p$fy)
      tit = paste0(tit, "<table> <tr> <th>x</th> <th>y</th> </tr>")
      tit = paste0(tit, paste("<tr>", apply(do.call(cbind, args = obj$regions[[p$region]][c("x","y")]), 1, FUN=function(r) paste("<td>",r,"</td>",collapse=" ")), "</tr>", collapse = "\n"), "</table>")
    }
    c("label"=p$name, "id"=p$name, "type"=typ, "title"=tit, "value"=val, 
      "color.background"=ifelse(p$name%in%highlight, col, "lightgrey"),
      "color.border"=ifelse(!p$name%in%highlight, col, "lightgrey"),
      "color.highlight.background"=col, 
      "color.highlight.border"=ifelse(!p$name%in%highlight, col, "lightgrey"),
      "color.hover.background"=col,
      "color.hover.border"=ifelse(!p$name%in%highlight, col, "lightgrey"))
  })), stringsAsFactors = FALSE)
  nodes$level = lev
  # create edges
  edges = lapply(obj$pops, FUN=function(p) {
    col = rgb(t(checkColor(p[c("color","lightModeColor")][[color_mode]])), maxColorValue = 255)
    if(p$type=="C") {
      from = p$names
    } else {
      from = p$base
    }
    return(data.frame(to=rep(p$name, length(from)), from=from, arrows="to", stringsAsFactors = FALSE))
  })
  visNetwork(nodes, do.call(rbind, edges[-1]), main = list(text="Populations network"), width = width, height = height) %>%
    visPhysics(enabled = TRUE, stabilization = FALSE) %>%
    visOptions(highlightNearest = list(enabled=TRUE, hover=TRUE, labelOnly=FALSE, degree=list(from=1,to=1), algorithm="hierarchical"), selectedBy="type") %>%
    visLayout(randomSeed = seed) %>%
    visInteraction(navigationButtons = TRUE) %>%
    visHierarchicalLayout(direction = direction, enabled = hierarchical) 
}

#' @title Print Network
#' @description
#' Print `IFC_network`
#' @param x object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @keywords internal
print.IFC_network <- function(x, ...) {
  dots  = list(...)
  q = c("count", "%_par", "%_tot")
  r = rownames(x)
  n = length(rownames(x))
  max = dots$max
  if(length(dots$justify) == 0) dots$justify = "right"
  if(length(x) == 0L) {
    cat(sprintf("data frame with 0 columns and %d row%s", n), ifelse(n > 1, "s", ""), "\n", sep = "")
  } else if(n == 0L) {
    print.default(names(x), quote = FALSE)
    cat(gettext("<0 rows> (or 0-length row.names)\n"))
  } else {
    if(is.null(max)) max <- getOption("max.print", 99999L)
    if(!is.finite(max)) stop("invalid 'max' / getOption(\"max.print\"): ", max)
    mod = attr(x, "color_mode")
    fg = tolower(mod)
    omit <- (n0 <- max%/%length(x[,1:3])) < n
    m <- as.matrix(do.call(format, c(list(if(omit) x[seq_len(n0), 1:3, drop = FALSE] else x[, 1:3, drop = TRUE]), dots)))
    if(nrow(m)) {
      m <- apply(m, 2, FUN = function(x) {
        w = suppressWarnings(max(nchar(x)))
        sapply(x, FUN = function(s) resize_string(trimws(s), w, dots$justify))
      })
      w = nchar(m[1, ])
      cat(sapply(1:3, FUN = function(i) resize_string(q[i], w[i], dots$justify)), "col", "tree", "\n", sep = " ")
    } else {
      cat(q[1:3], "col", "tree", "\n", sep = " ")
    }
    sapply(seq_len(nrow(m)), FUN = function(i) {
      tree = strsplit(r[i], split="\u00a0", fixed=TRUE)[[1]]
      nam = paste0(tail(tree,-1),collapse="\u00a0")
      if(any(nam %in% attr(x,"highlight"))) nam = ansi_bold(ansi_underline(nam))
      nod = paste0(head(tree,-1),collapse="\u00a0")
      bg = tolower(x[i, 4 + match(mod, c("black","white"))])
      if(bg == fg) fg = c("white","black")[match(fg, c("black", "white"))]
      sty = x[i, 4]
      cat(m[i, 1:3],
          ansi_colorize(paste0("[",ansi_underline("\u00af"),"]"), bg = map_ansi(bg), fg = map_ansi(fg)),
          paste(nod, ansi_colorize(map_utf8(sty)), nam),
          "\n", sep = " ")
    })
  }
  return(invisible(x))
}
