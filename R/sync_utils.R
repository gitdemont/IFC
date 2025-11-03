################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2025 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2025                                                                  #
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

#' @title Attributes Keeper
#' @description 
#' Helper to keep attributes while using IFC build functions.
#' @param x list to pass to \code{'what'}.
#' @param what a build function.
#' @details all attributes of \code{'x'} will be used with the exception of "names", "dims", and "dimnames".
#' @return result of applying \code{'what'} on \code{'x'} with preserved attributes from \code{'x'}.
#' @keywords internal
keep_attributes <- function(x, what) {
  ans <- do.call(what, x)
  for(k in setdiff(names(attributes(x)),c("names","dim","dimnames"))) { attr(ans, k) <- attr(x, k) }
  ans
}

#' @title Sync Name Decomposition
#' @description
#' Helper that will split text into chunks.
#' @param x string
#' @param name whether to return name only. Default is \code{TRUE}.
#' @param split string used for splitting. Default is \code{"|"}.
#' @details \code{'x'} is supposed to be a string composed of a \code{'name'} followed by \code{'split'} and ending with a string that does not contain \code{'split'}.
#' @return when \code{'x'} is of length 0, \code{""} is returned. Otherwise, the \code{'name'} when \code{'name'} is \code{TRUE} or the \code{'split'} non-containing remaining string when \code{FALSE}.
#' @keywords internal
splits <- function(x, name = TRUE, split = "|") {
  if(length(x) == 0) return("")
  if(length(x) != 1) stop("'x' should be scalar")
  if(length(split) != 1) stop("'split' should be scalar")
  if(!grepl(split, x, fixed = TRUE)) stop("can't find 'split' in 'x'")
  s = strsplit(x, split, fixed = TRUE)[[1]]
  if(name) {
    if(substr(x, nchar(x), nchar(x)) == split) return(substr(x, 1, nchar(x) - 1))
    return(paste0(head(s, -1), collapse = split))
  }
  if(substr(x, nchar(x), nchar(x)) == split) return(character())
  return(tail(s, 1))
}

#' @title Sync Limits
#' @description
#' Sync Limits.
#' @param lims limits. Default is \code{c(-Inf,Inf)}.
#' @return \code{range(lims,na.rm=TRUE)}
#' @keywords internal
sync_lims <- function(lims = c(-Inf, Inf)) {
  return(suppressWarnings(range(lims, na.rm=TRUE)))
}

#' @title Sync Clip
#' @description
#' Helper to clip vector to sync limits.
#' @param x a numeric vector
#' @param lims limits. Default is missing to use \code{sync_lims()}, otherwise it will be \code{sync_lims(lims)}.
#' @details \code{'x'} will be clip to \code{'lims'} excluding \code{'lims'}
#' @keywords internal
sync_clip <- function(x, lims) {
  if(missing(lims)) {lims = sync_lims()} else {lims=sync_lims(lims)}
  return(x[x > min(lims, na.rm = TRUE) & x < max(lims, na.rm = TRUE)])
}

#' @title Sync Type
#' @description
#' Helper to retrieve sync type.
#' @param reg a region, though it can be any R object
#' @param what a string. If provided returned sync type should be identical to \code{'what'}, otherwise \code{""} is returned.
#' @details sync type is the result of \code{substr(attr(reg, "sync"), 1, 4)}.
#' In addition, when this result is \code{"dual"} or \code{"quad"} and \code{'what'} is provided, \code{attr(reg, "sync")}:\cr
#' -should end with \code{"|1"} or \code{"|2"} for \code{"dual"}, otherwise \code{""} is returned.\cr
#' -should end with \code{"|1"}, \code{"|2"}, \code{"|3"} or \code{"|4"} for \code{"quad"} otherwise \code{""} is returned.
#' @return a string containing sync type.
#' @keywords internal
sync_type <- function(reg, what) {
  A = attr(reg, "sync")
  pre = try(splits(A, name = FALSE, split = "|"), silent = TRUE)
  if(inherits(pre, "try-error") || (length(pre) == 0)) return("")
  typ = substring(A, 1, 4)
  if(missing(what)) what = typ
  if(identical(typ, "quad") && identical(what, typ) && any(pre %in% c("1","2","3","4"))) return(typ)
  if(identical(typ, "dual") && identical(what, typ) && any(pre %in% c("1","2"))) return(typ)
  return("")
}

#' @title Sync Part
#' @description
#' Helper to retrieve sync part.
#' @param reg a region, though it can be any R object
#' @param what a string. If provided \code{sync_type(reg)} should be identical to \code{'what'}, otherwise \code{""} is returned.
#' @details sync part is the result of retrieving number after the last \code{"|"} encountered in \code{attr(reg, "sync")}.
#' @return a string containing sync part.
#' @keywords internal
sync_part <- function(reg, what) {
  A = attr(reg, "sync")
  pre = try(splits(A, name = FALSE, split = "|"), silent = TRUE)
  if(inherits(pre, "try-error") || (length(pre) == 0)) return("")
  if(missing(what)) return(pre)
  if(identical(what, sync_type(reg, what))) return(pre)
  return("")
}

#' @title Sync Name
#' @description
#' Helper to retrieve sync name.
#' @param reg a region, though it can be any R object
#' @param what a string. If provided \code{sync_type(reg)} should be identical to \code{'what'}, otherwise \code{""} is returned.
#' @details sync name is the result of splitting \code{attr(reg, "sync")} with \code{"|"} and keeping all that is before the last \code{"|"} encountered.\cr
#' Note that \code{NA_character_} will be returned in case of error during splitting i.e. no \code{"|"} found or \code{length(attr(reg, "sync")) != 1}.
#' @return a string containing sync name
#' @keywords internal
sync_name <- function(reg, what) {
  A = attr(reg, "sync")
  pre = try(splits(A, name = TRUE, split = "|"), silent = TRUE)
  if(inherits(pre, "try-error")) return(NA_character_)
  if(missing(what)) return(pre)
  if(identical(what, sync_type(reg, what))) return(pre)
  return("")
}

#' @title Sync Parse
#' @description
#' Helper to serialize sync name.
#' @param regions an `IFC_regions` object, though it can be any R object
#' @param what a string
#' @return a character vector of sync name.
#' @keywords internal
sync_parse <- function(regions, what) {
  if(missing(what)) return(unlist(recursive = FALSE, sapply(regions, simplify = FALSE, USE.NAMES = TRUE, sync_name)))
  return(unlist(recursive = FALSE, sapply(regions, what, simplify = FALSE, USE.NAMES = TRUE, sync_name)))
}

#' @title Sync Check
#' @description
#' Helper to check sync attribute.
#' @param regions an `IFC_regions` object.
#' @return data.frame with \code{'name'} of invalid sync attribute and \code{'error'} description
#' @keywords internal
sync_check <- function(regions) {
  if(length(regions) == 0) return(data.frame(name = NULL, error = NULL, row.names = NULL, stringsAsFactors = FALSE))
  N = names(regions)
  sync_names <- sync_parse(regions)
  has_sync = sapply(regions, FUN = function(x) length(attr(x, "sync")) != 0)
  is_sync = sapply(sync_names, FUN = function(x) !is.na(x) && !identical(x, ""))
  to_rem = N[has_sync & !is_sync]
  to_rem = data.frame(name = to_rem, error = rep(1, length(to_rem)), row.names = NULL, stringsAsFactors = FALSE)
  
  # validate dual/quad, check
  # -for dual that there are 2 sync regions "|1", and  "|2", and that type is "line"
  # -for quad that there are 4 sync regions "|1", "|2", "|3" and "|4", and that type is "rect"
  # -all x[1] are identical and same for y[1]
  d = data.frame(row.names = NULL,
                 name = N[is_sync],
                 sync = sync_names[is_sync],
                 part = sapply(regions[is_sync], FUN = function(r) sync_part(r)),
                 type = sapply(regions[is_sync], FUN = function(r) r$type),
                 x = sapply(regions[is_sync], FUN = function(r) r$x[1]),
                 y = sapply(regions[is_sync], FUN = function(r) r$y[1]),
                 stringsAsFactors = FALSE)
  rbind(to_rem, do.call(what = rbind, by(d, d$sync, simplify = FALSE, FUN = function(dd) {
    sync_typ = substring(dd$sync[1], 1, 4)
    err = integer(0)
    msg = switch(sync_typ,
                 "dual" = {
                   if(!identical(unique(dd$type), "line")) err = c(err, 2)
                   if(!(length(dd$name) == 2)) err = c(err, 3)
                   if(!identical(sort(dd$part), c("1","2"))) err = c(err, 11)
                   if(!((length(unique(dd$x)) == 1))) err = c(err, 12)
                   err
                 },
                 "quad" = {
                   if(!identical(unique(dd$type), "rect")) err = c(err, 4)
                   if(!(length(dd$name) == 4)) err = c(err,3)
                   if(!identical(sort(dd$part), c("1","2","3","4"))) err = c(err, 13)
                   if(!((length(unique(dd$x)) == 1) && (length(unique(dd$y)) == 1))) err = c(err, 14)
                   err
                 },
                 {
                   # err = c(err, 5) # TODO add err on unsupported sync
                   err
                 }
    )
    if(length(msg) == 0) return(NULL)
    data.frame(name = dd$name, error = paste0(msg, collapse = "|"), row.names = NULL, stringsAsFactors = FALSE)
  })))
}

#' @title Sync Validate
#' @description
#' Helper to validate sync regions.
#' @param regions an `IFC_regions` object.
#' @param partial a bool whether number of regions and their coordinates match sync attribute. Default is \code{FALSE} to only validate type.
#' @param quietly a bool whether to silent warning when a sync attribute is not valid. Default is \code{TRUE}.
#' @param caller a string, name of the calling function. Default is \code{character(0)}.
#' @return an `IFC_regions` with sync validated attributes
#' @keywords internal
sync_validate <- function(regions, partial = FALSE, quietly = TRUE, caller = character(0)) {
  N = names(regions)
  to_rem = sync_check(regions)
  return(structure(names = N, class = class(regions), lapply(seq_along(regions), FUN = function(i) {
    if(!any(N[i] %in% to_rem$name)) return(regions[[i]])
    code = as.integer(strsplit(to_rem$error[to_rem$name == N[i]], "|", fixed = TRUE)[[1]])
    if(partial) code = code[code <= 10]
    if(length(code) == 0) return(regions[[i]])
    msg = structure(c(
      "'sync' attribute is malformed",        # 1
      "regions should share same type 'line'",# 2
      "'sync' name should not be duplicated", # 3
      "regions should share same type 'rect'",# 4
      "'sync' type is not supported",         # 5
      "should be composed of exactly 2 'sync' regions", # 11
      "1st x coordinate should be identical among 'sync' regions", # 12
      "should be composed of exactly 4 'sync' regions", # 13
      "1st x/y coordinates should be identical among 'sync' regions" # 14
    ), names = c("1","2","3","4","5","11","12","13","14"))
    msg = msg[as.character(code)]
    r = regions[[i]]
    if(!quietly) {
      sync_typ = sync_type(r)
      warning(paste(caller, "can't validate 'sync'[\"",sep = ": "),paste0(attr(r, "sync"), collapse="\",\""),"\"] for ", r$label, "\n\t-",
              paste0(paste(paste0("\"", sync_typ, "\"", collapse =""), msg), collapse="\n\t-"), call. = FALSE, immediate. = TRUE)
    }
    attr(r, "sync") <- NULL
    return(r)
  })))
}

#' @title Dual Coordinates
#' @description
#' Helper to construct coordinates for "dual" sync.
#' @param coords a list containing x and y numeric vectors of at least 2 elements each, if not 1st element will be repeated twice.
#' @param part a string referencing sync part, as returned by \code{sync_part()}.
#' @param xlim x-limits. Default is missing to use \code{sync_lims()}, otherwise it will be \code{sync_lims(xlim)}.
#' @param ylim y-limits. Default is missing to use \code{sync_lims()}, otherwise it will be \code{sync_lims(ylim)}.
#' @details 1st x should be the coordinate where parts 1 and 2 are separated.
#' 1st y should be the y-coordinate for part 1 and 2nd y the one for part 2.
#' Then min/max xlim will be distributed to each part:\cr
#' -1: x = x[1], min(xlim); y = y[1], y[1]\cr
#' -2: x = x[1], max(xlim); y = y[2], y[2]\cr
#' otherwise, x = x[1], x[2]; y = y[1], y[2].\cr
#' In addition, \code{-Inf} x values will be replaced by \code{min(xlim)} and \code{+Inf} x values by \code{max(xlim)}, same for y values and \code{'ylim'}.
#' @return coordinates for the corresponding \code{part}
#' @keywords internal
dual_coords <- function(coords, part, xlim, ylim) {
  if(missing(xlim)) {xlim = sync_lims()} else {xlim = sync_lims(xlim)}
  if(missing(ylim)) {ylim = sync_lims()} else {ylim = sync_lims(ylim)}
  coords$x = rep_len(coords$x, 2)
  coords$y = rep_len(coords$y, 2)
  coords$x[coords$x == -Inf] <- min(xlim)
  coords$x[coords$x == Inf] <- max(xlim)
  coords$y[coords$y == -Inf] <- min(ylim)
  coords$y[coords$y == Inf] <- max(ylim)
  ans = list(x = c(coords$x[1], switch(part, "1"=min(xlim), "2"=max(xlim), coords$x[2])),
             y = switch(part, "1"=rep(coords$y[1], 2),"2"=rep(coords$y[2], 2), c(coords$y[1],coords$y[2])))
  return(ans)
}

#' @title Quad Coordinates
#' @description
#' Helper to construct coordinates for "quad" sync.
#' @param coords a list containing x and y numeric vectors of at least 2 elements each, if not 1st element will be repeated twice.
#' @param part a string referencing sync part, as returned by \code{sync_part()}.
#' @param xlim x-limits. Default is missing to use \code{sync_lims()}, otherwise it will be \code{sync_lims(xlim)}.
#' @param ylim y-limits. Default is missing to use \code{sync_lims()}, otherwise it will be \code{sync_lims(ylim)}.
#' @details 1st x and 1st y should be the x/y coordinates where parts 1-4 are separated.
#' Then, min/max xlim min/max ylim will be distributed to each part:\cr
#' -1: x = x[1], min(xlim); y = y[1], max(ylim)\cr
#' -2: x = x[1], max(xlim); y = y[1], max(ylim)\cr
#' -3: x = x[1], max(xlim); y = y[1], min(ylim)\cr
#' -4: x = x[1], min(xlim); y = y[1], min(ylim)\cr
#' otherwise, x = x[1], x[2]; y = y[1], y[2].\cr
#' In addition, \code{-Inf} x values will be replaced by \code{min(xlim)} and \code{+Inf} x values by \code{max(xlim)}, same for y values and \code{'ylim'}.
#' @return coordinates for the corresponding \code{part}
#' @keywords internal
quad_coords <- function(coords, part, xlim, ylim) {
  if(missing(xlim)) {xlim = sync_lims()} else {xlim = sync_lims(xlim)}
  if(missing(ylim)) {ylim = sync_lims()} else {ylim = sync_lims(ylim)}
  coords$x = rep_len(coords$x, 2)
  coords$y = rep_len(coords$y, 2)
  coords$x[coords$x == -Inf] <- min(xlim)
  coords$x[coords$x == Inf] <- max(xlim)
  coords$y[coords$y == -Inf] <- min(ylim)
  coords$y[coords$y == Inf] <- max(ylim)
  ans = list(x = c(coords$x[1], switch(part, "1"=min(xlim), "2"=max(xlim), "3"=max(xlim),"4"=min(xlim), coords$x[2])),
             y = c(coords$y[1], switch(part, "1"=max(ylim), "2"=max(ylim), "3"=min(ylim),"4"=min(ylim), coords$y[2])))
  return(ans)
}

#' @title Sync Coordinates
#' @description
#' Helper to construct coordinates for sync.
#' @param type type of sync. Expecting \code{"dual"} or \code{"quad"}.
#' @param coords a list containing x and y numeric vectors of at least 2 elements each, if not 1st element will be repeated twice. 
#' @param part a string referencing sync part, as returned by \code{sync_part()}.
#' @param xlim x-limits. Default is missing to use \code{sync_lims()}, otherwise it will be \code{sync_lims(xlim)}.
#' @param ylim y-limits. Default is missing to use \code{sync_lims()}, otherwise it will be \code{sync_lims(ylim)}.
#' @details see \code{dual_coords()} or \code{quad_coords()}.
#' @return coordinates for the corresponding \code{part}
#' @keywords internal
sync_coords <- function(type, coords, part, xlim, ylim) {
  if(missing(xlim)) {xlim = sync_lims()} else {xlim = sync_lims(xlim)}
  if(missing(ylim)) {ylim = sync_lims()} else {ylim = sync_lims(ylim)}
  ans = switch(type,
               "dual" = dual_coords(coords, part, xlim, ylim),
               "quad" = quad_coords(coords, part, xlim, ylim),
               coords)
}

#' @title Create Sync
#' @description
#' Helper to construct sync regions.
#' @param reg a pre built region see \code{\link{buildRegion}} except that \code{reg$type} should be \code{"dual"} or \code{"quad"} instead of "line" or "rect" respectively.
#' @param xlim x-limits. Default is missing to use \code{sync_lims()}, otherwise it will be \code{sync_lims(xlim)}.
#' @param ylim y-limits. Default is missing to use \code{sync_lims()}, otherwise it will be \code{sync_lims(ylim)}.
#' @param forbidden a character vector of already given sync names, see \code{sync_name()}. Default is \code{character()}.
#' @details if \code{reg$type} is not \code{"dual"} or \code{"quad"}, it is returned as list(reg).
#' /!\ Note that when \code{'xlim'} is missing or is \code{c(-Inf,Inf)} label position \code{'reg$cx'} is hard to define so it falls back to \code{'reg$x'}. Same for \code{'ylim'} and \code{'reg$cy'} for \code{"quad"}.
#' @return a list of region(s)
#' @keywords internal
sync_create <- function(reg, xlim, ylim, forbidden = character()) {
  if(missing(xlim)) {xlim = sync_lims()} else {xlim = sync_lims(xlim)}
  if(missing(ylim)) {ylim = sync_lims()} else {ylim = sync_lims(ylim)}
  sync_typ = reg$type
  if(identical(sync_typ, "dual") || identical(sync_typ, "quad")) {
    reg_back = reg
    A = random_name(n = 10, prefix = paste0(sync_typ,"_"), special = FALSE, forbidden = forbidden)
    reg$type = ifelse(identical(reg$type, "dual"), "line", "rect")
    if(sync_typ == "dual") { reg$ylogrange = "P"; reg$ytrans=NULL }
    if(sync_typ == "dual") reg$y[2] = reg$y[1]
    reg$x = rep(reg$x, length.out = 2)
    reg$y = rep(reg$y, length.out = 2)
    reg = keep_attributes(reg, what=buildRegion)
    if(sync_typ == "dual") reg$y[2] = reg_back$y[2]
    reg[setdiff(names(reg_back), names(reg))] <- reg_back[setdiff(names(reg_back), names(reg))]
    n = switch(sync_typ, "dual" = 2, "quad" = 4, 0)
    lapply(seq_len(n), FUN = function(i) {
      attr(reg, "sync") = paste0(A, "|", i)
      reg$label = paste0(reg$label, "_", i)
      Xtrans=reg$xtrans; if(length(Xtrans) == 0) Xtrans=reg$xlogrange
      trans_x=parseTrans(Xtrans)
      Ytrans=reg$ytrans; if(length(Ytrans) == 0) Ytrans=reg$ylogrange
      trans_y=parseTrans(Ytrans)
      coords = list(x=reg$x,y=reg$y)
      if(n == 2) {
        if(i == 2) {
          if(length(unique(coords$y)) == 1) coords$y = ifelse(coords$y[1] > 0.5, coords$y[1] - 0.1, coords$y[1] + 0.1)
          coords$y = rep(coords$y, 2)
        }
        coords = dual_coords(coords, as.character(i), xlim, ylim)
      }
      if(n == 4) coords = quad_coords(coords, as.character(i), xlim, ylim)
      reg$x = coords$x
      reg$y = coords$y
      reg$cx = applyTrans(mean(applyTrans(sync_clip(reg$x), trans_x)), trans_x, inverse = TRUE) # TODO if xlim is missing or c(-Inf,+Inf) cx is hard to define
      reg$cy = applyTrans(mean(applyTrans(sync_clip(reg$y), trans_y)), trans_y, inverse = TRUE) # TODO if ylim is missing or c(-Inf,+Inf) cy is hard to define
      reg
    })
  } else {
    list(reg)
  }
}
