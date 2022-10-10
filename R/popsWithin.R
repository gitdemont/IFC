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

#' @title IFC_pops Checker for Object Membership to Populations
#' @description
#' Helper that will be used by popsCompute to determine which objects are within populations or not.
#' @param pops list of populations.
#' @param regions list of regions.
#' @param features dataframe of features.
#' @param pnt_in_poly_algorithm algorithm used to determine if object belongs to a polygon region or not. Default is 1.\cr
#' Note that for the moment only 1(Trigonometry) is available.
#' @param pnt_in_poly_epsilon epsilon to determine if object belongs to a polygon region or not. It only applies when algorithm is 1. Default is 1e-12.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param title_progress character string, giving the title of the progress bar. Default is "".
#' @param ... other arguments to be passed.
#' @source For pnt_in_poly_algorithm, Trigonometry, is an adaptation of Jeremy VanDerWal's code \url{https://github.com/jjvanderwal/SDMTools}
#' @keywords internal
popsWithin <- function(pops, regions, features, pnt_in_poly_algorithm = 1, pnt_in_poly_epsilon = 1e-12,
                       display_progress = TRUE, title_progress = "", ...) {
  dots = list(...)
  # several checks
  assert(pops, cla = c("IFC_pops","Affiliated","Ordered"))
  assert(regions, cla = "IFC_regions")
  assert(features, cla = "IFC_features")
  pnt_in_poly_algorithm = as.integer(pnt_in_poly_algorithm); assert(pnt_in_poly_algorithm, len = 1, alw = 1)
  pnt_in_poly_epsilon = as.numeric(pnt_in_poly_epsilon); pnt_in_poly_epsilon = pnt_in_poly_epsilon[pnt_in_poly_epsilon>0]; pnt_in_poly_epsilon = pnt_in_poly_epsilon[is.finite(pnt_in_poly_epsilon)]
  assert(pnt_in_poly_epsilon, len = 1, typ = "numeric")
  display_progress = as.logical(display_progress); assert(display_progress, len = 1, alw = c(TRUE, FALSE))
  assert(title_progress, len = 1, typ = "character")
  
  # variables used
  K = class(pops)
  L = length(pops)
  alw_fun = sapply(c("&","|","!","("), USE.NAMES = TRUE, simplify = FALSE,
                   FUN = function(x) getFromNamespace(x, asNamespace("base")))
  obj_number = nrow(features)
  if(display_progress) {
    pb = newPB(session = dots$session, min = 0, max = L, initial = 0, style = 3)
    on.exit(endPB(pb))
  }
  for(i in 1:L) {
    fx_pos = NULL
    fy_pos = NULL
    pop=pops[[i]]
    # changes styles to R compatible
    pops[[i]]$style = map_style(pops[[i]]$style, toR=TRUE)
    # changes colors to R compatible
    pops[[i]]$color = map_color(pops[[i]]$color)
    pops[[i]]$lightModeColor = map_color(pops[[i]]$lightModeColor)
    switch(pop$type,
           "B" = { 
             pops[[i]]$obj=rep(TRUE,obj_number)
           }, 
           "G" = {
             pop_pos=which(names(regions)==pop$region) # here there should be only one !
             if(length(pop_pos)!=1) stop(pop$name, ', trying to compute a graphical population with a non-defined region: ["', pop$region, '"]', call. = FALSE)
             fx_pos=which(names(features)==pop$fx)
             if(length(fx_pos)!=1) stop(pop$name, ', trying to compute a graphical population with an unknown fx ["', pop$fx, '"]', call. = FALSE)
             x=features[,fx_pos]
             xlim=as.numeric(regions[[pop_pos]]$x)
             if(regions[[pop_pos]]$type == "line") {
               xlim=range(xlim)
               pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & x>=xlim[1] & x<=xlim[2] & !is.na(x)
             } else {
               fy_pos=which(names(features)==pop$fy)
               if(length(fy_pos)!=1) stop(pop$name, ', trying to compute a graphical population with an unknown fy ["', pop$fy, '"]', call. = FALSE)
               y=features[,fy_pos]
               ylim=as.numeric(regions[[pop_pos]]$y)
               Xtrans = regions[[pop_pos]]$xtrans; if(length(Xtrans) == 0) Xtrans = regions[[pop_pos]]$xlogrange
               trans_x = parseTrans(Xtrans)
               x = applyTrans(x, trans_x)
               xlim = applyTrans(xlim, trans_x)
               Ytrans = regions[[pop_pos]]$ytrans; if(length(Ytrans) == 0) Ytrans = regions[[pop_pos]]$ylogrange
               trans_y = parseTrans(Ytrans)
               y = applyTrans(y, trans_y)
               ylim = applyTrans(ylim, trans_y)
               switch(regions[[pop_pos]]$type, 
                      "oval" = {
                        pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & cpp_pnt_in_gate(pnts=cbind(x,y), gate = cbind(xlim,ylim), algorithm = 3)
                      },
                      "poly" = {
                        pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & cpp_pnt_in_gate(pnts=cbind(x,y), gate = cbind(xlim,ylim), algorithm = pnt_in_poly_algorithm, epsilon = pnt_in_poly_epsilon)
                      },
                      "rect" = {
                        pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & cpp_pnt_in_gate(pnts=cbind(x,y), gate = cbind(xlim,ylim), algorithm = 2)
                      })
             }
           }, 
           "C" = {
             pop_def_tmp=pop$split
             pop_def_tmp[pop_def_tmp=="And"] <- "&"
             pop_def_tmp[pop_def_tmp=="Or"] <- "|"
             pop_def_tmp[pop_def_tmp=="Not"] <- "!"
             replace_with=gen_altnames(pop$names,forbidden=c(pop_def_tmp))
             for(i_def in seq_along(pop$names)) pop_def_tmp[pop$names[i_def] == pop_def_tmp] <- rep(paste0("`",replace_with[i_def],"`"), sum(pop$names[i_def] == pop_def_tmp))
             e = lapply(pops[pop$names], FUN=function(i_pop) i_pop$obj)
             names(e) = replace_with
             pops[[i]]$obj=pops[[which(names(pops)==pop$base)]]$obj & eval(expr=parse(text=paste0(pop_def_tmp,collapse="")),envir=c(e, alw_fun),enclos=emptyenv())
           }, 
           "T" = {
             if(length(pop$obj) != obj_number) {
               if(anyNA(pop$obj)) stop(pop$name, ", trying to compute a tagged population containing NA/NaN")
               Kp = typeof(pop$obj)
               if(Kp%in%c("double","integer")) {
                 if((obj_number <= max(pop$obj)) ||
                    (min(pop$obj) < 0) ||
                    any(duplicated(pop$obj))) stop(pop$name, ", trying to compute a tagged population with element(s) outside of objects acquired")
                 pops[[i]]$obj=rep(FALSE,obj_number)
                 pops[[i]]$obj[pop$obj+1]=TRUE
               } else {
                 if(!Kp%in%"logical") stop(pop$name, ', trying to compute a tagged population of unknown type ["',Kp,'"]')
               }
             }
             if(sum(pops[[i]]$obj)==0) stop(pop$name, ", trying to compute a tagged population of length = 0")
             if(obj_number != length(pops[[i]]$obj)) stop(pop$name, ", trying to compute a tagged population with more element(s) than total number of objects acquired")
           })
    if(display_progress) {
      setPB(pb, value = i, title = title_progress, label = "extacting populations")
    }
  }
  class(pops) = c(setdiff(K, "IFC_pops"), "IFC_pops", "Processed")
  return(pops)
}
