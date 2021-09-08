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

#' @title `IFC_data` Object Check 
#' @description 
#' Ensures `IFC_data` obj$features correctly reflects actual number of collected objects.
#' When 'obj' is from rif or cif, number of objects can be different from actual number of collected objects
#' e.g. when it comes from merged of subset of file. There is no way to link an object with its feature value,
#' so the only solution is to remove all features and only keep "Object Number".
#' @param obj an `IFC_data` object.
#' @return an `IFC_data` object.
#' @keywords internal
checkObj <- function(obj) {
  assert(obj,cla="IFC_data")
  obj_count=as.integer(obj$description$ID$objcount)
  if((nrow(obj$features) != 0) && (obj_count == nrow(obj$features))) {
    if(!("Object Number" %in% names(obj$features))) {
      if(ncol(obj$features) == 0) {
        obj$features=data.frame("Object Number"= 0:(obj_count-1),check.names=FALSE)
        obj$features_def=structure(list(list(name="Object Number",type="single",userfeaturetype="No Parameters",def="Object Number")))
      } else {
        obj$features=cbind(obj$features,"Object Number"=0:(obj_count-1))
        obj$features_def=c(obj$features_def,structure(list(list(name="Object Number",type="single",userfeaturetype="No Parameters",def="Object Number"))))
      }
      class(obj$features)=c("data.frame","IFC_features")
      class(obj$features_def)=c("list","IFC_features_def")
    }
  } else {
    warning("Actual number of collected objects differs from features recorded within 'obj'", call. = FALSE, immediate. = TRUE)
    obj$features=structure(data.frame("Object Number"=0:(obj_count-1),check.names = FALSE), 
                           class = c("data.frame","IFC_features"))
    obj$features_def=structure(list(buildFeature(name="Object Number",val=0:(obj_count-1),def="Object Number")[1:4]),
                               names = "Object Number",class = c("list","IFC_features_def"))
    obj$pops=structure(list(list(name="All",type="B",base="",color="White",lightModeColor = "Black",
                                 style = 20,names="",obj=rep(TRUE, obj_count))),
                       names="All",class=c("Affiliated","Ordered","IFC_pops","Processed"))
    obj$graphs=structure(list(),class="IFC_graphs")
    obj$regions=structure(list(),class="IFC_regions")
  }
  # we ensure features are in a nice order
  if(ncol(obj$features) != 1) {
    obj$features=structure(obj$features[, order(names(obj$features))],class=c("data.frame","IFC_features"))
    obj$features_def=structure(obj$features_def[order(names(obj$features_def))],class=c("list","IFC_features_def"))
  }
  attr(obj$pops[["All"]], "reserved")=TRUE
  stats=data.frame(stringsAsFactors=FALSE,check.rows=FALSE,check.names=FALSE,t(sapply(names(obj$pops),FUN=function(p){
    count=sum(obj$pops[[p]]$obj)
    base=obj$pops[[p]]$base
    type=obj$pops[[p]]$type
    if(base == "") base="All"
    parent=sum(obj$pops[[base]]$obj)
    c("type"=type,"parent"=base,"count"=count,"perc_parent"=count/parent*100,"perc_tot"=count/obj_count*100)
  })))
  stats[,3]=as.numeric(stats[,3])
  stats[,4]=as.numeric(stats[,4])
  stats[,5]=as.numeric(stats[,5])
  obj$stats=stats
  return(obj)
}