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

#' @title Add Feature to IFC_data Object
#' @description
#' Adds features to an already existing `IFC_data` object.
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param features a list of features to add to obj. Each element of this list will be coerced by \code{\link{buildFeature}}.
#' @details A warning will be thrown if a provided feature is already existing in obj.\cr
#' In such a case this feature will not be added to obj.\cr
#' If any input feature is not well defined and can't be created then an error will occur.
#' @param ... Other arguments to be passed.
#' @examples
#' if(requireNamespace("IFCdata", quietly = TRUE)) {
#'   ## use a daf file
#'   file_daf <- system.file("extdata", "example.daf", package = "IFCdata")
#'   daf <- ExtractFromDAF(fileName = file_daf)
#'   ## copy 1st feature found in daf
#'   feat_def <- daf$features_def[[1]]
#'   if(length(feat_def) != 0) {
#'     feat_def_copy <- feat_def
#'     ## modify name and value of copied features
#'     feat_def_copy$name <- "copied_feature"
#'     feat <- daf$features[, feat_def$name]
#'     feat_copy <- feat
#'     feat_copy <- feat_copy * 10
#'     ## create new object with this new feature
#'     dafnew <- data_add_features(obj = daf, features = list(c(feat_def_copy, list(val = feat_copy))))
#'   }
#' } else {
#'   message(sprintf('Please run `install.packages("IFCdata", repos = "%s", type = "source")` %s',
#'                   'https://gitdemont.github.io/IFCdata/',
#'                   'to install extra files required to run this example.'))
#' }
#' @return an IFC_data object with features added.
#' @export
data_add_features <- function(obj, features, ...) {
  assert(obj, cla = "IFC_data")
  
  # try to coerce features
  features = lapply(features, FUN=function(x) do.call(what=buildFeature, args=x))
  names(features) = sapply(features, FUN=function(x) x$name)
  
  # removes duplicated inputs
  tmp = duplicated(names(features))
  if(any(tmp)) {
    warning(paste0("duplicated features automatically removed: ", names(features)[tmp]), immediate. = TRUE, call. = FALSE)
    features = features[!tmp]
  }
  
  # defines available parameters
  obj_count = as.integer(obj$description$ID$objcount)
  
  # check that new features are not duplicated and well defined
  # send warning on duplicated
  # send error on bad definition
  all_names_comb = c(names(obj$features), names(features))
  alt_names_comb = gen_altnames(all_names_comb)
  all_msk = strsplit(obj$description$masks$def[obj$description$masks$name == "MC"], split = "|Or|", fixed = TRUE)[[1]]
  all_names_avl = unname(c("true", "false", "True", "False",
                           all_msk, obj$description$masks$name,
                           obj$description$Images$name,
                           unlist(featureIFC(), recursive = FALSE, use.names = FALSE)))
  alt_names_avl = gen_altnames(all_names_avl)
  
  exported_feats = sapply(features, FUN=function(feat) {
    if(feat$name%in%names(obj$features)) {
      warning(paste0(feat$name, "\nnot exported: trying to export an already defined feature"), immediate. = TRUE, call. = FALSE)
      return(FALSE)
    }
    if(feat$type == "combined") {
      splitn(definition = feat$def,
             all_names = all_names_comb,
             alt_names = alt_names_comb,
             operators = c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT"),
             scalar = TRUE) # will error if definition is not ok
    } else {
      splitn(definition = feat$def,
             all_names = all_names_avl,
             alt_names = alt_names_avl,
             operators = character(), 
             scalar = TRUE, dsplit = TRUE) # will error if definition is not ok
    }
    return(TRUE)
  })
  exported_feats = features[exported_feats] # remove duplicated
  if(length(exported_feats) == 0) return(obj)
  
  # append new features to obj$features
  names(exported_feats) = sapply(exported_feats, FUN=function(x) x$name)
  obj$features = fastCbind(obj$features, sapply(exported_feats, simplify = FALSE, FUN=function(x) {
    if(x$type == "combined") { # new "combined" features are initially filled with 0 
      return(rep(0, obj_count))
    } else {
      if(length(x$val) != obj_count) stop(x$name, "\nbad feature value length, expected: ",  obj_count, ", but is: ", length(x$val))
      return(x$val)
    }
  }))
  
  # record features names
  features_names = names(obj$features)
  
  # append new definition to obj$features_def
  K = class(obj$features_def)
  obj$features_def = c(obj$features_def, sapply(exported_feats, USE.NAMES = TRUE, simplify = FALSE,
                                                FUN=function(x) x[c("name", "type", "userfeaturetype", "def")]))
  class(obj$features_def) = c(setdiff(K, "IFC_features_def"), "IFC_features_def")
  
  # recompute new combined feature(s) and order obj$features
  K = class(obj$features)
  obj$features = getFeaturesValues(features_def = exported_feats[sapply(exported_feats, FUN = function(f_def) f_def$type == "combined")],
                                   features = obj$features)[, features_names]
  class(obj$features) = c(setdiff(K, "IFC_features"), "IFC_features")
  return(obj)
}
