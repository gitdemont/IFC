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

#' @title Feature Values Computation
#' @name get_feat_value
#' @description
#' Computes feature values from feature definition
#' @param features a data.frame of features, typically an object of class `IFC_features`.
#' @param feat_def a feature definition as created by \code{\link{buildFeature}}.
#' @param operators operators used. Default is c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT").
#' @param split string used for splitting. Default is "|".
#' @details if 'feat_def$type' is "combined" returned value will be computed according in the scope of 'features' according to 'feat_def$def'.
#' Otherwise, 'feat_def$name' will be searched in 'features' to return value, and if not found 'feat_def$val' will be returned.
#' @return a vector of feature values.
#' @keywords internal
get_feat_value <- function(feat_def,
                           features,
                           operators = c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT"),
                           split = "|") {
  if(length(feat_def) == 0) return(NULL)
  # if feature is not combined value is either already in features or we take it from val
  if(feat_def$type != "combined") {
    if(feat_def$name %in% names(features)) return(features[, feat_def$name])
    return(feat_def$val)
  } 
  
  # identify features names and operators in feature definition
  def_tmp = splitn(definition = feat_def$def, all_names = names(features), operators = operators, scalar = TRUE)
  def_names = setdiff(def_tmp, operators)

  # variables used
  not_fun = setdiff(operators, c("+", "-", "*", "/", "(", ")"))
  alw_fun = sapply(setdiff(tolower(operators), c(")","sqr")), USE.NAMES = TRUE, simplify = FALSE,
                   FUN = function(x) getFromNamespace(x, asNamespace("base")))
  alw_fun = c(alw_fun, list("sqr" = function(x) x^2))
  
  # initialize bracket counters
  n = 0; w = 0; def_str = c()
  
  # add necessary bracket to string definition 
  for(i in seq_along(def_tmp)) {
    foo = def_tmp[i] %in% not_fun
    if(any(foo)) {
      def_str = c(def_str, tolower(def_tmp[i]))
    } else {
      def_str = c(def_str, def_tmp[i])
    }
    if(def_tmp[i] == "(") {
      w = w + 1
      next
    }
    if(def_tmp[i] == ")") {
      w = w - 1
      next
    }
    if(def_tmp[i] %in% not_fun) {
      n = n + 1
      def_str = c(def_str, "(")
    } else {
      if(w < n) {
        n = n - 1
        def_str = c(def_str, ")")
      }
    }
  }
  
  # terminate string definition correction with remaining brackets to close
  replicate(n, { def_str <<- c(def_str, ")") })
  
  # replace features names by their values and compute result according to corrected feature definition
  def_names=def_names[is.na(suppressWarnings(as.numeric(def_names)))]
  replace_with=c()
  for(i_def in seq_along(def_names)) replace_with=c(replace_with,random_name(n=10,special=NULL,forbidden=c(replace_with,def_str)))
  for(i_def in seq_along(def_names)) def_str[def_names[i_def]==def_str] <- rep(paste0("`",replace_with[i_def],"`"),sum(def_names[i_def]==def_str))
  e = lapply(def_names, FUN = function(x) features[ , x, drop = TRUE])
  names(e)=replace_with
  suppressWarnings(eval(expr=parse(text=paste0(def_str,collapse=" ")),envir=c(e, alw_fun),enclos=emptyenv()))
}

#' @title Features Values Extraction
#' @name getFeaturesValues
#' @description
#' Extracts features values according to features definitions
#' @param features a data.frame of features, typically an object of class `IFC_features`.
#' @param features_def a list of features definitions, typically an object of class `IFC_features_def`.
#' @param operators operators used. Default is c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT").
#' @param split string used for splitting. Default is "|".
#' @return a data.frame of features values.
#' @keywords internal
getFeaturesValues <- function(features,
                              features_def,
                              operators = c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT"),
                              split = "|", ...) {
  if(length(features_def) == 0) return(features)
  f_names = names(features)
  d_names = sapply(features_def, FUN = function(f_def) f_def$name)
  names(features_def) = d_names
  defs = lapply(features_def, FUN = function(f_def) {
    if(f_def$type != "combined") return(NULL)
    def_tmp = splitn(definition = f_def$def, all_names = c(f_names, d_names), operators = operators, scalar = TRUE)
    setdiff(def_tmp, operators)
  })
  names(defs) = d_names

  # order features
  i = 1
  l = length(defs)
  while (i < l) {
    index = defs[[i]]
    index = unlist(lapply(index, function(x) which(x == names(defs))))
    index = index[index > i]
    if(length(index) != 0) {
      defs = c(defs[index], defs[setdiff(1:l, index)])
      i = 1
    } else {
      i = i + 1
    }
  }
  
  # get features values
  for(i_name in unique(c(f_names, names(defs)))) {
    N = names(features)
    if(!any(i_name == N)) {
      # i_name does not exists yet in features
      # so we compute features values according to definition
      # and add it to features
      v = get_feat_value(features = features,
                         feat_def = features_def[[i_name]],
                         operators = operators,
                         split = split)
      if(length(v) == nrow(features)) {
        features = fastCbind(features,
                             structure(list(v), names = features_def[[i_name]]),
                             add_id = FALSE) 
      } else {
        stop("can't extract value for feature [",i_name,"]")
      }
    } else {                # recomputes features values
      if((i_name %in% d_names) && (features_def[[i_name]]$type == "combined")) {
        features[, i_name] = get_feat_value(features = features,
                                                     feat_def = features_def[[i_name]],
                                                     operators = operators,
                                                     split = split)
      }
    }
  }
  features
}
