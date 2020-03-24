#' @title IFC Feature Coercion
#' @description
#' Helper to build a list to allow feature export.
#' @param name feature's name. If missing, it will be determined thanks to def.
#' @param type feature's type. Default is "single". Allowed are "single", "combined", "computed".
#' @param def definition of the feature. Default is "Camera Line Number".
#' @param val a coercible to numeric vector of feature values. Default is NULL.\cr
#' Note that although not mandatory for \code{\link{buildFeature}}it has to be provided to
#' allow feature export in \code{\link{ExportToDAF}} and \code{\link{data_add_features}}.
#' @param ... Other arguments to be passed.
#' @return a list containing all feature information.
#' @export
buildFeature <- function(name, type = c("single","combined","computed")[1], def = "Camera Line Number", val = NULL, ...) {
  dots = list(...)
  assert(def, len = 1, typ = "character")
  val = as.numeric(unlist(val))
  if(!missing(name)) assert(name, len = 1, typ = "character")
  assert(type, len = 1, alw = c("single","combined","computed"))
  
  feat_comp = strsplit(def, split="|", fixed = TRUE)[[1]]
  if(type == "single") {
    if(!(typeof(val) == "double" | typeof(val) == "integer")) stop("'val' should be coercible to numeric vector")
    feat = feat_comp[1]
    userfeaturetype_avl = list("Mask Only"=c("Area", "Aspect Ratio", "Length", "Width", "Height", "Angle", "Centroid X", "Centroid Y", "Circularity", "Diameter",
                                             "Elongatedness", "Major Axis", "Minor Axis", "Perimeter", "Shape Ratio", 'Spot Area Min', "Spot Count", "Spot Distance Min",
                                             "Thickness Max", "Thickness Min"),
                               "Mask and Image"=c("Aspect Ratio Intensity", "Modulation", "Contrast", "Gradient RMS", "Intensity", "Mean Pixel", "Median Pixel",
                                                  "Max Pixel", "Raw Max Pixel", "Raw Min Pixel", "Saturation Count", "Saturation Percent", "Bright Detail Intensity R3",
                                                  "Bright Detail Intensity R7", "Angle Intensity", "Centroid X Intensity", "Centroid Y Intensity", "Compactness","Gradient Max", "Internalization", 
                                                  "Lobe Count","Major Axis Intensity", "Max Contour Position", "Min Pixel", "Minor Axis Intensity", "Raw Intensity", "Raw Mean Pixel",
                                                  "Raw Median Pixel", "Spot Intensity Max", "Spot Intensity Min", "Std Dev", "Symmetry 2", "Symmetry 3", "Symmetry 4", 
                                                  "Uncompensated Intensity", "Valley X", "Valley Y"),
                               "Image Only"=paste0("Bkgd ", c("Mean", "Std")),
                               "No Parameters"=c("Time", "Object Number", "Raw Centroid X", "Raw Centroid Y", "Flow Speed", "Camera Line Number", "Camera Timer", "Objects per mL", "Objects per sec"),
                               "Mask, Image and Scalar"=paste0(paste0("H ", rep(c("Contrast ","Correlation ","Energy ", "Entropy ", "Homogeneity ", "Variance "), each=2)), c("Mean", "Std")),
                               "Mask and Scalar"=c("Spot Count"),
                               "Mask and Three Images"=c("Bright Detail Colocalization 3"),
                               "Similarity"=c("Bright Detail Similarity R3", "Shift X", "Shift Y", "Similarity", "XCorr"),
                               "Delta Centroid"=paste0("Delta Centroid ", c("X","Y","XY")),
                               "Two Masks and Image"=c("Intensity Concentration Ratio"),
                               "Image and Scalar"=c("Ensquared Energy", "Diameter:"))
    type_avl = names(userfeaturetype_avl)
    feat_type = type_avl[sapply(type_avl, FUN=function(i) {any(feat==userfeaturetype_avl[[i]])})]
    if(feat == "Spot Count") {
      isScalar = feat_comp[length(feat_comp)-2]
      feat_type = "Mask Only"
      if(length(isScalar)!=0) if(isScalar == "True" | isScalar == "False") feat_type = "Mask and Scalar"
    }
  } else {
    feat_type = "Combined"
  }
  if(length(feat_type)==0) stop(paste0("Oops... 'def' is not correct, allowed are: ",paste(unlist(userfeaturetype_avl), collapse = ", "), "\n"))
  if(missing(name)) name = paste0(feat_comp, collapse = "_")
  return(list("name"=name, "type"=type, "userfeaturetype"=feat_type, "def"=def, "val"=val))
}
