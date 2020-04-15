#' @title Object Transformation
#' @description
#' Function to normalize, colorize and add background to images.
#' @param mat a finite numeric matrix.
#' @param msk a finite numeric matrix (mask identifying abnormalities). If missing, the default no cleansing will be done.
#' @param color a color.
#' @param input_range a finite numeric vector of 2 values, sets the 'input_range' of the input intensity values; values exceeding this 'input_range' are clipped.
#' @param mode color mode export. Either "rgb", "gray" or "raw".
#' @param type image object type.
#' @param add_noise logical, if TRUE adds normal noise to background using rnorm(), from \pkg{Rcpp}. Default is TRUE.
#' @param random_seed a single value, interpreted as an integer, or NULL to be used with set.seed() from \pkg{base} when 'add_noise' is set to TRUE. Default is NULL.
#' @param size a length 2 integer vector of final dimensions of the image, height 1st and width 2nd. Default is c(0,0) for no change.
#' @param bg_mean mean value of the background added. Default is 0.
#' @param bg_sd standard deviation of the background added. Default is 0.
#' @param full_range logical, only apply when mode is not "raw", if 'full_range' is TRUE, then input_range will be set to c(0, 4095) and 'gamma' forced to 1. Default is FALSE.
#' @param force_range logical, only apply when mode is not "raw", if 'force_range' is TRUE, then input_range will be adjusted to object range in [-4095, +Inf] and 'gamma' forced to 1. Default is FALSE.\cr
#' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
#' @param gamma gamma correction. Default is 1, for no correction.
#' @details When 'add_noise' is FALSE and 'msk', `removal` attribute has value "masked" or "MC",
#' backgound will be automatically set to minimal pixel value.
#' @return the matrix transformed according to input parameters
#' @export
objectTransform <- function(mat, msk, color, input_range, mode, type, 
                            add_noise = TRUE, random_seed = NULL, size = c(0,0),
                            bg_mean = 0, bg_sd = 0, full_range = FALSE, force_range = FALSE, gamma = 1) {
  foo = mat
  bg_2 = bg_mean
  sd_2 = bg_sd
  if(!missing(msk)) {
    if((attr(msk, "removal") %in% c("masked","raw")) && (!add_noise)) {
      bg_2 = -4096
      sd_2 = 0
    }
    foo = objectCleanse(mat = foo, msk = msk, add_noise = add_noise, random_seed = random_seed, bg = bg_mean, sd = bg_sd)
  }
  foo = objectResize(mat = foo, size = size, add_noise = add_noise, random_seed = random_seed, bg = bg_mean, sd = bg_sd)
  switch(mode,
         "gray" = {
           foo = objectNormalize(mat = foo, input_range = input_range, full_range = full_range, force_range = force_range, gamma = gamma)
         },
         "rgb" = {
           foo = objectColorize(mat = objectNormalize(mat = foo, input_range = input_range, full_range = full_range, force_range = force_range, gamma = gamma), color)
         })
  attr(foo, "input_range") <- input_range
  attr(foo, "full_range") <- full_range
  attr(foo, "force_range") <- force_range
  attr(foo, "gamma") <- gamma
  attr(foo, "color") <- color
  attr(foo, "mode") <- mode
  attr(foo, "raw") <- mat
  attr(foo, "BG_MEAN") <- bg_mean
  attr(foo, "BG_STD") <- bg_sd
  if(type == 2) {
    attr(foo, "class") <- "IFC_img"
  } else {
    attr(foo, "class") <- "IFC_msk"
  }
  return(foo)
}
