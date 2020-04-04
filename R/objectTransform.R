#' @title Object Transformation
#' @description
#' Function to normalize, colorize and add background to images.
#' @param mat a finite numeric matrix.
#' @param msk a finite numeric matrix (mask identifying abnormalities).
#' @param color a color.
#' @param input_range a finite numeric vector of 2 values, sets the 'input_range' of the input intensity values; values exceeding this 'input_range' are clipped.
#' @param mode color mode export. Either "rgb", "gray" or "raw".
#' @param type image object type.
#' @param cleanse logical, whether to cleanse the image or not. Default is FALSE.
#' @param add_noise logical, if TRUE adds normal noise to background using rnorm(), from \pkg{Rcpp}. Default is TRUE.
#' @param random_seed a single value, interpreted as an integer, or NULL to be used with set.seed() from \pkg{base} when 'add_noise' is set to TRUE. Default is NULL.
#' @param size a length 2 integer vector of final dimensions of the image. Default is NULL for no change.\cr
#' if 'size' is of length one then only width is changed.
#' @param bg_mean mean value of the background added. Default is 0.
#' @param bg_sd standard deviation of the background added. Default is 0.
#' @param full_range logical, only apply when mode is not "raw", if 'full_range' is TRUE, then input_range will be set to c(0, 4095) and 'gamma' forced to 1. Default is FALSE.
#' @param force_range logical, only apply when mode is not "raw", if 'force_range' is TRUE, then input_range will be adjusted to mat range and 'gamma' forced to 1. Default is FALSE.\cr
#' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
#' @param gamma gamma correction. Default is 1, for no correction.
#' @return the matrix transformed according to input parameters
#' @export
objectTransform <- function(mat, msk, color, input_range, mode, type, cleanse = FALSE,
                            add_noise = TRUE, random_seed = NULL, size = NULL,
                            bg_mean = 0, bg_sd = 0, full_range = FALSE, force_range = FALSE, gamma = 1) {
  foo = mat
  if(cleanse) foo = objectCleanse(mat = foo, msk = msk, add_noise = add_noise, random_seed = random_seed, bg = bg_mean, sd = bg_sd)
  if(!is.null(size)) foo = objectResize(mat = foo, size = size, add_noise = add_noise, random_seed = random_seed, bg = bg_mean, sd = bg_sd)
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
  attr(foo, "RAW") <- mat
  attr(foo, "BG_MEAN") <- ifelse(bg_mean >= 0, bg_mean, 0) # negative values (i.e. -1) are used for removal of non masked objects
  attr(foo, "BG_STD") <- bg_sd
  if(type == 2) {
    attr(foo, "class") <- "IFC_img"
  } else {
    attr(foo, "class") <- "IFC_msk"
  }
  return(foo)
}
