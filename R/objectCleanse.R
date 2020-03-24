#' @title IFC_object Cleanser
#' @description
#' Removes abnormalities (clipped/debris) from image
#' @param mat a numeric matrix (image).
#' @param msk a numeric matrix (mask identifying abnormalities).
#' @param add_noise if TRUE adds normal noise to background using rnorm(), from \pkg{Rcpp}. Default is TRUE.
#' @param random_seed a single value, interpreted as an integer when add_noise is set to TRUE. Default is 1.
#' @param bg mean value of the background added if add_noise is TRUE. Default is 0.
#' @param sd standard deviation of the background added if add_noise is TRUE. Default is 0.
#' @return According to msk, pixel values in mat are substituted by either bg [add_noise == FALSE] or rnorm(n = prod(dim(mat), mean=bg, sd=sd)) [add_noise == TRUE].
#' @export
objectCleanse = function(mat, msk, add_noise = TRUE, random_seed = 1, bg = 0, sd = 0) {
  if(add_noise) {
    set.seed(random_seed)
    on.exit(set.seed(NULL))
  }
  foo = cpp_cleanse(mat, msk, add_noise = add_noise, bg = bg, sd = sd)
  attr(foo, "msk_cleanse") <- msk
  return(foo)
}
