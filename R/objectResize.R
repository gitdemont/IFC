#' @title Object Resizing
#' @description
#' Resizes mat to new dimensions.
#' @param mat a numeric matrix.
#' @param size a length 2 integer vector of final dimensions of the image, height 1st and width 2nd. Default is c(0,0) for no change.
#' @param add_noise if TRUE adds normal noise when size is larger than mat dimensions using rnorm(), from \pkg{Rcpp}. Default is TRUE.
#' @param random_seed a single value, interpreted as an integer, or NULL to be used with set.seed() from \pkg{base} when 'add_noise' is set to TRUE. Default is NULL.
#' @param bg mean value of the background added if add_noise is TRUE. Default is 0.
#' @param sd standard deviation of the background added if add_noise is TRUE. Default is 0.
#' @return a resized matrix with padding background if desired size is larger than original mat dimensions.
#' @export
objectResize <- function(mat, size = c(0,0), add_noise = TRUE, random_seed = NULL, bg = 0, sd = 0) {
  if(length(size) != 2) return(mat)
  if(add_noise) {
    set.seed(random_seed)
    on.exit(set.seed(NULL))
  }
  return(cpp_resize2(mat = mat, new_height = size[1], new_width = size[2], add_noise = add_noise, bg = bg, sd= sd))
}
