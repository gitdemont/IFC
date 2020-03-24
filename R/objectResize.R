#' @title IFC_object Resizing
#' @description
#' Resizes mat to new dimensions.
#' @param mat a numeric matrix.
#' @param size a vector of final dimensions (Height, Width) of the desired returned mat. Default is NULL for no change.
#' if 'size' is of length one then only Width is changed.
#' @param add_noise if TRUE adds normal noise when size is larger than mat dimensions using rnorm(), from \pkg{Rcpp}. Default is TRUE.
#' @param random_seed a single value, interpreted as an integer when add_noise is set to TRUE. Default is 1.
#' @param bg mean value of the background added if add_noise is TRUE. Default is 0.
#' @param sd standard deviation of the background added if add_noise is TRUE. Default is 0.
#' @return a resized matrix with padding background if desired size is larger than original mat dimensions.
#' @export
objectResize <- function(mat, size = NULL, add_noise = TRUE, random_seed = 1, bg = 0, sd = 0) {
  if(length(size) == 0) return(mat)
  if(add_noise) {
    set.seed(random_seed)
    on.exit(set.seed(NULL))
  }
  if(length(size) == 1) {
    return(cpp_resize(mat = mat, new_height = 0, new_width = size, add_noise = add_noise, bg = bg, sd= sd))
  } else {
    return(cpp_resize(mat = mat, new_height = size[1], new_width = size[2], add_noise = add_noise, bg = bg, sd= sd))
  }
}
