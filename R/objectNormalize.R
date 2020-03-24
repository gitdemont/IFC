#' @title IFC_object Intensity Normalization
#' @description
#' Normalizes a matrix to [0,1].
#' @param mat a finite numeric matrix.
#' @param input_range a finite numeric vector of 2 values, sets the range of the input intensity values.\cr
#' Values exceeding this range are clipped. Default is c(0, 4095).
#' @param full_range if 'full_range' is TRUE, then 'input_range' will be set to c(0, 4095) and 'gamma' forced to 1. Default is FALSE.
#' @param force_range if 'force_range' is TRUE, then 'input_range' will be adjusted to mat range and 'gamma' forced to 1. Default is FALSE.\cr
#' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
#' @param gamma gamma correction. Default is 1, for no correction.
#' @details Note that negative values are used internally for removal of unmasked objects.
#' @return a [0,1] normalized matrix
#' @export
objectNormalize <- function(mat, input_range=c(0,4095), full_range=FALSE, force_range=FALSE, gamma=1) {
  cpp_normalize(mat = mat, input_range = input_range, full_range = full_range, force_range = force_range, gamma = gamma)
}
