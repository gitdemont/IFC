#' @title Smooth LinLog Transformation
#' @description
#' Transforms values in lin-log
#' @param x A numeric vector.
#' @param hyper value where transition between Lin/Log is applied.
#' @param base base of Log scale.
#' @param lin_comp value that is used to smooth transition between Lin/Log. Default is log(base).
#' @export
smoothLinLog <- function(x, hyper=1000, base=10, lin_comp=log(base)) {
  stopifnot(hyper > 0, base > 0, lin_comp >0)
  return(cpp_smoothLinLog(x, hyper, base, lin_comp))
}
