#' @title Object File Export
#' @description
#' Exports images to various types.
#' @param x a numeric matrix.
#' @param type image type. Supported values are: "bmp", "jpeg", "png", and "tiff".
#' @param ... other arguments to be passed.
#' @keywords internal
objectWrite <- function(x, type, ...) {
  dots=list(...)
  switch(type, 
         png = png::writePNG(image = x, ...),
         tiff = tiff::writeTIFF(what = x, ...),
         jpeg = jpeg::writeJPEG(image = x, ...),
         bmp = {
           if((length(dots) != 0)) {
             switch(typeof(dots[[1]]),
                    raw = return(cpp_writeBMP(image = x)),
                    character = {
                      towrite = file(description = dots[[1]], open = "wb")
                      on.exit(close(towrite))
                      return(writeBin(object = cpp_writeBMP(image = x), con = towrite))
                    },
                    stop("objectWrite: target should be raw() or a character file path")
                    )
           }
           stop("objectWrite: no target to write \"bmp\" to")
         },
         stop(paste0("objectWrite: can only write to \"bmp\", \"jpeg\", \"png\" and \"tiff\"; [",type,"] is not supported")))
}
