#' @title Progress Bar Initializer
#' @description
#' Initializes a progress bar.
#' @param session the shiny session object, as provided by shinyServer to the server function. Default is missing, to use "txtProgressBar" or "winProgressBar" (on Windows).
#' @param title,label character strings, giving the 'title'(='message' for shiny progress bar) and the 'label'(='detail' for shiny progress bar).
#' @param min,max (finite) numeric values for the extremes of the progress bar. Must have 'min' < 'max'.
#' @param initial initial value for the progress bar.
#' @param width only apply when 'session' is missing,the width of the progress bar. If missing, the default, will be NA for "txtProgressBar" and 300 for "winProgressBar".
#' @param style does not apply for "winProgressBar", the ‘style’ of the bar. If missing, the default, will be 3 "txtProgressBar" and getShinyOption("progress.style", default = "notification") for shiny progress bar
#' @param char only apply for "txtProgressBar", the character (or character string) to form the progress bar.
#' @param file only apply for "txtProgressBar", an open connection object or "" which indicates the console: stderr() might be useful here. Default is "".
#' @return an object of class "txtProgressBar" "winProgressBar" or "Progress".
#' @keywords internal
newPB <- function(session = getDefaultReactiveDomain(),
                  title, label, 
                  min = 0, max = 1,
                  initial = 0,
                  width,
                  style,
                  char = "=",
                  file = "") {
  fun = stop
  args = list("newPB: can't create progress bar")
  if(length(session) == 0) {
    if(.Platform$OS.type == "windows") {
      args = list(min = min,
                  max = max,
                  initial = initial)
      if(missing(title)) {
        args = c(args, list(title = "R progress bar"))
      } else {
        args = c(args, list(title = title))
      }
      if(missing(label)) {
        args = c(args, list(label = ""))
      } else {
        args = c(args, list(label = label))
      }
      if(missing(width)) {
        args = c(args, list(width = 300))
      } else {
        args = c(args, list(width = width))
      }
      fun = winProgressBar
    } else {
      args = list(min = min,
                  max = max, 
                  initial = initial, 
                  char = char,
                  title = title, 
                  label = label, 
                  file = file)
      if(missing(style)) {
        args = c(args, list(style = 3))
      } else {
        args = c(args, list(style = style))
      }
      if(missing(width)) {
        args = c(args, list(width = NA))
      } else {
        args = c(args, list(width = width))
      }
      fun = txtProgressBar
    }
  } else {
    args = list(session = session,
                min = min,
                max = max)
    if(missing(style)) {
      args = c(args, list(style = getShinyOption("progress.style", default = "notification")))
    } else {
      if(style[1] %in% c("old", "notification")) {
        args = c(args, list(style = style))
      } else {
        args = c(args, list(style = getShinyOption("progress.style", default = "notification")))
      }
    }
    fun = Progress$new
  }
  return(do.call(what = fun, args = args))
}

#' @title Progress Bar Modifyer
#' @description
#' Changes a progress bar.
#' @param pb an object of class "txtProgressBar" "winProgressBar" or "Progress".
#' @param value new value for the progress bar.
#' @param title,label character strings, giving the 'title'(='message' for shiny progress bar) and the 'label'(='detail' for shiny progress bar).
#' @return an object of class "txtProgressBar" "winProgressBar" or "Progress".
#' @keywords internal
setPB <- function(pb, value = NULL, title = NULL, label = NULL) {
  K = class(pb)
  if("txtProgressBar" %in% K) utils::setTxtProgressBar(pb = pb, value = value, title = title, label = label)
  if("txtProgressBar" %in% K) utils::setWinProgressBar(pb = pb, value = value, title = title, label = label)
  if("Progress" %in% K) pb$set(value = value, message = title, detail = label)
}

#' @title Progress Bar Terminator
#' @description
#' Closes a progress bar.
#' @param pb an object of class "txtProgressBar" "winProgressBar" or "Progress".
#' @return an object of class "txtProgressBar" "winProgressBar" or "Progress".
#' @keywords internal
endPB <- function(pb) {
  K = class(pb)
  if("Progress" %in% class(pb)) {
    pb$close()
  } else {
    close(con = pb)
  }
}