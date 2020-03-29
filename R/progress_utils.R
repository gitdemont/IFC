#' @title Progress Bar Initializer
#' @description
#' Initializes a progress bar.
#' @param session the shiny session object, as provided by shinyServer to the server function. Default is missing, to use "txtProgressBar" or "winProgressBar" (on Windows).
#' @param title,label character strings, giving the 'title'(='message' for shiny progress bar) and the 'label'(='detail' for shiny progress bar).
#' @param min,max (finite) numeric values for the extremes of the progress bar. Must have 'min' < 'max'.
#' @param initial initial value for the progress bar.
#' @param steps (finite) numeric value for the number of individual chunk of the progress bar. Default is 21.
#' @param width only apply when 'session' is missing,the width of the progress bar. If missing, the default, will be NA for "txtProgressBar" and 300 for "winProgressBar".
#' @param style does not apply for "winProgressBar", the ‘style’ of the bar. If missing, the default, will be 3 "txtProgressBar" and getShinyOption("progress.style", default = "notification") for shiny progress bar
#' @param char only apply for "txtProgressBar", the character (or character string) to form the progress bar.
#' @param file only apply for "txtProgressBar", an open connection object or "" which indicates the console: stderr() might be useful here. Default is "".
#' @return an object of class "IFC_progress" containing a progress bar of class "txtProgressBar" "winProgressBar" or "Progress".
#' @keywords internal
newPB <- function(session,
                  title, label, 
                  min = 0, max = 1,
                  initial = 0,
                  steps = 21,
                  width,
                  style,
                  char = "=",
                  file = "") {
  fun = stop
  args = list("newPB: can't create progress bar")
  if(length(session) == 0) {
    if(.Platform$OS.type == "windows") {
      args = list(min = 0,
                  max = steps - 1,
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
      bar = do.call(what = fun, args = args)
      setWinProgressBar(bar, value = min)
    } else {
      args = list(min = 0,
                  max = steps - 1, 
                  initial = initial, 
                  char = char,
                  file = file)
      mess = c()
      if(!missing(title)) mess = title
      if(!missing(label)) mess = c(mess, label)
      if(length(mess) != 0 )cat("\n", mess, "\n", sep = "\n", file = file)
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
      bar = do.call(what = fun, args = args)
      setTxtProgressBar(bar, value = min)
    }
  } else {
    args = list(session = session,
                min = 0,
                max = steps- 1)
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
    bar = do.call(what = fun, args = args)
    bar$set(value = min)
  }
  ans = list(bar = bar,
             seq = seq(min, max, length.out = steps))
  class(ans) = "IFC_progress"
  return(ans)
}

#' @title Progress Bar Updater
#' @description
#' Updates a progress bar.
#' @param pb an object of class "IFC_progress" containing a progress bar of class "txtProgressBar" "winProgressBar" or "Progress".
#' @param value new value for the progress bar.
#' @param title,label character strings, giving the 'title'(='message' for shiny progress bar) and the 'label'(='detail' for shiny progress bar).
#' @keywords internal
setPB <- function(pb, value = NULL, title = NULL, label = NULL) {
  K = class(pb$bar)
  if("txtProgressBar" %in% K) {
    if(pb$bar$getVal() == which(value <= pb$seq)[1]) return(NULL)
    return(setTxtProgressBar(pb = pb$bar, value = which(value <= pb$seq)[1], title = title, label = label))
  }
  if("winProgressBar" %in% K) {
    if(getWinProgressBar(pb$bar) == which(value <= pb$seq)[1]) return(NULL)
    return(setWinProgressBar(pb = pb$bar, value = which(value <= pb$seq)[1], title = title, label = label))
  }
  if("Progress" %in% K) {
    if(pb$bar$getValue() == which(value <= pb$seq)[1]) return(NULL)
    return(pb$bar$set(value = which(value <= pb$seq)[1], message = title, detail = label))
  }
  return(NULL)
}

#' @title Progress Bar Terminator
#' @description
#' Terminates a progress bar.
#' @param pb an object of class "IFC_progress" containing a progress bar of class "txtProgressBar" "winProgressBar" or "Progress".
#' @keywords internal
endPB <- function(pb) {
  K = class(pb$bar)
  if("Progress" %in% class(pb$bar)) {
    pb$bar$close()
  } else {
    close(con = pb$bar)
  }
}