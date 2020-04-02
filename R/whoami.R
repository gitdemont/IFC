#' @title Input Identification
#' @description
#' Helper that identifies input arguments thanks to their IFC classes even if they are not named.
#' @param entries arguments from the function \code{\link{whoami}} is called.
#' @return a list whose members are:
#' -fileName: value of fileName if provided as a named argument,\cr
#' -info; value of `IFC_info` object identified,\cr
#' -param: value of `IFC_param` object identified,\cr
#' -offsets: value of `IFC_offset` object identified,\cr
#' -ifd: value of `IFC_ifd` object identified,\cr
#' -data: value of `IFC_data` object identified.\cr
#' @keywords internal
whoami = function(entries) {
  entries1 = entries[[1]]
  if(typeof(entries1) == "closure") {
    from = "closure"
  } else {
    from = as.character(entries1)
  }
  args = entries[-1]
  L = length(args)
  
  # retrieve arguments values
  val = lapply(1:L, FUN=function(i) {
    switch(typeof(args[[i]]),
           "symbol" = return(get(as.character(args[[i]]))),
           "language" = return(eval(args[[i]]))
    )
    return(args[[i]])
  })
  
  # identify arguments of classes we look for
  # classes = list("fileName",
  #                c("IFC_info", "IFC_display"), # to use former IFC_display 
  #                "IFC_param", 
  #                c("IFC_offset", "IFC_offsets"), # to use former IFC_offsets
  #                "IFC_ifd", 
  #                "IFC_data")
  # found = sapply(1:L, FUN = function(i) {
  #   sapply(classes, FUN = function(k) inherits(x = val[[i]], what = k, which = FALSE))
  # })
  classes = c("fileName", "IFC_info", "IFC_param", "IFC_offset", "IFC_ifd", "IFC_data")
  found = sapply(1:L, FUN = function(i) {
    inherits(x = val[[i]], what = classes, which = TRUE)
  })
  
  # count how many times a classes is present
  # error if at least 2 arguments are with same class
  times = apply(found, 1, sum)
  foo = times > 1
  if(any(foo)) {
    N = names(args)
    NN = sapply(1:L, FUN=function(i) {
      return(ifelse(N[i] == "", "<unk>", N[i]))
    })
    dup = which(foo)
    msg = sapply(1:length(dup), FUN = function(i) {
      idx = which(found[dup[i], ]==1)
      nam = NN[idx]
      paste0("too many elements of class `", classes[dup[i]], "`: [", paste0(paste(idx, nam, sep = "="), collapse = ",") ,"]")
    })
    stop(paste0(msg, collapse = "\n"))
  }
  
  # identify who is who and what it was
  iam = lapply(1:length(classes), FUN=function(i) {
    return(which(found[i, ] == 1))
  })
  was = rep(0, times = length(classes))
  was[times == 1] <- unlist(iam)
  new = lapply(iam, FUN = function(i) {
    if(length(i) != 0) return(val[[i]])
  })
  names(new) = c("fileName", "info", "param", "offsets", "ifd", "data")
  
  # add fileName if present
  fil = names(args) %in% "fileName"
  if(any(fil)) {
    new$fileName = args$fileName
    was = c(which(fil), was)
  } else {
    was = c(0, was)
  }
  attr(new, "was") <- was 
  attr(new, "from") <- from
  return(new)
}
