#' @title Input Identification
#' @description
#' Helper that identifies input arguments thanks to their IFC classes even if they are not or mis named.
#' @param entries arguments from the function \code{\link{whoami}} is called.
#' /!\ \code{\link{whoami}} MUST be called explicitly this way: whoami(entries = as.list(match.call())).
#' @param search a named list of classes to search for entries.
#' @param reinit whether to reinitialize arguments to their default values in called environment. Default is TRUE.
#' @return a list whose members are 'fileName': value of fileName if provided as a named argument 
#' in entries and all classes defined in 'search'
#' @keywords internal
whoami = function(entries = as.list(match.call()),
                  search = list(info = "IFC_info", 
                                 param = "IFC_param",
                                 offsets = "IFC_offset"),
                  reinit = TRUE) {
  eval_from = parent.frame(2)
  entry1 = entries[[1]]
  if(typeof(entry1) == "closure") {
    from = entry1
  } else {
    from = as.character(entry1)
  }
  args = entries[-1]
  L = length(args)
  
  # empty 
  if(L == 0) {
    new = list(fileName == list())
    attr(new, "was") <- 0
    attr(new, "from") <- from
    return(new)
  }
  
  # retrieve arguments values
  val = lapply(1:L, FUN=function(i) {
    switch(typeof(args[[i]]),
           "symbol" = return(get(as.character(args[[i]]), eval_from)),
           "language" = return(eval(args[[i]], eval_from))
    )
    return(args[[i]])
  })
  
  classes = c(list("fileName" = "fileName"), search)
  # identify arguments of classes we look for
  found = sapply(1:L, FUN = function(i) {
    sapply(classes, FUN = function(k) inherits(x = val[[i]], what = k, which = FALSE))
  })
  
  # count how many times a classes is present
  times = apply(found, 1, sum)
  
  # error if at least 2 arguments are with same class
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
    stop(paste0("\nIn ", entry1, ":\n", msg, collapse = "\n"))
  }
  
  # if not error identify who is who and where it was
  iam = lapply(1:length(classes), FUN=function(i) {
    return(which(found[i, ] == 1))
  })
  was = rep(0, times = length(classes))
  was[times == 1] <- unlist(iam)
  new = lapply(iam, FUN = function(i) {
    if(length(i) != 0) return(val[[i]])
  })
  names(new) = names(classes)
  
  # add fileName if it is a named argument
  fil = names(args) %in% "fileName"
  if(any(fil)) {
    new$fileName = args$fileName
    was = c(which(fil), was)
  } else {
    was = c(0, was)
  }
  attr(new, "was") <- was 
  attr(new, "from") <- from
  
  # reinit to arguments of searched classes found to their default value
  if(reinit) {
    call_from = parent.frame(1)
    form = formals(fun = from)
    mism = na.omit(names(args)[was])
    sapply(mism, FUN = function(x) assign(x = x, value = form[[x]], inherits = FALSE, envir = call_from, immediate = TRUE))
  }
  return(new)
}
