.pkgenv <- new.env(emptyenv())

.onLoad <- function(libname, pkgname) {
  data_avl <- requireNamespace("IFCdata", quietly = TRUE)
  .pkgenv[["data_avl"]] <- data_avl
  invisible()
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- c(paste("Package 'IFC' version", packageVersion("IFC")),
           "\nType `citation(\"IFC\")` for citing this R package in publications.")
  if(!.pkgenv[["data_avl"]]) {
    msg <- c(msg,
             "\nTo use examples in this package, you should install 'IFCdata' package.",
             "\nTo install 'IFCdata' package, run `install.packages('IFCdata', 'https://gitdemont.github.io/IFCdata', type = 'source')`.")
  }
  packageStartupMessage(paste(strwrap(msg), collapse = "\n"))      
  invisible()
}