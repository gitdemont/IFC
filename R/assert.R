#' @title Assert that Certain Conditions are Met
#' @description
#' Ensures that a variable respects several parameters
#' @param x variable to test
#' @param len integer vector of allowed length for x. Default is NULL, for not checking this parameter.
#' @param cla character vector of allowed classes of x. Default is NULL, for not checking this parameter.
#' @param typ character vector of allowed types of x. Default is NULL, for not checking this parameter.
#' @param alw allowed values for x. Default is NULL, for not checking this parameter.
#' @param fun function to execute when mandatory parameters are not met. Default is "stop". Allowed are "stop","warning","message","return".
#' @details /!\ alw parameter when used should be coercible to a logical, integer, numeric, complex or character vector. Otherwise, an error will be thrown.
#' @keywords internal
assert = function(x, len=NULL, cla=NULL, typ=NULL, alw=NULL, fun="stop") {
  foo = cpp_assert(x, len=len, cla=cla, typ=typ, alw=alw, fun=fun)
  if(!any(foo)) return(invisible(NULL))
  # since some conditions are not met everything is recomputed to allow user to check what unmet parameter(s).
  ele = c(len = NULL, cla = NULL, typ = NULL, alw = NULL)
  if(foo[1]) {
    tmp = len%in%length(x)
    rejected = x[!tmp]
    ele["len"] = paste0("len: is of length [",paste0(length(x), collapse=","),"]. Allowed length",ifelse(length(len)>1,"s are"," is"),": ", paste0(len,collapse=","))
  }
  if(foo[2]) {
    tmp = cla%in%class(x)
    ele["cla"] = paste0("cla: ['",paste0(class(x), collapse="','"),"'] not of required class",ifelse(length(cla)>1,"es",""),": ", paste0(paste0("'",cla[!tmp],"'"), collapse=" & "))
  }
  if(foo[3]) {
    tmp = typ%in%typeof(x)
    ele["typ"] = paste0("typ: [",paste0(typeof(x), collapse=","),"] not of required type",ifelse(length(typ)>1,"s",""),": ", paste0(paste0("'",typ[!tmp],"'"), collapse=" & "))
  }
  if(foo[4]) {
    tmp = x%in%alw
    rejected = x[!tmp]
    ele["alw"] = paste0("alw: [",paste0(head(rejected,5), collapse=","),ifelse(length(rejected)>5,", ...",""),"] ",ifelse(length(rejected)>1,"are","is")," not allowed. Allowed values are: ", paste0(alw,collapse=","))
  }
  args = list(paste0(paste0("'", paste0(as.character(substitute(x)),collapse="$"),"':\n"), paste(" -", ele, collapse = "\n")))
  if(fun == "warning") args = c(args, call.=FALSE, immediate.=TRUE)
  do.call(what = fun, args = args)
}
