#' @title IFC Population Coercion
#' @description
#' Helper to build a list to allow population export.
#' @param name name of the population.
#' @param type type of population. Either "B", "C", "G" or "T" for Base, Combined, Graphical or Tagged, respectively.\cr
#' If missing, the default, 'type' will be deduced from other parameters.
#' If 'name' is "All" type will be "B". Otherwise, if 'fx' is given type will be "G".
#' Otherwise, "T", if 'definition' is missing but not 'obj' or "C" if 'definition' is not missing.
#' @param base which population is base on. Default is base='All'. Only needed when type = "T".
#' @param color color of the population. See \code{\link{paletteIFC}} for allowed colors. If not provided, will be sampled. 
#' @param lightModeColor lightModeColor of the population. See \code{\link{paletteIFC}} for allowed colors. If not provided, will be sampled.
#' @param style style of the population. Either 20, 4, 3, 1, 5, 0, 2, 18, 15, 17, respectively for:
#' "Simple Dot", "Cross", "Plus", "Empty Circle", "Empty Diamond", "Empty Square",
#' "Empty Triangle", "Solid Diamond", "Solid Square", "Solid Triangle".
#' @param region Only if type='G'. Name of the region defining the population.
#' @param fx Only needed if type='G'. Name of the x-feature defining the population.
#' @param fy Only needed if type='G' and only if region is defined in 2D. Name of the y-feature defining the population. 
#' @param definition Only needed if type='C'. Parameters defining the population.
#' @param obj Only needed if type='T'. Either a:\cr
#' -Logical vector of same length as 'All' population indicating if a cell belongs to the population or not.\cr
#' -Numeric Vector of indices of cells that belongs to the population.
#' @param ... Other arguments to be passed.
#' @return a list containing all population information.
#' @export
buildPopulation <- function(name, type, base="All", color, lightModeColor, style, region, fx, fy, definition, obj, ...) {
  dots = list(...)
  if(missing(name)) stop("'name' can't be missing")
  assert(name, len=1, typ="character")
  if(missing(type)) {
    type = ""
    if(missing(fx)) {
      # if(missing(definition) & !missing(obj)) type="T"
      # if(!missing(definition) & missing(obj)) type="C"
      if(missing(definition)) {
        if(!missing(obj)) type="T"
      } else{
        type="C"
      }
    } else {
      if(missing(obj) & missing(definition)) type="G"
    }
    if(name=="All") type="B"
  }
  assert(type, len=1, alw=c("B","C","G","T"))
  if(name=="All" && type!="B") stop("when 'name' is \"All\", 'type' has to be \"B\"")
  if(type=="B") {
    return(list("name"="All","type"="B","base"="","color"="White","lightModeColor"="Black","style"="Simple Dot"))
  }
  if(missing(color)) {
    if(missing(lightModeColor)) {
      tmp = sample(nrow(paletteIFC("")),1)
    } else {
      assert(lightModeColor, len=1, alw=unlist(paletteIFC("")))
      tmp = which(paletteIFC("")%in%lightModeColor, arr.ind=TRUE)[1]
      if(is.na(tmp)) tmp = sample(nrow(paletteIFC("")),1) 
    }
    color = paletteIFC("")$color[tmp]
    lightModeColor = paletteIFC("")$lightModeColor[tmp]
  } else {
    if(color%in%paletteIFC("")$color_R) color = paletteIFC("")$color[color==paletteIFC("")$color_R][1]
    assert(color, len=1, alw=paletteIFC("palette"))
  }
  if(missing(lightModeColor)) {
    if(missing(color)) {
      tmp = sample(nrow(paletteIFC("")),1)
    } else {
      assert(color, len=1, alw=unlist(paletteIFC("")))
      tmp = which(color==paletteIFC(""), arr.ind=TRUE)[1]
      if(is.na(tmp)) tmp = sample(nrow(paletteIFC("")),1)
    }
    color = paletteIFC("")$color[tmp]
    lightModeColor = paletteIFC("")$lightModeColor[tmp]
  } else {
    if(lightModeColor%in%paletteIFC("")$lightModeColor_R) lightModeColor = paletteIFC("")$lightModeColor[lightModeColor==paletteIFC("")$lightModeColor_R][1]
    assert(lightModeColor, len=1, alw=paletteIFC("palette"))
  }
  tmp_style = c(20, 4, 3, 1, 5, 0, 2, 18, 15, 17)
  names(tmp_style)=c("Simple Dot","Cross","Plus","Empty Circle","Empty Diamond","Empty Square","Empty Triangle","Solid Diamond","Solid Square","Solid Triangle")
  if(missing(style)) {
    style = names(sample(tmp_style, 1))
  } else {
    if(style%in%tmp_style) style = names(tmp_style[which(style==tmp_style)][1])
    if(style%in%names(tmp_style)) {
      style = names(tmp_style[which(style==names(tmp_style))][1])
    } else {
      style = names(sample(tmp_style, 1))
    }
  }
  assert(style, len=1, alw=names(tmp_style))
  if(type=="T") {
    if(missing(obj)) stop("'obj' can't be missing when type='T'")
    if(length(obj)!=0) if(!any(class(obj)%in%c("logical","numeric","integer"))) stop("when type='T', 'obj' has to be a logical or a numeric/integer vector")
    return(list("name"=name,"type"=type,"base"="All","color"=color,"lightModeColor"=lightModeColor,"style"=style,"obj"=obj))
  }
  if(type=="G") {
    if(missing(region)) stop("'region' can't be missing when type='G'")
    assert(region, len=1, typ="character")
    if(missing(fx)) stop("'fx' can't be missing when type='G'")
    assert(fx, len=1, typ="character")
    ret = list("name"=name,"type"=type,"base"=base,"color"=color,"lightModeColor"=lightModeColor,"style"=style,"region"=region,"fx"=fx)
    if(missing(fy)) return(ret)
    if(length(fy)==0) return(ret)
    assert(fy, len=1, typ="character")
    return(c(ret, list("fy"=fy)))
  }
  if(type=="C") {
    if(missing(definition)) stop("'definition' can't be missing when type='C'")
    assert(definition, len=1, typ="character")
    return(list("name"=name,"type"=type, "base"="All", "color"=color,"lightModeColor"=lightModeColor,"style"=style,"definition"=definition))
  }
}
