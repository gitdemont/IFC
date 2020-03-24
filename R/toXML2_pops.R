#' @title IFC_pops XML Conversion
#' @description 
#' Helper to convert populations (IFC_pops object) to XML nodes.
#' @param pops an IFC_pops object.
#' @param verbose whether to display message about current action. Default is FALSE.
#' @param display_progress whether to display a progress bar. Default is TRUE.
#' @param pb_title character string, giving the window title on the dialog box. Default is "".
#' @return a xml_node.
#' @keywords internal
toXML2_pops = function(pops, verbose = FALSE, display_progress = TRUE, pb_title = "") {
  assert(verbose, alw = c(TRUE, FALSE))
  if(verbose) message("creating pops node")
  assert(pops, cla = "IFC_pops")
  if(length(pops)==0) return(xml_new_node(name = "Pops", text = ""))
  assert(display_progress, alw = c(TRUE, FALSE))
  assert(pb_title, len = 1, typ = "character")
  tmp_style = c(20, 4, 3, 1, 5, 0, 2, 18, 15, 17)
  names(tmp_style)=c("Simple Dot","Cross","Plus","Empty Circle","Empty Diamond","Empty Square","Empty Triangle","Solid Diamond","Solid Square","Solid Triangle")
  L = length(pops)
  if(display_progress) {
    if(.Platform$OS.type == "windows") {
      pb = winProgressBar(title = pb_title, label = "information in %", min = 0, max = 100, initial = 1)
      pb_fun = setWinProgressBar
    } else {
      pb = txtProgressBar(title = pb_title, label = "information in %", min = 0, max = 100, initial = 1, style = 3)
      pb_fun = setTxtProgressBar
    }
    on.exit(close(pb))
    pops_nodes = xml_new_node(name = "Pops", .children = lapply(1:L, FUN=function(i_pop) {
      k=i_pop/L*100
      pb_fun(pb, value = k, label = sprintf("Creating pops nodes %.0f%%", k))
      pop = pops[[i_pop]]
      if(pop$color=="Cyan4") pop$color <- "Teal"
      if(pop$lightModeColor=="Cyan4") pop$lightModeColor <- "Teal"
      if(pop$color=="Green4") pop$color <- "Green"
      if(pop$lightModeColor=="Green4") pop$lightModeColor <- "Green"
      if(pop$color=="chartreuse") pop$color <- "Lime"
      if(pop$lightModeColor=="chartreuse") pop$lightModeColor <- "Lime"
      pop$style <- names(which(pop$style == tmp_style))[1]
      switch(pop$type,
             "B" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style")])
             },
             "C" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "definition")])
             },
             "G" = {
               if(length(pop$fy)==0) {
                 xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "region", "fx")])
               } else { 
                 xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "region", "fx", "fy")])
               }
             },
             "T" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style")],
                            .children = lapply(cpp_num_to_string(which(pop$obj)-1), FUN=function(o) {
                              paste0('ob O="', o, '"')
                            }))
             })
    }))
  } else {
    pops_nodes = xml_new_node(name = "Pops", .children = lapply(1:L, FUN=function(i_pop) {
      pop = pops[[i_pop]]
      if(pop$color=="Cyan4") pop$color <- "Teal"
      if(pop$lightModeColor=="Cyan4") pop$lightModeColor <- "Teal"
      if(pop$color=="Green4") pop$color <- "Green"
      if(pop$lightModeColor=="Green4") pop$lightModeColor <- "Green"
      if(pop$color=="chartreuse") pop$color <- "Lime"
      if(pop$lightModeColor=="chartreuse") pop$lightModeColor <- "Lime"
      pop$style <- names(which(pop$style == tmp_style))[1]
      switch(pop$type,
             "B" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style")])
             },
             "C" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "definition")])
             },
             "G" = {
               if(length(pop$fy)==0) {
                 xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "region", "fx")])
               } else {
                 xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style", "region", "fx", "fy")])
               }
             },
             "T" = {
               xml_new_node(name = "Pop", attrs = pop[c("name", "type", "base", "color", "lightModeColor", "style")],
                            .children = lapply(cpp_num_to_string(which(pop$obj)-1), FUN=function(o) {
                              paste0('ob O="', o, '"')
                            }))
             })
    }))
  }
  return(pops_nodes)
}
