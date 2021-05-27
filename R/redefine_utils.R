################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2021 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2021                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Jean-Pierre Marolleau, Loïc Garçon,       #
#                     CHU Amiens                                               #
#                                                                              #
# DISCLAIMER:                                                                  #
# -You are using this package on your own risk!                                #
# -We do not guarantee privacy nor confidentiality.                            #
# -This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or #
# contributors be liable for any direct, indirect, incidental, special,        #
# exemplary, or consequential damages (including, but not limited to,          #
# procurement of substitute goods or services; loss of use, data, or profits;  #
# or business interruption) however caused and on any theory of liability,     #
# whether in contract, strict liability, or tort (including negligence or      #
# otherwise) arising in any way out of the use of this software, even if       #
# advised of the possibility of such damage.                                   #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with IFC. If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################


################################################################################
# This file contains several functions under development which are intended    #
# to rename and propagate renaming into IFC_data objects                       #
################################################################################

#' @title IFC Mask Coercion
#' @description
#' Helper to build a mask.
#' @param masks an `IFC_masks` object or a data.frame containing masks definition and name. Default is missing.
#' @param images a data.frame containing images definition. Default is missing.
#' @param definition whether to compute masks definition or masks names. Default is TRUE, to compute definition.
#' @param to_match_mask a vector of masks names to use for matching 'masks' names. Default is NULL
#' @param to_replace_mask a vector of masks names to use for replacing 'masks' names. Default is NULL
#' @param to_match_image a vector of images names to use for matching 'images' names. Default is NULL
#' @param to_replace_image a vector of images names to use for replacing 'images' names. Default is NULL
#' @param ... Other arguments to be passed.
#' @details function that can either change definition or name
#' it will be used in a loop to incorporate new definition.
#' causing name to be changed. allowing new redefinition of mask 
#' dependent on other mask to be changed, and so on
#' when a mask or an image name is not found because not yet defined 
#' an error is raised and catched
#' TODO maybe this error catching causes some overhead to be improved
#' @return a vector of masks definition or names depending on 'definition' parameter.
#' @keywords internal
buildMask = function(masks,
                     images, 
                     definition = TRUE, 
                     to_match_mask = NULL,
                     to_replace_mask = NULL,
                     to_match_image = NULL,
                     to_replace_image = NULL, ...) {
  mask_avl = list("Threshold" = list(def = list("mask", "image", structure("integer", range = c(0,100))), 
                                     ord = c(1,2,3)), # mask, image, integer [0-100]
                  "Dilate" = list(def = list("mask", structure("integer", range = c(0,19))), 
                                  ord = c(1,2)), # mask, integer [0-19]
                  "Erode" = list(def = list("mask", structure("integer", range = c(0,19))), 
                                 ord = c(1,2)), # mask, integer [0-19]
                  "Morphology" = list(def = list("mask", "image"), 
                                      ord = c(1,2)),  # mask, image
                  "Spot" = list(def = list("mask", "image",
                                           structure("integer", map_val = c(0,1), map_txt = c("Bright", "Dark")),
                                           structure("numeric", range = c(0,100), round = 2),
                                           structure("integer", range = c(0,30)),
                                           structure("integer", range = c(0,30))), 
                                ord = c(1,2,3,4,5,6)), # mask, image, integer* [0-1], numeric [0.00-100.00], integer [0-30], integer[0-30] #* 0=Bright, 1=Dark
                  "AdaptiveErode"= list(def = list("mask", "image", structure("integer", range = c(0,100))), 
                                        ord = c(1,2,3)), #mask, image, integer [0-100]
                  "Component" =  list(def = list("mask", "image", "image",
                                                 structure("integer", 
                                                           map_val = c(0,6,61,64,63,
                                                                       24,25,34,26,
                                                                       20,17,15,16,18,53,
                                                                       29),
                                                           map_txt = c("Area", "Aspect Ratio", "Circularity", "Thickness Max", "Thickness Min",
                                                                       "Bright Detail Intensity R3", "Bright Detail Intensity R7", "Contrast", "Gradient RMS", 
                                                                       "Intensity", "Max Pixel", "Mean Pixel", "Median Pixel", "Min Pixel", "Saturation Count",
                                                                       "Bright Detail Similarity R3")),
                                                 structure("integer", map_val = c(0,1), map_txt = c("Descending", "Ascending")),
                                                 structure("integer", range = c(0,100)),
                                                 structure("integer", range = c(4095))), 
                                      ord = c(6,4,1,2,3,5)),  #mask, image, image, integer* [], integer** [0-1], integer***[1-100], integer[4095]
                  #*  0=Area
                  #*  6=Aspect Ratio
                  #* 24=Bright Detail Intensity R3
                  #* 25=Bright Detail Intensity R7
                  #* 29=Bright Detail Similarity R3
                  #* 61=Circularity
                  #* 34=Contrast
                  #* 26=Gradient RMS
                  #* 20=Intensity
                  #* 17=Max Pixel
                  #* 15=Mean pixel
                  #* 16=Median Pixel
                  #* 18=Min Pixel
                  #* 53=Saturation Count
                  #* 64=Thickness Max
                  #* 63=Thickness Min
                  #** 0=Descending, 1=Ascending
                  #*** rank
                  "Fill" = list(def = "mask", 
                                ord = 1),# mask
                  "Inspire" = list(def = list("mask","image"), 
                                   ord = c(1,2)),# mask
                  "Intensity" = list(def = list("mask", "image", structure("integer", range = c(0,4095)), structure("integer", range = c(0,4095))), 
                                     ord = c(1,2,3,4)), #mask, image, integer [0-4095], integer [0-4095]
                  "Interface" = list(def = list("mask", "mask", structure("integer", range = c(0,19))), 
                                     ord = c(1,2,3)), #mask, mask, integer [0-19]
                  "LevelSet" = list(def = list("mask", "image", 
                                               structure("integer", map_val = c(1,2,3,4), map_txt = c("Combined","Dim","Middle","Bright")),
                                               structure("numeric", range = c(0,10), round = 0)), 
                                    ord = c(1,2,3,4)), #mask, image, integer* [1-4], numeric [0-10] #* 1=Combined, 2=Dim, 3=Middle, 4=Bright
                  "Object"= list(def = list("mask","image",
                                            structure("integer", map_val = c(1), map_txt = c("Tight"))), 
                                 ord = c(1,2,3)),# mask #mask, image, integer* [1] #*1=Tight
                  "Peak"= list(def = list("mask","image",
                                          structure("integer", map_val = c(0,1), map_txt = c("Bright","Dark")),
                                          structure("numeric", range = c(0,100), round = 2)), 
                               ord = c(1,2,3,4)), #mask, image, integer* [0-1], numeric [0.00-100.00] #* 0=Bright, 1=Dark
                  "Range"= list(def = list("mask",
                                           structure("integer", range = c(0,5000)),
                                           structure("integer", range = c(0,5000)),
                                           structure("numeric", range = c(0,1), round = 1),
                                           structure("numeric", range = c(0,1), round = 1)), 
                                ord = c(1,2,3,4,5)), #mask, integer [0-500], numeric [0.0-1.0]
                  "Skeleton"= list(def = list("mask","image",
                                              structure("integer", map_val = c(0,1), map_txt = c("Thin","Tight"))), 
                                   ord = c(1,2,3)), #mask, image, integer* [0-1] #* 0=Thin, 1=Thick
                  "System"= list(def = list("mask","image",
                                            structure("integer", range = c(0,100))), 
                                 ord = c(1,2,3)), #mask, image, integer [0-100]
                  "Valley"= list(def = list("mask","image",
                                            structure("integer", range = c(0,10))), 
                                 ord = c(1,2,3)), #mask, image, integer [1-10]
                  "Watershed"= list(def = list("mask","image",
                                               structure("integer", map_val = c(0,1), map_txt = c("none","IW")),
                                               structure("integer", map_val = c(0,1), map_txt = c("none","Inv")),
                                               structure("integer", range = c(0,10)),
                                               structure("integer", range = c(0,20))), 
                                    ord = c(1,2,3,4,5,6)) #mask, image, integer* [0-1], integer** [0-1], integer*** [0-10], integer**** [1-20] 
                  #* 0=none, 1=IW
                  #** 0=none, 1= Inv, is for Valleys, 1 is for Peaks
                  #*** smoothing
                  #**** thickness
  )
  M = masks
  I = images
  
  sapply(1:nrow(M), FUN = function(i_row) {
    comp = strsplit(M$def[i_row], split="|", fixed = TRUE)[[1]]
    first = comp[1]
    tmp = which(first == names(mask_avl))
    if(length(tmp) == 0) {
      bool = try(validate_bool(x = comp, all_names = c(sprintf("M%02i", I$physicalChannel), M$name)), silent = TRUE)
      if(definition) {
        comp = sapply(comp, FUN = function(x) {
          foo = x == to_match_mask
          if(any(foo)) return(to_replace_mask[foo])
          return(x)
        })
        return(paste0(comp, collapse = "|"))
      } else{
        if(inherits(bool, what = "try-error")) return("")
        return(paste0(bool, collapse = " ")) 
      }
    } else {
      comp = comp[-1]
      def = mask_avl[[tmp]]
      check_length = try(assert(comp, len = length(def$def)), silent = TRUE)
      if(inherits(x = check_length, what = "try-error")) {
        if(first == "Spot") {
          tmp = rep(TRUE, length(def$def))
          tmp[4] <- FALSE
          def$def = def$def[tmp]
          def$ord = def$ord[-length(def$ord)]
          assert(comp, len = length(def$def))
        } else {
          stop("bad mask definition length: ",M$def[i_row])
        }
      }
      if(definition) {
        val = sapply(1:length(def$def), FUN =function(i) {
          typ = def$def[[i]]
          switch(typ, 
                 "image" = {
                   tmp = comp[i] == to_match_image
                   if(any(tmp)) return(to_replace_image[tmp]) 
                 },
                 "mask" = {
                   tmp = comp[i] == to_match_mask
                   if(any(tmp)) return(to_replace_mask[tmp]) 
                 })
          return(comp[i])
        })
        val = paste0(c(first,val), collapse = "|")
      } else {
        val = try(sapply(1:length(def$def), FUN =function(i) {
          typ = def$def[[i]]
          switch(typ,
                 "mask"= {
                   if(!(comp[i] %in% c(sprintf("M%02i", I$physicalChannel), M$name))) stop("can't find mask: ", comp[i])
                   return(comp[i])
                 },
                 "image" = {
                   if(!(comp[i] %in% c(sprintf("M%02i", I$physicalChannel), I$name))) stop("can't find image: ", comp[i])
                   return(comp[i])
                 },
                 "integer" = {
                   foo = as.integer(comp[i])
                   tmp = attr(typ, "range")
                   if(length(tmp) != 0) {
                     foo = foo[(foo >= tmp[1]) & (foo <= tmp[2])]
                     assert(foo, len = 1)
                   }
                   tmp = attr(typ, "map_val")
                   if(length(tmp) != 0) {
                     assert(foo, alw = tmp)
                     foo = attr(typ, "map_txt")[tmp == foo]
                   }
                   return(foo)
                 },
                 "numeric" = {
                   foo = as.numeric(comp[i])
                   tmp = attr(typ, "range")
                   if(length(tmp) != 0) {
                     foo = foo[(foo >= tmp[1]) & (foo <= tmp[2])]
                     assert(foo, len = 1)
                   }
                   tmp = attr(typ, "round")
                   if(length(tmp) != 0) {
                     foo = round(foo, tmp)
                   }
                   return(foo)
                 },
                 stop("not compatible type")
          )
        }), silent = TRUE)
        if(inherits(val, what = "try-error")) {
          return("")
        } else {
          val = val[def$ord]
          switch(first,
                 "Intensity" = {
                   val = c(val[1:2], paste0(val[3], "-", val[4]))
                 },
                 "Range" = {
                   val = c(val[1], paste0(val[2], "-", val[3]), paste0(val[4], "-", val[5]))
                 },
                 "Watershed" = {
                   if(val[4] == "none") val = val[-4]
                   # TODO strange behavior of Watershed which does not report line width in name when
                   # mask is not intensity weighted despite being well kept in description
                   if(val[3] == "none") val = val[1]
                 },
                 "Component" = {
                   if(comp[4] != 29) val = val[-5]
                   if(comp[4] %in% c(0,6,61,64,63)) val = val[-4]
                 }
          )
          return(paste0(first, "(", paste0(val, collapse = ", "), ")"))
        }
      }
    }
  })
}

#' @title IFC_masks Mask Redefinition
#' @description
#' Helper to rename a mask within IFC_masks.
#' @param masks an `IFC_masks` object or a data.frame containing masks definition and name. Default is missing.
#' @param images a data.frame containing images definition. Default is missing.
#' @param to_match_mask a string with a mask name to use for matching 'masks' names. Default is NULL
#' @param to_replace_mask a string of mask name to use for replacing 'masks' names. Default is NULL
#' @param ... Other arguments to be passed.
#' @return a vector of masks definition or names depending on definition.
#' @keywords internal
redefine_masks_mask <- function(masks, images, to_match_mask = NULL, to_replace_mask = NULL, ...) {
  if(length(to_match_mask) != 1 || length(to_replace_mask) != 1)  stop("'to_match_mask' should be of length 1")
  M = masks
  M$name[M$name == to_match_mask] <- to_replace_mask
  old_names = M$name

  masks_names = buildMask(masks, images, definition = FALSE)
  # at this step if any computed masks_names if "" it means that M$def is not correct
  tmp = masks_names == ""
  if(any(tmp)) stop("mask",ifelse(sum(tmp) > 1,"s are"," is")," not well defined:\n\t -",paste0(M$def[tmp],collapse="\n\t -"))
  
  # some masks name may have been customized, so we don't want
  # to change them whatever happens, so we store which we can modify
  to_modify = M$name == masks_names
  
  M$name <- ""
  M$name[!to_modify]  <- old_names[!to_modify]
  new_names <- NULL
  
  # loop until all names are defined
  while(any(M$name == "")) {
    M$def = buildMask(M, images, TRUE, to_match_mask = to_match_mask, to_replace_mask = to_replace_mask)
    new_names = buildMask(M, images, FALSE, to_match_mask = to_match_mask, to_replace_mask = to_replace_mask)
    # we store mapping beetwen initial masks names and new masks names that have been replaced
    # thanks to new masks definition taking into account images names replacement
    to_match_mask = old_names[to_modify & (new_names != "")]
    to_replace_mask = new_names[to_modify & (new_names != "")]
    # just to be sure that while loop can be escaped
    # in case M$name can't be correctly computed more than one time
    if(identical(M$name[to_modify], new_names[to_modify])) break 
    M$name[to_modify] = new_names[to_modify]
  }
  M$name[to_modify] = new_names[to_modify]
  class(M) <- class(masks)
  return(list(masks = M, images = images))
}

#' @title IFC_masks Image Redefinition
#' @description
#' Helper to rename images within masks definition.
#' @param masks an `IFC_masks` object or a data.frame containing masks definition and name. Default is missing.
#' @param images a data.frame containing images definition. Default is missing.
#' @param new_images_names a vector of image name to use for replacing 'images' names. Default is images$name
#' @param ... Other arguments to be passed.
#' @return a list whose members are:
#' -masks, an `IFC_masks` object or data.frame containing masks definition and name.
#' -images, a data.frame containing images definition. 
#' @keywords internal
redefine_masks_image <- function(masks, images, new_images_names = images$name, ...) {
  if(length(new_images_names) != length(images$name))  stop("'new_images_names' should be of same length of images$name")
  
  dots = list(...)
  
  # all allowed mask functions, with indication on how they are defined + 
  # how they are ordered to compose regular name
  
  # masks ordering
  # mask_comp = lapply(M$def, FUN = function(m) {
  #   all_names = c(sprintf("M%02i", I$physicalChannel), M$name)
  #   operators = c("And", "Or", "Not", "(", ")")
  #   split = splitn(definition = m, all_names = all_names , operators = operators)
  #   split[split %in% all_names]
  # })
  # names(mask_comp) = M$name
  # 
  # i=1; l=length(mask_comp)
  # while(i<l) {
  #   m = mask_comp[[i]]
  #   index = setdiff(c(sprintf("M%02i", I$physicalChannel), m),"")
  #   index = unlist(lapply(index, function(x) which(x==names(mask_comp))))
  #   index = index[index>i]
  #   if(length(index)!=0) {
  #     mask_comp = c(mask_comp[index],mask_comp[setdiff(1:l,index)])
  #     i=1
  #   } else {
  #     i=i+1
  #   }
  # }
  # 
  # M = M[sapply(1:nrow(M), FUN = function(i_row) which(M$name == names(mask_comp)[i_row])), ]
  
  M = masks
  I = images
  
  # we store current M$name
  old_names = M$name
  # we compute expected M$name based on M$def
  masks_names = buildMask(M, I, FALSE)
  # at this step if any computed masks_names if "" it means that M$def is not correct
  tmp = masks_names == ""
  if(any(tmp)) stop("mask",ifelse(sum(tmp) > 1,"s are"," is")," not well defined:\n\t -",paste0(M$def[tmp],collapse="\n\t -"))
  
  # some masks name may have been customized, so we don't want
  # to change them whatever happens, so we store which we can modify
  to_modify = masks$name == masks_names
  
  to_match = images$name
  to_replace = new_images_names
  definition = !identical(to_match, to_replace)
  if(definition) {
    to_keep = to_match != to_replace
    to_match_image = to_match[to_keep]
    to_replace_image = to_replace[to_keep]
  }
  
  # -when change_definition is FALSE, the above function f() returns computed mask name
  # -when change_definition is TRUE, the function returns mask definition
  ## where images in definition are replaced with to_replace input from the user
  ## where masks in definition are replaced with to_replace_masks which are changed by looping f(FALSE) f(TRUE)
  if(definition) {
    M$name <- ""
    M$name[!to_modify]  <- old_names[!to_modify]
    I$name <- to_replace
    new_names <- NULL
    to_match_mask <- NULL
    to_replace_mask <- NULL
    M$def = buildMask(M, I, TRUE,
                      to_match_mask = to_match_mask, to_replace_mask = to_replace_mask,
                      to_match_image = to_match_image, to_replace_image = to_replace_image)
    # loop until all names are defined
    while(any(M$name == "")) {
      M$def = buildMask(M, I, TRUE, 
                        to_match_mask = to_match_mask, to_replace_mask = to_replace_mask,
                        to_match_image = to_match_image, to_replace_image = to_replace_image)
      new_names = buildMask(M, I, FALSE, to_match_mask = to_match_mask, to_replace_mask = to_replace_mask,
                            to_match_image = to_match_image, to_replace_image = to_replace_image)
      # we store mapping beetwen initial masks names and new masks names that have been replaced
      # thanks to new masks definition taking into account images names replacement
      to_match_mask = old_names[to_modify & (new_names != "")]
      to_replace_mask = new_names[to_modify & (new_names != "")]
      # just to be sure that while loop can be escaped
      # in case M$name can't be correctly computed more than one time
      if(identical(M$name[to_modify], new_names[to_modify])) break 
      M$name[to_modify] = new_names[to_modify]
    }
    M$name[to_modify] = new_names[to_modify]
    # if any M$name is "" it means that the loop was escaped because new definition can't be applied
    if(any( M$name == "")) stop("can't modify masks with new names")
    class(M) <- class(masks)
    class(I) <- class(images)
    return(list(masks = M, images = I))
  } else {
    M$name[to_modify] = masks_names[to_modify]
    class(M) <- class(masks)
    class(I) <- class(images)
    return(list(masks = M, images = I))
  }
}

#' @title IFC_masks Redefinition
#' @description
#' Helper to modify features_def according to masks or images redefinition
#' @param masks an `IFC_masks` object or a data.frame containing masks definition and name. Default is missing.
#' @param images a data.frame containing images definition. Default is missing.
#' @param new_images_names a vector of image name to use for replacing 'images' names. Default is images$name
#' @param to_match_mask a string with a mask name to use for matching 'masks' names. Default is NULL
#' @param to_replace_mask a string of mask name to use for replacing 'masks' names. Default is NULL
#' @param ... Other arguments to be passed.
#' @return a list whose members are:
#' -masks, an `IFC_masks` object or data.frame containing masks definition and name.
#' -images, a data.frame containing images definition. 
#' @keywords internal
redefine_masks <- function(masks, images, new_images_names = images$name, to_match_mask = NULL, to_replace_mask = NULL, ...) {
  ans = redefine_masks_image(masks, images, new_images_names = new_images_names)
  if(length(to_match_mask) != 0) {
    for(i in 1:length(to_match_mask)) {
      ans = redefine_masks_mask(ans$masks, ans$images, to_match_mask = to_match_mask[i], to_replace_mask = to_replace_mask[i])
    }
  }
  return(ans)
}

#' @title Feature Default Name Computation
#' @description
#' Helper to compute default name of a feature.
#' @param feat_def list containing feature definition
#' @return a string with default name
#' @keywords internal
feature_namer <- function(feat_def) {
  def_def <- feat_def
  def_def$def <- gsub("\\|Diameter:\\|1\\|20$", "", def_def$def) # for Ensquared
  def_def$def <- gsub("\\|Granularity:\\|1\\|20$", "", def_def$def) # for Haralick features
  def_def$def <- gsub("\\|\\|True\\|.*$", "", def_def$def) # for Spot count
  def_def$def <- gsub(paste0("(\\|?)(true)(\\|?)"), paste0("\\1","IntensityWeighted","\\3"), def_def$def) # for Delta centroid features
  def_def$def <- gsub(paste0("(\\|?)(false)(\\|?)"), paste0("\\1","\\3"), def_def$def) # for Delta centroid features
  def_def$def <- gsub(paste0("\\|\\|"), "|", def_def$def) # for Delta centroid features
  if(feat_def$type == "combined") {
    gsub(" \\)", ")", gsub("\\( ", "(", gsub("|", " ", def_def$def, fixed = TRUE)))
  } else {
    gsub("|", "_", def_def$def, fixed = TRUE)
  }
}

#' @title IFC_features_def Feature Redefinition
#' @description
#' Helper to rename a feature within IFC_features_def
#' @param masks an `IFC_features_def` objector a list containing features definition. Default is missing.
#' @param to_match_feat a string with a features_def name to use for matching 'features_def' names. Default is NULL
#' @param to_replace_feat a string of features_def name to use for replacing 'features_def' names. Default is NULL
#' @param ... Other arguments to be passed.
#' @return an `IFC_features_def` object, or a list containing features definition
#' @keywords internal
redefine_features_def_feat <- function(features_def, to_match_feat = NULL, to_replace_feat = NULL, ...) {
  if(length(to_match_feat) != 1 || length(to_replace_feat) != 1)  stop("'to_match_feat' should be of length 1")
  comb_operators = comb_operators = c("+", "-", "*", "/", "ABS", "COS", "SIN", "SQR", "SQRT")
  def = features_def
  old_names = sapply(features_def, FUN = function(x) x$name)
  def[[which(old_names == to_match_feat)]]$name <- to_replace_feat
  
  exp_names = sapply(features_def, feature_namer)
  to_modify = old_names == exp_names
  
  # we extract names
  LL = length(def)
  new_names = sapply(def, FUN = function(x) x$name)
  cur_names = rep("", LL)
  while(any("" == cur_names)) {
    def = lapply(1:LL, FUN = function(i_feat) {
      def_def = def[[i_feat]]
      if((def_def$type == "combined")) {
        splitted = splitn(definition = def[[i_feat]]$def, all_names = old_names, operators = comb_operators)
        foo = sapply(splitted, USE.NAMES = FALSE, FUN = function(x) {
          if(any(x == comb_operators)) return(x)
          tmp = x == old_names
          if(any(tmp)) return(new_names[tmp])
          return("")
        })
        if(any(foo == "")) {
          def_def$name = ""
        } else {
          def_def$def = paste0(foo, collapse = "|")
          if(!identical(splitted, foo) && to_modify[i_feat]) {
            new_names[i_feat] <<- c(paste0(foo, collapse = " "))
            def_def$name = new_names[i_feat]
          }
        }
      }
      return(def_def)
    })
    def
    new_names = sapply(def, FUN = function(def_def) def_def$name)
    if(identical(new_names, cur_names)) break
    cur_names = new_names
  }
  names(def) = sapply(def, feature_namer)
  class(def) = class(features_def)
  return(def)
}

#' @title IFC_features_def Mask or Image Redefinition
#' @description
#' Helper modify features_def according to masks or images redefinition
#' @param features_def an `IFC_features_def` object, or a list containing features definition. Default is missing.
#' @param masks an `IFC_masks` object or a data.frame containing masks definition and name. Default is missing.
#' @param images a data.frame containing images definition. Default is missing.
#' @param ... Other arguments to be passed.
#' @return a list whose members are:\cr
#' -features_def, an `IFC_features_def` object, or a list containing features definition\cr
#' -masks, an `IFC_masks` object or a data.frame containing masks definition and name.\cr
#' -images, a data.frame containing images definition.
#' @keywords internal
redefine_features_def_msk_img <- function(features_def, masks, images, ...) {
  dots = list(...)
  # first we check that masks names are valid with input parameters
  # and recompute new masks names applying these parameters
  foo = do.call(what = redefine_masks, args = c(list(masks = masks, images = images), dots))
  M = foo$masks
  I = foo$images
  
  # we define names to change in features definition
  # names to changed can be new masks names but also new images names
  has_changed = M$name != masks$name 
  to_match = images$name
  to_replace = foo$images$name
  to_keep = to_match != to_replace
  to_match_image = c(to_match[to_keep], masks$name[has_changed])
  to_replace_image = c(to_replace[to_keep], M$name[has_changed])
  
  # matches are encapsulated within | |
  to_find = paste0("|",to_match_image,"|")
  to_replace_image = paste0("|",to_replace_image,"|")
  
  # ordering should ensure that names consisted of repeated pattern should be treated 1st
  order_ = order(nchar(to_find), decreasing = TRUE)
  to_find = to_find[order_]
  to_replace_image = to_replace_image[order_]
  
  # we store current features names
  old_names = names(features_def)
  
  # some masks name may have been customized, so we don't want
  # to change them whatever happens, so we compute expected name and compare with actual names
  # if same, it means that they were not customized and can they be modified
  exp_names = sapply(features_def, feature_namer)
  to_modify = old_names == exp_names
  
  # we modify features definition with new masks and new images
  # in addition we modify $name when original name was not customized
  def = features_def
  L = length(to_find)
  LL = length(features_def)
  def = lapply(1:LL, FUN = function(i_feat) {
    def_def <- features_def[[i_feat]]
    sapply(1:L, FUN = function(i_pat) {
      # features definition are encapsulated within | | 
      # to ensure that we can capture matches at the start or end of def_def$def
      # and then 1st and last | are removed
      tmp = gsub(to_find[i_pat], to_replace_image[i_pat], x = paste0("|",def_def$def,"|"), fixed = TRUE) 
      def_def$def <<- substr(tmp, 2, nchar(tmp)-1)
    })
    if(to_modify[i_feat]) def_def$name <- feature_namer(def_def)
    def_def
  })
  
  # now we need to modify features which depend on other features if any
  # we check which features are combined
  comb_operators = c("+", "-", "*", "/", "ABS", "COS", "SIN", "SQR", "SQRT")
  
  # we extract names
  new_names = sapply(def, FUN = function(def_def) def_def$name)
  cur_names = rep("", LL)
  while(any("" == cur_names)) {
    def = lapply(1:LL, FUN = function(i_feat) {
      def_def = def[[i_feat]]
      if((def_def$type == "Combined")) {
        splitted = splitn(definition = def[[i_feat]]$def, all_names = old_names, operators = comb_operators)
        foo = sapply(splitted, USE.NAMES = FALSE, FUN = function(x) {
          if(any(x == comb_operators)) return(x)
          tmp = x == old_names
          if(any(tmp)) return(new_names[tmp])
          return("")
        })
        if(any(foo == "")) {
          def_def$name = ""
        } else {
          def_def$def = paste0(foo, collapse = "|")
          if(!identical(splitted, foo) && to_modify[i_feat]) {
            new_names[i_feat] <<- c(paste0(foo, collapse = " "))
            def_def$name = new_names[i_feat]
          }
        }
      }
      return(def_def)
    })
    new_names = sapply(def, FUN = function(def_def) def_def$name)
    if(identical(new_names, cur_names)) break
    cur_names = new_names
  }
  if(any(new_names == "")) stop("can't rename features")
  names(def) = new_names
  class(def) <- class(features_def)
  return(list(features_def = def,
              masks = M,
              images = I))
}

#' @title IFC_features_def Redefinition
#' @description
#' Helper modify features_def according to masks or images redefinition
#' @param features_def an `IFC_features_def` object, or a list containing features definition. Default is missing.
#' @param masks an `IFC_masks` object or a data.frame containing masks definition and name. Default is missing.
#' @param images a data.frame containing images definition. Default is missing.
#' @param to_match_feat a string with a features_def name to use for matching 'features_def' names. Default is NULL
#' @param to_replace_feat a string of features_def name to use for replacing 'features_def' names. Default is NULL
#' @param ... Other arguments to be passed.
#' @return a list whose members are:\cr
#' -features_def, an `IFC_features_def` object, or a list containing features definition\cr
#' -masks, an `IFC_masks` object or a data.frame containing masks definition and name.\cr
#' -images, a data.frame containing images definition.
#' @keywords internal
redefine_features_def <- function(features_def, masks, images, to_match_feat = NULL, to_replace_feat = NULL, ...) {
  ans = redefine_features_def_msk_img(features_def, masks, images, ...)
  if(length(to_match_feat) != 0) {
    for(i in 1:length(to_match_feat)) {
      ans$features_def = redefine_features_def_feat(features_def = ans$features_def, to_match_feat = to_match_feat[i], to_replace_feat = to_replace_feat[i])
    }
  }
  return(ans)
}

#' @title IFC_data Redefinition
#' @description
#' Helper to modify images, masks, features in a IFC data object when features_def is changed
#' @param obj an `IFC_data` object. Default is missing.
#' @param new_feat_def a list with new definitions from modify_features().
#' @param ... Other arguments to be passed.
#' @return an `IFC_data` object.
#' @keywords internal
redefine_obj <- function(obj, new_feat_def) {
  assert(obj, cla = "IFC_data")

  # we check which features names have been modified
  old_names = names(obj$features_def)
  new_names = names(new_feat_def$features_def)
  to_modify = old_names != new_names
  
  # modify features names
  names(obj$features) = sapply(names(obj$features), FUN = function(x) {
    tmp = (x == old_names) & to_modify
    if(any(tmp)) return(new_names[tmp][1])
    return(x)
  })
  
  # modify graphs
  K = class(obj$graphs)
  obj$graphs = lapply(obj$graphs, FUN = function(g) {
    tmp = (g$f1 == old_names) & to_modify
    if(any(tmp)) {
      if(g$f1 == g$xlabel) g$xlabel <- new_names[tmp][1]
      g$f1 <- new_names[tmp][1]
    }
    if(g$type != "histogram") {
      tmp = (g$f2 == old_names) & to_modify
      if(any(tmp)) {
        if(g$f2 == g$ylabel) g$ylabel <- new_names[tmp][1]
        g$f2 <- new_names[tmp][1]
      }
    }
    return(g)
  })
  class(obj$graphs) <- K
  
  # modify pops
  K = class(obj$pops)
  obj$pops = lapply(obj$pops, FUN = function(p) {
    if("fx" %in% names(p)) {
      tmp = (p$fx == old_names) & to_modify
      if(any(tmp)) p$fx <- new_names[tmp][1]
    }
    if("fy" %in% names(p)) {
      tmp = (p$fy == old_names) & to_modify
      if(any(tmp)) p$fy <- new_names[tmp][1]
    }
    return(p)
  })
  class(obj$pops) <- K
  
  obj$features_def <- new_feat_def$features_def
  obj$description$masks <- new_feat_def$masks
  obj$description$Images <- new_feat_def$images
  return(obj)
}