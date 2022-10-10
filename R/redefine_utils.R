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
  # at this step if any computed masks_names is"" it means that M$def is not correct
  tmp = masks_names == ""
  if(any(tmp)) stop("mask",ifelse(sum(tmp) > 1,"s are"," is")," not well defined:\n\t -",paste0(M$def[tmp],collapse="\n\t -"))
  
  # some masks name may have been customized, so we don't want
  # to change them whatever happens, so we store which we can modify
  to_modify = M$name == masks_names
  # to_modify = rep(T, nrow(M))
  
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
  # to_modify = rep(T, nrow(M))
  
  to_match = images$name
  to_replace = new_images_names
  definition = !identical(to_match, to_replace)
  if(definition) {
    to_keep = to_match != to_replace
    to_match_image = to_match[to_keep]
    to_replace_image = to_replace[to_keep]
  }
  
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
#' Helper to modify features_def according to masks or images redefinition.
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
  if(length(feat_def$split) != 0) {
    if(feat_def$type == "combined") {
      return(paste0(feat_def$split, collapse = " "))
    } else {
      ans = feat_def$split
      ans[ans == "true"] <- "IntensityWeighted"
      ans = ans[!(ans %in% c("false"))]
      ans = gseq(ans, c("Granularity:","1","20")) # TODO maybe modify things here
      ans = gseq(ans, c("Diameter:","1","20"))    #
      ans = gseq(ans, c("False","1","20"))        #
      ans = gseq(ans, c("True","8","4"))          #
      ans = gseq(ans, c("","1","20"))             #
      return(paste0(ans[ans != ""], collapse = "_"))
    }
  }
  ans = feat_def$def
  ans <- gsub("\\|((Diameter|Granularity):)?\\|1\\|20$", "", ans) # for Haralick features and Ensquared
  ans <- gsub("\\|\\|(True|False)?\\|.*$", "", ans) # for Spot count and Spot range
  ans <- gsub("(\\|?)(true)(\\|?)", paste0("\\1","IntensityWeighted","\\3"), ans) # for Delta centroid features
  ans <- gsub("(\\|?)(false)(\\|?)", paste0("\\1","\\3"), ans) # for Delta centroid features
  ans <- gsub("\\|\\|", "|", ans) # for Delta centroid features
  ans <- gsub("\\|+$", "", ans) # remove trailing |
  if(feat_def$type == "combined") {
    gsub(" \\)", ")", gsub("\\( ", "(", gsub("|", " ", ans, fixed = TRUE)))
  } else {
    gsub("|", "_", ans, fixed = TRUE)
  }
}

#' @title IFC_features_def Definition Splitting
#' @description
#' Helper to split a features definitions within IFC_features_def.
#' @param features_def an `IFC_features_def` object or a list containing features definition. Default is missing.
#' @param all_names the names of all allowed names.
#' @param alt_names vector of same length as 'all_names' to use for substitution. It can be used to speed up the process.
#' @param operators operators used. Default is c("And", "Or", "Not", "(", ")").
#' @param m_names names of masks.
#' @param i_names names of images.
#' @param comb_operators operators used for combined features. Default is c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT").
#' @param extr_operators operators used for non combined features. Default is c("true", "false", "True", "False").
#' Those will be added to 'm_names', 'i_names' and all possible features names (as given by the getFromNamespace("featureIFC", "IFC")())
#' @param split string used for splitting. Default is "|".
#' @param force whether to force splitting even if split is detected
#' @return an `IFC_features_def` object, or a list containing features definition
#' @keywords internal
split_feat <- function(features_def,
                       all_names = names(features_def),
                       alt_names = NULL,
                       m_names = sprintf("M%02i",1:12),
                       i_names = sprintf("Ch%02i",1:12),
                       comb_operators = c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT"),
                       extr_operators = c("true", "false", "True", "False"),
                       split = "|",
                       force = FALSE) {
  force = as.logical(force); assert(force, len=1, alw=c(TRUE,FALSE))
  def = features_def
  force = force || any(sapply(def, FUN  = function(x) length(x$split)) == 0)
  all_names_avl = c(unlist(featureIFC(), recursive = FALSE, use.names = FALSE), extr_operators, m_names, i_names)
  alt_names_avl = gen_altnames(all_names_avl)
  if(force) {
    if(missing(alt_names) || (length(alt_names) != length(all_names))) alt_names = gen_altnames(all_names)
    for(i_def in seq_along(def)) {
      if(def[[i_def]]$type == "combined") {
        def[[i_def]]$split <- splitn(definition = def[[i_def]]$def,
                                     all_names = all_names,
                                     alt_names = alt_names,
                                     operators = comb_operators,
                                     scalar = TRUE)
      } else {
        def[[i_def]]$split <- splitn(definition = def[[i_def]]$def,
                                     all_names = all_names_avl,
                                     alt_names = alt_names_avl,
                                     operators = character(), 
                                     scalar = TRUE, dsplit = TRUE)
      }
    }
  }
  class(def) = class(features_def)
  return(def)
}

#' @title IFC_features_def Feature Redefinition
#' @description
#' Helper to rename a feature within IFC_features_def.
#' @param features_def an `IFC_features_def` object or a list containing features definition. Default is missing.
#' @param to_match_feat a string with a features_def name to use for matching 'features_def' names. Default is NULL
#' @param to_replace_feat a string of features_def name to use for replacing 'features_def' names. Default is NULL
#' @param force_default whether to force default names for features (except the one defined by 'to_replace_feat').
#' This removes custom names and replaces them with default values. Default is FALSE.
#' @param ... Other arguments to be passed.
#' @return an `IFC_features_def` object, or a list containing features definition
#' @keywords internal
redefine_features_def_feat <- function(features_def,
                                       to_match_feat = NULL,
                                       to_replace_feat = NULL,
                                       force_default = FALSE, ...) {
  dots = list(...)
  force_default=as.logical(force_default); assert(force_default, len=1, alw=c(TRUE,FALSE))
  if(length(to_match_feat) != 1 || length(to_replace_feat) != 1)  stop("'to_match_feat' should be of length 1")
  assert(to_match_feat, typ = "character")
  assert(to_replace_feat, typ = "character")
  def = features_def
  # if(length(unique(names(def))) != length(def)) stop("'features_def' should have unique names")
  old_names = sapply(features_def, FUN = function(x) x$name)
  tmp = which(old_names == to_match_feat)
  if(length(tmp) == 0) {
    warning("can't find [",to_match_feat,"] in features to rename")
    return(features_def)
  }
  if(length(tmp) != 1) {
    stop("[",to_match_feat,"] matches more then once in features to rename")
  }
  # rename desired feature
  def[[tmp]]$name <- to_replace_feat
  names(def)[tmp] <- to_replace_feat
  
  # determine which features are combined (which may be impacted by initial feature renaming)
  cmb_feats = sapply(def, USE.NAMES = FALSE, FUN = function(x) x$type == "combined")
  if(any(cmb_feats)) {
    # get indices of combined features
    cmb_feats = which(cmb_feats)
    comb_operators = c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT")
    # split combined features definition if not already done
    if(any(sapply(def[cmb_feats], FUN = function(def_def) length(def_def$split)) == 0)) {
      alt_names = gen_altnames(old_names)
      for(i_def in cmb_feats) {
        def[[i_def]]$split <- unname(splitn(definition = def[[i_def]]$def,
                                            all_names = old_names,
                                            alt_names = alt_names,
                                            operators = comb_operators,
                                            scalar = TRUE))
      }
    }
    
    # determine which features names can be modified, basically all when force_default is TRUE
    # or when FALSE only the ones whose names correspond to regular default name otherwise
    if(force_default) {
      to_modify = rep(TRUE, length(cmb_feats)) & (names(def)[cmb_feats] != to_replace_feat)
    } else {
      exp_names = sapply(def[cmb_feats], USE.NAMES = FALSE, feature_namer)
      to_modify = exp_names == names(def)[cmb_feats] 
    }
    
    # recursively modify combined features definition which depend on initial feature
    has_changed = TRUE
    while(has_changed) {
      has_changed = FALSE
      for(i_comb in seq_along(cmb_feats)) {
        i_def = cmb_feats[i_comb]
        for(i in seq_along(to_match_feat)) {
          tmp = def[[i_def]]$split == to_match_feat[i]
          if(any(tmp)) {
            def[[i_def]]$split[tmp] <- to_replace_feat[i]
          }
          def[[i_def]]$def = paste0(def[[i_def]]$split, collapse = "|")
          if(to_modify[i_comb]) {
            name = feature_namer(def[[i_def]])
            if(name != def[[i_def]]$name) {
              has_changed = TRUE
              def[[i_def]]$name = name
              names(def)[i_def] = name
            }
          }
        }
      }
    }
  }
  class(def) = class(features_def)
  return(def)
}

#' @title IFC_features_def Mask or Image Redefinition
#' @description
#' Helper to modify features_def according to masks or images redefinition.
#' @param features_def an `IFC_features_def` object, or a list containing features definition. Default is missing.
#' @param masks an `IFC_masks` object or a data.frame containing masks definition and name. Default is missing.
#' @param images a data.frame containing images definition. Default is missing.
#' @param force_default whether to force default names for masks and features.
#' This removes custom names and replaces them with default values. Default is FALSE.
#' @param ... Other arguments to be passed.
#' @return a list whose members are:\cr
#' -features_def, an `IFC_features_def` object, or a list containing features definition\cr
#' -masks, an `IFC_masks` object or a data.frame containing masks definition and name.\cr
#' -images, a data.frame containing images definition.
#' @keywords internal
redefine_features_def_msk_img <- function(features_def, masks, images, force_default = FALSE, ...) {
  dots = list(...)
  force_default=as.logical(force_default); assert(force_default, len=1, alw=c(TRUE,FALSE))
  
  if(force_default) { # recover default masks names
    default_names = sapply(masks, FUN = function(x) buildMask(masks, images, FALSE))
    to_modify = masks$name[masks$name != default_names[, "name"]]
    to_modify = setdiff(to_modify, c("MC","None","NMC"))
    to_replace = default_names[, "name"][masks$name %in% to_modify]
    obj_default = redefine_features_def_msk_img(features_def = features_def,
                                                masks = masks,
                                                images = images,
                                                to_match_mask = to_modify,
                                                to_replace_mask = to_replace,
                                                force_default = FALSE)
    features_def = obj_default$features_def
    masks = obj_default$masks
    images = obj_default$images
  }
  
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
  to_find = to_match_image
  to_replace_image = to_replace_image
  
  # ordering should ensure that names consisted of repeated pattern should be treated 1st
  order_ = order(nchar(to_find), decreasing = TRUE)
  to_find = to_find[order_]
  to_replace_image = to_replace_image[order_]
  
  def = features_def
  # if(length(unique(names(def))) != length(def)) stop("'features_def' should have unique names")
  old_names = names(features_def)
  comb_operators = c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT")
  all_msk = unlist(strsplit(as.character(masks$def[masks$name %in% "MC"]), split = "|Or|", fixed = TRUE), recursive = FALSE, use.names = FALSE)
  
  # we check if def has already been split otherwise we do it
  def = split_feat(features_def = def,
                   all_names = names(def),
                   m_names = c(all_msk, masks$name),
                   i_names = images$name, 
                   comb_operators = comb_operators,
                   extr_operators = c("true", "false", "True", "False"),
                   split = "|",
                   force = FALSE)
  
  # determine which features names can be modified, basically all when force_default is TRUE
  # or when FALSE only the ones whose names correspond to regular default name otherwise
  if(force_default) {
    to_modify = rep(TRUE, length(def))
  } else {
    exp_names = sapply(def, USE.NAMES = FALSE, feature_namer)
    to_modify = exp_names == names(def)
  }
  
  # recursively modify features definition and keep track of modified names in
  # to_match_feat-> to_replace_feat to allow further combined features modifications
  to_match_feat = c()
  to_replace_feat = c()
  has_changed = TRUE
  while(has_changed) {
    has_changed = FALSE
    for(i_def in seq_along(def)) {
      if(def[[i_def]]$type == "combined") {
        for(i in seq_along(to_match_feat)) {
          tmp = def[[i_def]]$split == to_match_feat[i]
          if(any(tmp)) {
            def[[i_def]]$split[tmp] <- to_replace_feat[i]
          }
          def[[i_def]]$def = paste0(def[[i_def]]$split, collapse = "|")
        }
      } else {
        for(i in seq_along(to_find)) {
          tmp = def[[i_def]]$split == to_find[i]
          if(any(tmp)) {
            def[[i_def]]$split[tmp] <- to_replace_image[i]
          }
        }
        def[[i_def]]$def <- paste0(def[[i_def]]$split, collapse = "|")
      }
      if(to_modify[i_def]) {
        name = feature_namer(def[[i_def]])
        if(name != def[[i_def]]$name) {
          has_changed = TRUE
          names(def)[i_def] = name
          to_match_feat =  c(to_match_feat, def[[i_def]]$name)
          to_replace_feat = c(to_replace_feat, name)
          tmp = duplicated(to_match_feat)
          if(length(to_replace_feat[tmp]) > 1) {
            stop("features redefinition results in multiple matches",
                 paste0("\n-\t", paste0(to_match_feat[tmp], " -> ", to_replace_feat[tmp])))
          } else {
            to_match_feat = to_match_feat[!tmp]
            to_replace_feat = to_replace_feat[!tmp]
          }
          def[[i_def]]$name = name
        }
      }
    }
  }
  class(def) <- class(features_def)
  return(list(features_def = def,
              masks = M,
              images = I))
}

#' @title IFC_features_def Redefinition
#' @description
#' Helper modify features_def according to masks or images redefinition.
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
#' Helper to modify images, masks, features, pops, graphs in a `IFC_data` object when features_def is changed.
#' @param obj an `IFC_data` object. Default is missing.
#' @param new_feat_def a list with new definitions from redefine_features_def() or redefine_features_def_msk_img().
#' @param ... Other arguments to be passed.
#' @return an `IFC_data` object.
#' @keywords internal
redefine_obj <- function(obj, new_feat_def, ...) {
  assert(obj, cla = "IFC_data")
  
  # retrieve mapping if any
  map = attr(new_feat_def, "map")
  if(length(map) == 0) map = list(initial = names(obj$features_def), to = names(new_feat_def$features_def))
  
  # modify features definition
  K = class(obj$features_def)
  obj$features_def <- sapply(USE.NAMES = TRUE, simplify = FALSE, new_feat_def$features_def, FUN = function(x) {
    x[c("name","type","userfeaturetype","def")] # remove split
  })
  
  # identify and remove duplicated features
  names(obj$features_def) = sapply(obj$features_def, FUN = function(x) x$name)
  dup = duplicated(names(obj$features_def))
  if(any(dup)) {
    warning("'obj' redefinition led to features with same name which have been discarded:\n\t-", paste0(unique(names(obj$features_def)[dup]), collapse = "\n\t-"))
  }
  map = sapply(map, simplify = FALSE, FUN = function(x) x[!dup])
  obj$features_def = obj$features_def[!dup]
  class(obj$features_def) = K
  
  # modify features values, at 1st create empty and then copy with mapping
  new_feat = as.data.frame(matrix(NaN, nrow = nrow(obj$features), ncol = length(obj$features_def)),
                           stringsAsFactors = FALSE, make.names = FALSE)
  colnames(new_feat) = names(obj$features_def)
  done = c()
  for(i in names(obj$features_def)) {
    tmp = map$to == i
    if(any(tmp)) {
      j = map$initial[tmp]
      tmp = names(obj$features) %in% j
      if(any(tmp)) {
        k = names(obj$features)[tmp]
        if(!any(k %in% done)) new_feat[, i] <- obj$features[, k[1]]
        done = c(done, k)
      }
    }
  }
  obj$features = structure(new_feat, class = class(obj$features))
  
  # modify graphs
  K = class(obj$graphs)
  obj$graphs = lapply(obj$graphs, FUN = function(g) {
    tmp = g$f1 == map$initial
    if(any(tmp)) {
      if(g$f1 == g$xlabel) g$xlabel <- map$to[tmp][1]
      g$f1 <- map$to[tmp][1]
    }
    if(g$type != "histogram") {
      tmp = g$f2 == map$initial
      if(any(tmp)) {
        if(g$f2 == g$ylabel) g$ylabel <- map$to[tmp][1]
        g$f2 <- map$to[tmp][1]
      }
    }
    return(g)
  })
  class(obj$graphs) <- K
  
  # modify pops
  K = class(obj$pops)
  obj$pops = lapply(obj$pops, FUN = function(p) {
    if("fx" %in% names(p)) {
      tmp = p$fx == map$initial
      if(any(tmp)) p$fx <- map$to[tmp][1]
    }
    if("fy" %in% names(p)) {
      tmp = p$fy == map$initial
      if(any(tmp)) p$fy <- map$to[tmp][1]
    }
    return(p)
  })
  class(obj$pops) <- K
  
  # modify masks
  obj$description$masks <- new_feat_def$masks
  
  # modify images
  obj$description$Images <- new_feat_def$images
  return(obj)
}

#' @title IFC_data Default Naming
#' @description
#' Helper to reset masks and features names to their default values in a `IFC_data` object.
#' @param obj an `IFC_data` object. Default is missing.
#' @param ... Other arguments to be passed.
#' @return an `IFC_data` object.
#' @keywords internal
usedefault_obj <- function(obj, ...) {
  assert(obj, cla = "IFC_data")
  new_feat_def =  redefine_features_def_msk_img(features_def = obj$features_def,
                                                masks = obj$description$masks,
                                                images = obj$description$Images,
                                                force_default = TRUE)
  return(redefine_obj(obj = obj, new_feat_def = new_feat_def))
}

#' @title Channel Switch
#' @description
#' Switches Channel in `IFC_data` object
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param from,to an integer index of channel. 'from' and 'to' should be different.
#' @param BF should 'from' channel be considered as brightfield. Default is TRUE.
#' @param MODE collection mode (as retrieved by getInfo) determining the range. Default is 1.
#' @details 'BF' and 'MODE' will be used only if 'to' is not found in 'obj'.\cr
#' If switching channel results in duplicated features definition, e.g. "Intensity_M01_Ch01" and "Intensity_M04_Ch04" exists in 'from'
#' and user called switch_channel(obj, 4, 1). So, "Intensity_M01_Ch01" will become "Intensity_M01_Ch01" (the same) and,
#' "Intensity_M04_Ch04" will become "Intensity_M01_Ch01". Meaning that resulting "Intensity_M01_Ch01" can originate from
#' either "Intensity_M01_Ch01" or "Intensity_M04_Ch04". In such case duplicates will be collected in 'dup' member of attr(, "map").\cr
#' /!\ Note also that 'initial' member of attr(, "map") will always be the one of 'from', i.e. "Intensity_M01_Ch01" will be mapped with "Intensity_M04_Ch04" .
#' @return a list, intended to be passed to 'new_feat_def' argument of getFromNamespace("redefine_obj", "IFC") whose members are:\cr
#' -features_def, an 'IFC_features_def' object, or a list containing features definition\cr
#' -masks, an 'IFC_masks' object or a data.frame containing masks definition and name.\cr
#' -images, a data.frame containing images definition.
#' @keywords internal
switch_channel <- function(obj, from, to, BF = TRUE, MODE = 1) {
  # check that from and to can be found in obj
  from = as.integer(from); assert(from, len = 1, alw = 1:12) #obj$description$Images$physicalChannel)
  to = as.integer(to); assert(to, len = 1, alw = setdiff(1:12, from))
  from_msk = sprintf("M%02i", from)
  
  # set mask and image defaults for from
  if(any(obj$description$Images$physicalChannel == from)) {
    from_img = obj$description$Images[obj$description$Images$physicalChannel == from, "name"]
  } else {
    warning("'from'[",sprintf("%02i",from),"] does not exist in 'obj'. Be careful when using the result.", immediate. = TRUE)
    BF = as.logical(BF); assert(BF, len = 1, alw = c(TRUE, FALSE))
    MODE = as.integer(MODE); assert(MODE, len = 1)
    from_img = sprintf("Ch%02i", from)
    obj$description$masks[obj$description$masks$name == "MC", "def"] = paste0(sprintf("M%02i", sort(c(from,obj$description$Images$physicalChannel))), collapse = "|Or|")
    def_img = rbind(buildImage(physicalChannel = from, name = from_img, BF = BF, MODE = MODE))
    def_img = def_img[,colnames(def_img) %in% colnames(obj$description$Images)]
    obj$description$Images = rbind(obj$description$Images, def_img)
    obj$description$Images = obj$description$Images[order(obj$description$Images$physicalChannel), , drop = FALSE]
  }
  
  # set mask and image defaults for to
  to_msk = sprintf("M%02i", to)
  if(any(obj$description$Images$physicalChannel == to)) {
    to_img = obj$description$Images[obj$description$Images$physicalChannel == to, "name"]
  } else {
    warning("'to'[",sprintf("%02i",to),"] does not exist in 'obj'. Be careful when using the result.", immediate. = TRUE)
    BF = as.logical(BF); assert(BF, len = 1, alw = c(TRUE, FALSE))
    MODE = as.integer(MODE); assert(MODE, len = 1)
    to_img = sprintf("Ch%02i", to)
    obj$description$masks[obj$description$masks$name == "MC", "def"] = paste0(sprintf("M%02i", sort(c(to,obj$description$Images$physicalChannel))), collapse = "|Or|")
    def_img = rbind(buildImage(physicalChannel = to, name = to_img, BF = BF, MODE = MODE))
    def_img = def_img[,colnames(def_img) %in% colnames(obj$description$Images)]
    obj$description$Images = rbind(obj$description$Images, def_img)
    obj$description$Images = obj$description$Images[order(obj$description$Images$physicalChannel), , drop = FALSE]
  }
  
  obj_default = list(masks = obj$description$masks,
                     images = obj$description$Images,
                     features_def = obj$features_def)
  
  comb_operators = c("+", "-", "*", "/", "(", ")", "ABS", "COS", "SIN", "SQR", "SQRT")
  all_msk = unlist(strsplit(as.character(obj_default$masks$def[obj_default$masks$name %in% "MC"]), split = "|Or|", fixed = TRUE), recursive = FALSE, use.names = FALSE)
  
  # create split for obj_default$features_def
  obj_default$features_def = split_feat(features_def = obj_default$features_def,
                                        m_names = c(all_msk, obj_default$masks$name),
                                        i_names = obj_default$images$name, 
                                        comb_operators = comb_operators,
                                        extr_operators = c("true", "false", "True", "False"),
                                        split = "|",
                                        force = FALSE)
  
  # add an extra image and mask to allow mask and image modifications
  obj_default$masks = rbind.data.frame(obj_default$masks[1, ], obj_default$masks, stringsAsFactors = TRUE)
  obj_default$images = rbind.data.frame(obj_default$images[1, ],obj_default$images, stringsAsFactors = FALSE)
  obj_default$images[1, "name"] = gen_altnames("foo", forbidden = obj$description$Images$name)
  obj_default$images[1, "physicalChannel"] = 0
  
  # use default names
  obj_default = redefine_features_def_msk_img(obj_default$features_def,
                                              masks = obj_default$masks, 
                                              images = obj_default$images,
                                              force_default = TRUE)
  
  # replace all 'to' masks and images by 'from' into obj_default
  new_images_names = obj_default$images$name
  new_images_names[new_images_names == to_img] <- from_img
  obj_default$masks[1, "name"] = to_msk
  obj_default$masks[1, "def"] = from_msk
  
  obj_from = redefine_features_def_msk_img(obj_default$features_def,
                                           masks = obj_default$masks, 
                                           images = obj_default$images,
                                           to_match_mask = to_msk,
                                           to_replace_mask = from_msk,
                                           new_images_names = new_images_names,
                                           force_default = FALSE)
  
  # replace all 'from' masks and images by 'to' into obj_default
  new_images_names = obj_default$images$name
  new_images_names[new_images_names == from_img] <- to_img
  obj_default$masks[1, "name"] = from_msk
  obj_default$masks[1, "def"] = to_msk
  
  obj_to = redefine_features_def_msk_img(obj_default$features_def,
                                         masks = obj_default$masks, 
                                         images = obj_default$images,
                                         to_match_mask = from_msk,
                                         to_replace_mask = to_msk,
                                         new_images_names = new_images_names,
                                         force_default = FALSE)
  
  # compute all masks names that have changed
  msk_names_from = buildMask(obj_from$masks, obj_from$images, definition = FALSE)
  msk_names_to = buildMask(obj_to$masks, obj_to$images, definition = FALSE)
  
  # compute final obj
  obj_final = redefine_features_def_msk_img(obj_to$features_def,
                                            masks = obj_to$masks, 
                                            images = obj_to$images,
                                            to_match_mask = msk_names_from,
                                            to_replace_mask = msk_names_to,
                                            new_images_names = new_images_names,
                                            force_default = FALSE)
  
  # remove duplicated images
  obj_final$images[obj_final$images$physicalChannel == from, "physicalChannel"] <- to
  obj_final$images = obj_final$images[-1, ]
  obj_final$images[obj$description$Images$physicalChannel == to, "physicalChannel"] <- from
  obj_final$images[obj$description$Images$physicalChannel == to, "name"] <- from_img
  obj_final$images = obj_final$images[order(as.integer(obj_final$images$physicalChannel)), ]
  
  # recover default mask
  obj_final$masks = obj_final$masks[-1, ]
  obj_final$masks[obj_final$masks$name %in% c("MC","None","NMC"), ] <- obj$description$masks[obj$description$masks$name %in% c("MC","None","NMC"), ]
  
  # remove duplicated features names
  no_dup = !duplicated(names(obj_final$features_def))
  dup = names(obj_final$features_def)[!no_dup]
  names(dup) = names(obj$features_def)[!no_dup]
  dup_i = names(obj$features_def)[!no_dup]         # initial names, including custom
  dup_o = names(obj_default$features_def)[!no_dup] # initial names, using default naming (i.e. without custom)
  dup_f = names(obj_from$features_def)[!no_dup]
  dup_t = names(obj_final$features_def)[!no_dup]
  tmp = which(dup_o == dup_f)
  
  # create mapping between input and output names
  map = list(initial = names(obj$features_def)[no_dup],
             to = names(obj_final$features_def)[no_dup],
             dup = structure(dup_t[tmp], names = dup_i[tmp]))
  # replace initial names by their duplicated value
  for(i in tmp) map$initial[map$initial == dup_t[i]] <- dup_i[i]
  obj_final$features_def = obj_final$features_def[no_dup]
  
  # keep only necessary masks
  S = sapply(obj_final$features_def, FUN = function(x) x$split)
  m_used = names(which(sapply(obj_final$mask$name, FUN = function(x) any(sapply(S, FUN = function(y) x %in% y)))))
  m_deps = sapply(obj_final$mask$name, FUN = function(x) which(sapply(obj_final$mask$def, FUN = function(y) x %in% strsplit(y, "|", fixed = TRUE)[[1]])))
  m_need = which(obj_final$mask$name %in% c("MC","None","NMC", m_used))
  LL = -1
  while(length(m_need) != LL) {
    LL = length(m_need)
    m_need = sort(unique(c(m_need, which(sapply(m_deps, FUN = function(x) any(sapply(m_need, FUN = function(y) y %in% x)))))))
  }
  obj_final$masks = obj_final$masks[m_need, ]
  obj_final$masks = obj_final$masks[!duplicated(obj_final$masks$name), ]
  
  attr(x = obj_final, which = "map") <- map
  return(obj_final)
}

#' @title Channel Swap
#' @description
#' Swaps Channels within `IFC_data` object
#' @param obj an `IFC_data` object extracted with features extracted.
#' @param chan1,chan2 an integer index of channel. 'chan1' and 'chan2' should be different.
#' @param BF should 'from' channel be considered as brightfield. Default is TRUE.
#' @param MODE collection mode (as retrieved by getInfo) determining the range. Default is 1.
#' @details 'BF' and 'MODE' will be used only if 'chan1' or 'chan2' is not found in 'obj'.\cr
#' /!\ NOTE: In case of conflict between resulting names, 'chan2' will be preferred, maximizing number of features with 'chan2'.
#' @return a list, intended to be passed to 'new_feat_def' argument of getFromNamespace("redefine_obj", "IFC") whose members are:\cr
#' -features_def, an 'IFC_features_def' object, or a list containing features definition\cr
#' -masks, an 'IFC_masks' object or a data.frame containing masks definition and name.\cr
#' -images, a data.frame containing images definition.
#' @keywords internal
swap_channel <- function(obj, chan1, chan2, BF = TRUE, MODE = 1) {
  # compute the results of changing channel name chan1 -> chan2
  obj_FT = switch_channel(obj, chan1, chan2, BF, MODE)
  # compute the results of changing channel name chan2 -> chan1
  obj_TF = switch_channel(obj, chan2, chan1, BF, MODE)
  
  # determine features names mapping
  map_TF = attr(obj_TF, "map")
  map_FT = attr(obj_FT, "map")
  
  # names which remains identical in both version chan1 -> chan2 and chan2 -> chan1
  keep0 = intersect(map_TF$to, map_FT$to)
  names(keep0) = sapply(keep0, FUN = function(x) map_FT$initial[map_FT$to == x])
  
  # all different except both for chan1 -> chan2 and chan2 -> chan1
  tmp1 = (map_FT$to != map_FT$initial) & !(map_FT$to %in% keep0)
  keep1 = map_FT$to[tmp1]
  names(keep1) = map_FT$initial[tmp1]
  tmp2 = (map_TF$to != map_TF$initial) & !(map_FT$to %in% keep0)
  keep2 = map_TF$to[tmp2]
  names(keep2) = map_TF$initial[tmp2]
  
  # conflict between names chan1 -> chan2 and chan2 -> chan1
  conflict = intersect(names(keep1), names(keep2))
  conflict = map_FT$to[map_FT$initial %in% conflict]
  names(conflict) = sapply(conflict, FUN = function(x) map_FT$initial[map_FT$to == x])
  keep1 = keep1[!names(keep1) %in% names(conflict)]
  keep2 = keep2[!names(keep2) %in% names(conflict)]
  
  # build final object
  obj_final = list(features_def = list(),
                   masks = data.frame(),
                   images = data.frame())
  obj_final$masks = rbind(obj_TF$masks, obj_FT$masks)
  obj_final$masks = obj_final$masks[!duplicated(obj_final$masks$name), ]
  
  obj_final$images = obj_FT$images
  obj_final$features_def = c(obj_FT$features_def[keep0],
                             obj_FT$features_def[keep1],
                             obj_TF$features_def[keep2],
                             obj_FT$features_def[conflict]) # we keep chan2 in case of conflict
  keep = c(keep0, keep1, keep2, conflict)
  
  # keep only necessary masks
  S = sapply(obj_final$features_def, FUN = function(x) x$split)
  m_used = names(which(sapply(obj_final$mask$name, FUN = function(x) any(sapply(S, FUN = function(y) x %in% y)))))
  m_deps = sapply(obj_final$mask$name, FUN = function(x) which(sapply(obj_final$mask$def, FUN = function(y) x %in% strsplit(y, "|", fixed = TRUE)[[1]])))
  m_need = which(obj_final$mask$name %in% c("MC","None","NMC", m_used))
  LL = -1
  while(length(m_need) != LL) {
    LL = length(m_need)
    m_need = sort(unique(c(m_need, which(sapply(m_deps, FUN = function(x) any(sapply(m_need, FUN = function(y) y %in% x)))))))
  }
  obj_final$masks = obj_final$masks[m_need, ]
  obj_final$masks = obj_final$masks[!duplicated(obj_final$masks$name), ]
  
  # determine new mapping
  map = list(initial = names(keep),
             to = unname(keep))
  
  attr(x = obj_final, which = "map") <- map
  return(obj_final)
}
