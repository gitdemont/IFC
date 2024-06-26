################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2020 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2020                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,             #
#                     Jean-Pierre Marolleau, Loïc Garçon,                      #
#                     INSERM, UPD, CHU Amiens                                  #
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

#' @title Object File Export
#' @description
#' Exports images to various types.
#' @param x a numeric matrix.
#' @param type image type. Supported values are: "bmp", "jpeg", "png", and "tiff".
#' @param ... other arguments to be passed.
#' @keywords internal
objectWrite <- function(x, type, ...) {
  dots=list(...)
  switch(type, 
         png = {
           text = c("software" = paste0("IFC ",asNamespace("IFC")[[".pkgenv"]][["version"]]))
           text = c(text, dots$text)
           text = c(text, c("title" = dots$tags$`269`$map, "creation time" = dots$tags$`306`$map))
           dots = dots[!(names(dots) %in% c("text", "tags"))]
           do.call(what = png::writePNG, args = c(list(image = x, text = text), dots))
         },
         tiff = do.call(what = writetiff, args = c(list(image = x), dots)),
         multi = do.call(what = writemulti, args = c(list(image = x), dots)),
         jpeg = {
           dots = dots[!(names(dots) %in% "tags")]
           do.call(what = jpeg::writeJPEG, args = c(list(image = x), dots))
         },
         bmp = {
           if((length(dots) != 0)) {
             switch(typeof(dots[[1]]),
                    raw = return(cpp_writeBMP(image = x)),
                    character = return(writeBin(object = cpp_writeBMP(image = x), con = dots[[1]])),
                    stop("objectWrite: target should be raw() or a character file path")
                    )
           }
           stop("objectWrite: no target to write \"bmp\" to")
         },
         stop(paste0("objectWrite: can only write to \"bmp\", \"jpeg\", \"png\" and \"tiff\"; [",type,"] is not supported")))
}
