################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2022 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2022                                                                  #
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

#' @title Groups Combination
#' @description 
#' Combines data.frame groups into a unique one
#' @param df a data.frame whose columns represent grouping factor.
#' @param collapse a string used to collapse groups levels.
#' @param keepNAlevels whether to keep NA levels resulting from groups merging
#' @param ... other arguments to be passed.
#' @return a named list containing row indices of grouping combinations.
#' @keywords internal
group <- function(df, collapse = ".", keepNAlevels = FALSE, ...) {
  assert(keepNAlevels, alw = c(TRUE, FALSE))
  assert(collapse, typ = "character", len = 1)
  assert(df, cla = "data.frame")
  ret = cpp_group_df(df)
  v = order(ret)
  x = attr(ret, "table")
  y = 1
  n = character()
  drop = TRUE # always drop unused levels
  if(drop) x = x[x != 0]
  ans = structure(lapply(seq_along(x), FUN = function(i) {
    if(x[i] == 0) {
      n <<- c(n, "")
      return(NULL)
    }
    z = v[y:(y+x[i]-1)]
    y <<- y + x[i]
    x[i] <<- y
    n <<- c(n, paste0(df[z[1],], collapse = collapse))
    if(!keepNAlevels && anyNA(df[z[1],])) return(NULL)
    z
  }), names = n)
  ans[sapply(ans, length) != 0]
}
