################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2022 Yohann Demont                                             #
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

#' @title Population Sampling
#' @description
#' Creates a sample from a population
#' @param obj an `IFC_data` object extracted by ExtractFromDAF(extract_features = TRUE) or ExtractFromXIF(extract_features = TRUE).
#' @param pop name of the population to sample.
#' @param size a non-negative integer giving the number of items to choose.
#' @param new_name name of the exported population.
#' @param random_seed a list of elements to pass to \link[base]{set.seed} or a single value, interpreted as an integer, or NULL to be used when 'add_noise' is set to TRUE. Default is NULL.
#' Note that NA_integer_ or list(seed = NA_integer_) can be used to not call \link[base]{set.seed} at all.
#' @param ... Other arguments to be passed.
#' @details population is exported as tagged population.
#' @return an IFC_data object with sampled pop added.
#' @keywords internal
data_add_pop_sample = function(obj, pop, size, new_name, random_seed = NULL, ...) {
  dots = list(...)
  assert(obj, cla = "IFC_data")
  obj_number = obj$description$ID$objcount
  size = as.integer(size[!is.na(size) & (size >= 0)]); assert(size, len = 1)
  pop = as.character(pop); assert(pop, len = 1)
  if(!any(pop == names(obj$pops))) stop("can't find 'pop': \"",pop,"\" in 'obj$pops'")
  if(missing(new_name)) stop("'new_name' can't be missing")
  new_name = as.character(new_name); assert(new_name, len = 1, typ = "character")
  if(any(new_name == names(obj$pops))) stop("'new_name': \"",new_name,"\" should be different from already existing names of 'obj$pops'")
  SEED = fetch_seed(random_seed)
  f = function(x) { 
    with_seed(x, SEED$seed, SEED$kind, SEED$normal.kind, SEED$sample.kind)
  }
  idx = which(obj$pops[[pop]]$obj) - 1
  if(length(idx) < size) {
    warning("desired sample 'size' is larger than the population")
    foo = idx
  } else {
    foo = integer()
    if(length(idx) != 0) {
      foo = f(sample(x = idx, size = size, replace = FALSE))
    }
  }
  data_add_pops(pops = list(list(name = new_name, type = "T",
                                 color = obj$pops[[pop]]$color, lightModeColor = obj$pops[[pop]]$lightModeColor,
                                 obj = foo)),
                obj = obj, ...)
}
