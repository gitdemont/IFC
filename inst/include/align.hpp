/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2021 Yohann Demont                                              

  It is part of IFC package, please cite:                                       
  -IFC: An R Package for Imaging Flow Cytometry              
  -YEAR: 2021                                                                   
  -COPYRIGHT HOLDERS: Yohann Demont, Jean-Pierre Marolleau, Loïc Garçon,        
                      CHU Amiens                                                

  DISCLAIMER:                                                                   
  -You are using this package on your own risk!                                 
  -We do not guarantee privacy nor confidentiality.                             
  -This program is distributed in the hope that it will be useful, but WITHOUT  
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or         
  FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or  
  contributors be liable for any direct, indirect, incidental, special,         
  exemplary, or consequential damages (including, but not limited to,           
  procurement of substitute goods or services; loss of use, data, or profits;   
  or business interruption) however caused and on any theory of liability,      
  whether in contract, strict liability, or tort (including negligence or       
  otherwise) arising in any way out of the use of this software, even if        
  advised of the possibility of such damage.                                    

  You should have received a copy of the GNU General Public License             
  along with IFC. If not, see <http://www.gnu.org/licenses/>.                   
*/
  
#ifndef IFC_ALIGN_HPP
#define IFC_ALIGN_HPP
  
#include <Rcpp.h>
#include "affine.hpp"
using namespace Rcpp;

//' @title Spatial Offsets Image Correction for Image
//' @name cpp_align_img
//' @description
//' This function uses bilinear interpolation to apply spatial offset correction on image
//' @param mat, a NumericMatrix.
//' @param dx, a double x spatial offset. Default is NA_REAL for no change.
//' @param dy, a double y spatial offset. Default is NA_REAL for no change.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions. Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.0.
//' @details It is intended to be applied on raw images matrices from .rif files so has to generate spatial offset corrected image matrices.\cr
//' See William E. Ortyn et al. Sensitivity Measurement and Compensation in Spectral Imaging. Cytometry A 69 852-862 (2006).
//' \doi{10.1002/cyto.a.20306}
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_align_img(const Rcpp::NumericMatrix mat,
                                  const double dx = NA_REAL,
                                  const double dy = NA_REAL,
                                  const bool add_noise = true, 
                                  const double bg = 0.0,
                                  const double sd = 0.0) {
  if(any(is_na(Rcpp::NumericVector::create(dx, dy)))) {
    if(all(is_na(Rcpp::NumericVector::create(dx, dy)))) {
      return(Rcpp::clone(mat));
    }
    Rcpp::stop("hpp_align_img: bad offset value");
  }
  if(dx == 0.0 && dy == 0.0) return(mat);
  // if((std::abs(dx) >= 1.0) || (std::abs(dy) >= 1.0)) Rcpp::stop("hpp_align_img: offset should be ]-1,+1[");
  Rcpp::NumericMatrix tr = Rcpp::no_init_matrix(2, 3);
  tr[0] = 1.0; tr[1] = 0.0; tr[2] = 0.0; tr[3] = 1.0; tr[4] = dx; tr[5] = dy;
  return affine_T(mat, tr, Rcpp::IntegerVector::create(0,0), 2, add_noise, bg, sd);
}

//' @title Spatial Offsets Image Correction for Mask
//' @name cpp_align_msk
//' @description
//' This function uses bilinear interpolation to apply spatial offset correction on mask
//' @param msk, a IntegerMatrix.
//' @param dx, a double x spatial offset. Default is NA_REAL for no change.
//' @param dy, a double y spatial offset. Default is NA_REAL for no change.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions. Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.0.
//' @details It is intended to be applied on raw images matrices from .rif files so has to generate spatial offset corrected image matrices.\cr
//' See William E. Ortyn et al. Sensitivity Measurement and Compensation in Spectral Imaging. Cytometry A 69 852-862 (2006).
//' \doi{10.1002/cyto.a.20306}
//' @return a IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix hpp_align_msk(const Rcpp::IntegerMatrix msk,
                                  const double dx = NA_REAL,
                                  const double dy = NA_REAL,
                                  const bool add_noise = true, 
                                  const double bg = 0.0,
                                  const double sd = 0.0) {
  if(any(is_na(Rcpp::NumericVector::create(dx, dy)))) {
    if(all(is_na(Rcpp::NumericVector::create(dx, dy)))) {
      return(Rcpp::clone(msk));
    }
    Rcpp::stop("hpp_align_msk: bad offset value");
  }
  if(dx == 0.0 && dy == 0.0) return(msk);
  // if((std::abs(dx) >= 1.0) || (std::abs(dy) >= 1.0)) Rcpp::stop("hpp_align_msk: offset should be ]-1,+1[");
  Rcpp::NumericMatrix tr = Rcpp::no_init_matrix(2, 3);
  tr[0] = 1.0; tr[1] = 0.0; tr[2] = 0.0; tr[3] = 1.0; tr[4] = dx; tr[5] = dy;
  Rcpp::IntegerMatrix foo = affine_T(msk, tr, Rcpp::IntegerVector::create(0,0), 2, add_noise, bg, sd);
  if(msk.hasAttribute("removal")) foo.attr("removal") = msk.attr("removal");
  return foo;
}

//' @title Spatial Offsets Image Correction
//' @name cpp_align
//' @description
//' This function uses bilinear interpolation to apply spatial offset correction on image
//' @param mat, a NumericMatrix.
//' @param dx, a double x spatial offset. Default is NA_REAL for no change.
//' @param dy, a double y spatial offset. Default is NA_REAL for no change.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions. Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.0.
//' @details It is intended to be applied on raw images matrices from .rif files so has to generate spatial offset corrected image matrices.\cr
//' See William E. Ortyn et al. Sensitivity Measurement and Compensation in Spectral Imaging. Cytometry A 69 852-862 (2006).
//' \doi{10.1002/cyto.a.20306}
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_align(const Rcpp::NumericMatrix mat,
                              const double dx = NA_REAL,
                              const double dy = NA_REAL,
                              const bool add_noise = true, 
                              const double bg = 0.0,
                              const double sd = 0.0) {
  Rcpp::NumericMatrix out = hpp_align_img(mat, dx, dy, add_noise, bg, sd);
  if(mat.hasAttribute("mask")) {
    IntegerMatrix foo = mat.attr("mask");
    Rcpp::IntegerMatrix MM = hpp_align_msk(foo, dx, dy, add_noise, bg, sd);
    if(foo.hasAttribute("removal")) MM.attr("removal") = foo.attr("removal");
    out.attr("mask") = MM;
  }
  return out;
}

#endif
