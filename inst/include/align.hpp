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
using namespace Rcpp;

//' @title Spatial Offsets Image Correction for Image
//' @name cpp_align_img
//' @description
//' This function uses bilinear interpolation to apply spatial offset correction on image
//' @param mat, a NumericMatrix.
//' @param dx, a double x spatial offset. It has to be within ]-1,+1[. Default is NA_REAL for no change.
//' @param dy, a double y spatial offset. It has to be within ]-1,+1[. Default is NA_REAL for no change.
//' @details It is intended to be applied on raw images matrices from .rif files so has to generate spatial offset corrected image matrices.\cr
//' See William E. Ortyn et al. Sensitivity Measurement and Compensation in Spectral Imaging. Cytometry A 69 852-862 (2006).
//' \doi{10.1002/cyto.a.20306}
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_align_img(const Rcpp::NumericMatrix mat,
                                  const double dx = NA_REAL,
                                  const double dy = NA_REAL) {
  if(any(is_na(Rcpp::NumericVector::create(dx, dy)))) {
    if(all(is_na(Rcpp::NumericVector::create(dx, dy)))) {
      return(Rcpp::clone(mat));
    }
    Rcpp::stop("hpp_align_img: bad offset value");
  } 
  if((std::abs(dx) >= 1.0) || (std::abs(dy) >= 1.0)) Rcpp::stop("hpp_align_img: offset should be ]-1,+1[");
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  Rcpp::NumericMatrix out = Rcpp::no_init_matrix(mat_r, mat_c);
  if((dx == 0.0) && (dy == 0.0)) {
    for(R_len_t i_col = 1; i_col < mat_c - 1; i_col++) {
      for(R_len_t i_row = 2; i_row < mat_r - 1; i_row++) {
        out(i_row, i_col) = mat(i_row, i_col);
      }
    }
    return out(Rcpp::Range(2, mat_r - 2), Rcpp::Range(1, mat_c  - 2)); 
  } else {
    double sx = dx, sy = dy, ssx, ssy;
    uint8_t deltax = 0, deltay = 0;
    if(dx < 0) {
      sx = 1 - std::abs(dx);
      deltax = 1;
    }
    if(dy < 0) {
      sy = 1 - std::abs(dy);
      deltay = 1;
    }
    ssx = 1 - sx;
    ssy = 1 - sy;
    for(R_len_t i_col = 0; i_col < mat_c - 1; i_col++) {
      for(R_len_t i_row = 1; i_row < mat_r - 1; i_row++) {
        out(i_row, i_col) = (mat(i_row    , i_col) * ssx + mat(i_row    , i_col + 1) * sx) * ssy +
                            (mat(i_row + 1, i_col) * ssx + mat(i_row + 1, i_col + 1) * sx) * sy;
      }
    }
    return out(Rcpp::Range(2 - deltay, mat_r - 2 - deltay), Rcpp::Range(1 - deltax, mat_c  - 2 - deltax));
  }
  return R_NilValue;
}

//' @title Spatial Offsets Image Correction for Mask
//' @name cpp_align_msk
//' @description
//' This function uses bilinear interpolation to apply spatial offset correction on mask
//' @param msk, a IntegerMatrix.
//' @param dx, a double x spatial offset. It has to be within ]-1,+1[. Default is NA_REAL for no change.
//' @param dy, a double y spatial offset. It has to be within ]-1,+1[. Default is NA_REAL for no change.
//' @details It is intended to be applied on raw images matrices from .rif files so has to generate spatial offset corrected image matrices.\cr
//' See William E. Ortyn et al. Sensitivity Measurement and Compensation in Spectral Imaging. Cytometry A 69 852-862 (2006).
//' \doi{10.1002/cyto.a.20306}
//' @return a IntegerMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix hpp_align_msk(const Rcpp::IntegerMatrix msk,
                                  const double dx = NA_REAL,
                                  const double dy = NA_REAL) {
  if(any(is_na(Rcpp::NumericVector::create(dx, dy)))) {
    if(all(is_na(Rcpp::NumericVector::create(dx, dy)))) {
      return(Rcpp::clone(msk));
    }
    Rcpp::stop("hpp_align_msk: bad offset value");
  } 
  if((std::abs(dx) >= 1.0) || (std::abs(dy) >= 1.0)) Rcpp::stop("hpp_align_msk: offset should be ]-1,+1[");
  R_len_t mat_r = msk.nrow();
  R_len_t mat_c = msk.ncol();
  Rcpp::IntegerMatrix out = Rcpp::no_init_matrix(mat_r, mat_c);
  if((dx == 0.0) && (dy == 0.0)) {
    for(R_len_t i_col = 1; i_col < mat_c - 1; i_col++) {
      for(R_len_t i_row = 2; i_row < mat_r - 1; i_row++) {
        out(i_row, i_col) = msk(i_row, i_col);
      }
    }
    return out(Rcpp::Range(2, mat_r - 2), Rcpp::Range(1, mat_c  - 2));
  } else{
    double sx = dx, sy = dy, ssx, ssy;
    uint8_t deltax = 0, deltay = 0;
    if(dx < 0) {
      sx = 1 - std::abs(dx);
      deltax = 1;
    }
    if(dy < 0) {
      sy = 1 - std::abs(dy);
      deltay = 1;
    }
    ssx = 1 - sx;
    ssy = 1 - sy;
    for(R_len_t i_col = 0; i_col < mat_c - 1; i_col++) {
      for(R_len_t i_row = 1; i_row < mat_r - 1; i_row++) {
        out(i_row, i_col) = (msk(i_row    , i_col) * ssx + msk(i_row    , i_col + 1) * sx) * ssy +
                            (msk(i_row + 1, i_col) * ssx + msk(i_row + 1, i_col + 1) * sx) * sy;
      }
    }
    return out(Rcpp::Range(2 - deltay, mat_r - 2 - deltay), Rcpp::Range(1 - deltax, mat_c  - 2 - deltax));
  }
  return R_NilValue;
}

//' @title Spatial Offsets Image Correction
//' @name cpp_align
//' @description
//' This function uses bilinear interpolation to apply spatial offset correction on image
//' @param mat, a NumericMatrix.
//' @param dx, a double x spatial offset. It has to be within ]-1,+1[. Default is NA_REAL for no change.
//' @param dy, a double y spatial offset. It has to be within ]-1,+1[. Default is NA_REAL for no change.
//' @details It is intended to be applied on raw images matrices from .rif files so has to generate spatial offset corrected image matrices.\cr
//' See William E. Ortyn et al. Sensitivity Measurement and Compensation in Spectral Imaging. Cytometry A 69 852-862 (2006).
//' \doi{10.1002/cyto.a.20306}
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix hpp_align(const Rcpp::NumericMatrix mat,
                              const double dx = NA_REAL,
                              const double dy = NA_REAL) {
  Rcpp::NumericMatrix out = hpp_align_img(mat, dx, dy);
  if(mat.hasAttribute("mask")) {
    IntegerMatrix foo = mat.attr("mask");
    Rcpp::IntegerMatrix MM = hpp_align_msk(foo, dx, dy);
    if(foo.hasAttribute("removal")) MM.attr("removal") = foo.attr("removal");
    out.attr("mask") = MM;
  }
  return out;
}

#endif
