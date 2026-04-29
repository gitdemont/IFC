/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2026 Yohann Demont                                              

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
  
#ifndef IFC_AFFINE_HPP
#define IFC_AFFINE_HPP

#include <Rcpp.h>  
#include "resize.hpp"
using namespace Rcpp;

template <typename type>
double cubic   (type *m, type *n, type *o, type *p, double *s, double a) {
  //  |   0     1       0    0 |
  //  |   a     0      -a    0 |
  //  |-2*a  -a-3   2*a+3    a |
  //  |   a   a+2    -a-2   -a |
  double u = *m * (   a) + *n * ( 2+a) + *o * ( -2-a) + *p * (-a);
  double v = *m * (-2*a) + *n * (-3-a) + *o * (3+2*a) + *p * ( a);
//double w = *m * (   a) + *n * (   0) + *o * (   -a) + *p * ( 0);
  double w = *m * (   a) +             + *o * (   -a)            ;
//double z = *m * (   0) + *n * (   1) + *o * (    0) + *p * ( 0);
//return u * std::pow(*s,3) + v * std::pow(*s,2) + w * std::pow(*s,1) + z * std::pow(*s,0);
  return u * *s * *s * *s + v * *s * *s + w * *s + *n;
}
template <typename type>
double cubic05 (type *m, type *n, type *o, type *p, double *s) { // Catmull-Rom, a = -0.5
  return *n + 0.5 * *s * (*o - *m + *s * (2.0 * *m - 5.0 * *n + 4.0 * *o - *p + *s * (3.0 * (*n - *o) + *p - *m)));
}

template <int RTYPE>  
Rcpp::Matrix<RTYPE> affine_T(Rcpp::Matrix<RTYPE> mat,
                             const Rcpp::NumericMatrix tr,
                             const Rcpp::IntegerVector size = Rcpp::IntegerVector::create(0,0),
                             const uint8_t interpolation = 2,
                             const bool add_noise = true, 
                             const double bg = 0.0,
                             const double sd = 0.0) {
  if(tr.nrow() != 2 || tr.ncol() != 3) Rcpp::stop("hpp_affine: 'tr' should be a length [2,3] matrix");
  if(size.size() != 2) Rcpp::stop("hpp_affine: 'size' should be a length 2 vector");
  R_len_t mat_r = mat.nrow(), mat_c = mat.ncol();
  
  // define output size
  R_len_t new_height = size[0] == 0 ? mat_r : size[0];
  R_len_t new_width =  size[1] == 0 ? mat_c : size[1];
  
  // instantiate out
  Rcpp::Matrix<RTYPE> out = Rcpp::no_init_matrix(new_height, new_width);
  
  // initialise bkg with bg or noisy bg
  Rcpp::Vector<RTYPE> bkg = Rcpp::no_init_vector(new_height * new_width);
  if(add_noise) {
    bkg = Rcpp::rnorm(new_height * new_width, bg, sd); 
  } else {
    bkg.fill(bg);
  }
  
  switch(interpolation) {
    case 2: { // bilinear
      R_len_t mm1_r = mat_r - 1, mm1_c = mat_c - 1;
      for(R_len_t i_col = 0, i = 0; i_col < new_width; i_col++) {
        for(R_len_t i_row = 0; i_row < new_height; i_row++, i++) {
          double x = tr[0]*(i_col+0.5) + tr[2]*(i_row+0.5) + tr[4];
          double y = tr[1]*(i_col+0.5) + tr[3]*(i_row+0.5) + tr[5];
          x -= 0.5;
          y -= 0.5;
          R_len_t fy = std::floor(y);
          R_len_t fx = std::floor(x);
          double sx = x-fx;
          double sy = y-fy;
          out[i] = (((fy< 0 || fy>=mat_r || fx< 0 || fx>=mat_c) ? bkg[i] : mat(fy  , fx  )) * (1-sx) +
                    ((fy< 0 || fy>=mat_r || fx<-1 || fx>=mm1_c) ? bkg[i] : mat(fy  , fx+1)) *    sx) * (1-sy) +
                   (((fy<-1 || fy>=mm1_r || fx< 0 || fx>=mat_c) ? bkg[i] : mat(fy+1, fx  )) * (1-sx) +
                    ((fy<-1 || fy>=mm1_r || fx<-1 || fx>=mm1_c) ? bkg[i] : mat(fy+1, fx+1)) *    sx) *    sy;
        }
      }
    }
      break;
  case 3: { // bicubic, a = -0.5
      Rcpp::Vector<RTYPE> u = Rcpp::no_init_vector(4);
      Rcpp::NumericVector v = Rcpp::no_init_vector(4);
      for(R_len_t i_col = 0, i = 0; i_col < new_width; i_col++) {
        for(R_len_t i_row = 0; i_row < new_height; i_row++, i++) {
          double x = tr[0]*(i_col+0.5) + tr[2]*(i_row+0.5) + tr[4];
          double y = tr[1]*(i_col+0.5) + tr[3]*(i_row+0.5) + tr[5];
          x -= 0.5;
          y -= 0.5;
          R_len_t fx = std::floor(x);
          R_len_t fy = std::floor(y);
          double sx = x-fx;
          double sy = y-fy;
          fx--;
          for(short m = 0; m < 4; m++) {
            bool x_in = fx>=0 && fx<mat_c;
            R_len_t yy = std::floor(y) - 1;
            for(short n = 0; n < 4; n++) {
              u[n] = (x_in && yy>=0 && yy<mat_r) ? mat(yy, fx) : bkg[i];
              yy++;
            }
            v[m] = cubic05(&u[0], &u[1], &u[2], &u[3], &sy);
            fx++;
          }
          out[i] = cubic05(&v[0], &v[1], &v[2], &v[3], &sx);
        }
      }
    }
      break;
  // case 4: { // bicubic, a = -0.75, TODO pass 'a' param
  //     Rcpp::Vector<RTYPE> u = Rcpp::no_init_vector(4);
  //     Rcpp::NumericVector v = Rcpp::no_init_vector(4);
  //     for(R_len_t i_col = 0, i = 0; i_col < new_width; i_col++) {
  //       for(R_len_t i_row = 0; i_row < new_height; i_row++, i++) {
  //         double x = tr[0]*(i_col+0.5) + tr[2]*(i_row+0.5) + tr[4];
  //         double y = tr[1]*(i_col+0.5) + tr[3]*(i_row+0.5) + tr[5];
  //         x -= 0.5;
  //         y -= 0.5;
  //         R_len_t fx = std::floor(x);
  //         R_len_t fy = std::floor(y);
  //         double sx = x-fx;
  //         double sy = y-fy;
  //         fx--;
  //         for(short m = 0; m < 4; m++) {
  //           bool x_in = fx>=0 && fx<mat_c;
  //           R_len_t yy = std::floor(y) - 1;
  //           for(short n = 0; n < 4; n++) {
  //             u[n] = (x_in && yy>=0 && yy<mat_r) ? mat(yy, fx) : bkg[i];
  //             yy++;
  //           }
  //           v[m] = cubic(&u[0], &u[1], &u[2], &u[3], &sy, 0.75);
  //           fx++;
  //         }
  //         out[i] = cubic(&v[0], &v[1], &v[2], &v[3], &sx, 0.75);
  //       }
  //     }
  //   }
  //     break;
    default: { // Default Nearest 
      for(R_len_t i_col = 0, i = 0; i_col < new_width; i_col++) {
        for(R_len_t i_row = 0; i_row < new_height; i_row++, i++) {
          R_len_t fx = std::floor(tr[0]*(i_col+0.5) + tr[2]*(i_row+0.5) + tr[4]);
          R_len_t fy = std::floor(tr[1]*(i_col+0.5) + tr[3]*(i_row+0.5) + tr[5]);
          out[i] = (fy>=0 && fy<mat_r && fx>=0 && fx<mat_c) ? mat(fy, fx) : bkg[i];
        }
      }
    }
  }
  return resize_T(out, size[0], size[1], add_noise, bg, sd);
}

//' @title Affine Transformation
//' @name cpp_affine
//' @description
//' Function that applies affine transformation on 2D image
//' @param mat Matrix (Raw, Integer, Logical or Numeric).
//' @param tr NumericMatrix, the affine transformation matrix.
//' @param size IntegerVector, final dimension of resulted transformation. Default is c(0,0) for no change. Negative values will remove rows/cols from output.
//' @param interpolation uint8_t, type of interpolation applied on output (2 = "bilinear", 3 = "bicubic" a = -0.5). Default is 2. Values other than 2 or 3 will result in nearest neighbor interpolation.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions. Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.0.
//' @return the result of the affine transformation of mat.
//' @keywords internal
// [[Rcpp::export]]
SEXP hpp_affine (SEXP mat,
                 const Rcpp::NumericMatrix tr,
                 const Rcpp::IntegerVector size = Rcpp::IntegerVector::create(0,0),
                 const uint8_t interpolation = 2,
                 const bool add_noise = true, 
                 const double bg = 0.0,
                 const double sd = 0.0) {
  switch(TYPEOF(mat)) {
  // case NILSXP : return mat;
  case RAWSXP : return affine_T<RAWSXP>(mat, tr, size, interpolation, add_noise, bg, sd);
  case LGLSXP : return affine_T<LGLSXP>(mat, tr, size, interpolation, add_noise, bg, sd);
  case INTSXP : return affine_T<INTSXP>(mat, tr, size, interpolation, add_noise, bg, sd);
  case REALSXP : return affine_T<REALSXP>(mat, tr, size, interpolation, add_noise, bg, sd);
  default: Rcpp::stop("hpp_affine: not supported type in 'mat'");
  }
}

template <int RTYPE> 
SEXP affrescale_T (Rcpp::Matrix<RTYPE> mat,
                   const double scale = 1.0,
                   const bool resize = false,
                   const uint8_t interpolation = 2,
                   const bool add_noise = true,
                   const double bg = 0.0,
                   const double sd = 0.0) {
  if(!R_finite(scale)) Rcpp::stop("affrescale_T: 'scale' should be finite");
  if(scale <= 0.0)  Rcpp::stop("affrescale_T: 'scale' should be > 0.0");
  R_len_t mat_r = mat.nrow();
  R_len_t mat_c = mat.ncol();
  double k = 1 / scale;
  Rcpp::NumericMatrix tr = Rcpp::no_init_matrix(2, 3);
  tr[0] = k; tr[1] = 0.0; tr[2] = 0.0; tr[3] = k; tr[4] = 0.0; tr[5] = 0.0;
  Rcpp::IntegerVector size = Rcpp::IntegerVector::create(mat_r * scale, mat_c * scale);
  if(resize) return resize_T(affine_T(mat, tr, size, 2, add_noise, bg, sd), mat_r, mat_c, add_noise, bg, sd);
  return affine_T(mat, tr, size, 2, add_noise, bg, sd);
}
                 
//' @title Image Upscaling
//' @name cpp_upscale
//' @description
//' Function that upscales mat according to scale.
//' @param mat Matrix (Raw, Integer Logical or Numeric).
//' @param scale double, giving the scaling factor. Default is 1.0 for no change.
//' @param resize bool, whether to resize upscaled mat to its original dimension. Default is false.
//' @param interpolation uint8_t, type of interpolation applied on output (2 = "bilinear", 3 = "bicubic" a = -0.5). Default is 2. Values other than 2 or 3 will result in nearest neighbor interpolation.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions. Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.0.
//' @return a upscaled matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
SEXP hpp_upscale (SEXP mat,
                  const double scale = 1.0,
                  const bool resize = false,
                  const uint8_t interpolation = 2,
                  const bool add_noise = true,
                  const double bg = 0.0,
                  const double sd = 0.0) {
  switch(TYPEOF(mat)) {
  // case NILSXP : return mat;
  case RAWSXP : return affrescale_T<RAWSXP>(mat, scale, resize, interpolation, add_noise, bg, sd);
  case LGLSXP : return affrescale_T<LGLSXP>(mat, scale, resize, interpolation, add_noise, bg, sd);
  case INTSXP : return affrescale_T<INTSXP>(mat, scale, resize, interpolation, add_noise, bg, sd);
  case REALSXP : return affrescale_T<REALSXP>(mat, scale, resize, interpolation, add_noise, bg, sd);
  default: Rcpp::stop("hpp_upscale: not supported type in 'mat'");
  }
}

//' @title Image Downscaling
//' @name cpp_downscale
//' @description
//' Function that downscales mat according to scale.
//' @param mat Matrix (Raw, Integer Logical or Numeric).
//' @param scale double, giving the scaling factor. Default is 1.0 for no change.
//' @param resize bool, whether to resize downscaled mat to its original dimension. Default is false.
//' @param interpolation uint8_t, type of interpolation applied on output (2 = "bilinear", 3 = "bicubic" a = -0.5). Default is 2. Values other than 2 or 3 will result in nearest neighbor interpolation.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions. Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.0.
//' @return a downscaled matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
SEXP hpp_downscale (SEXP mat,
                    const double scale = 1.0,
                    const bool resize = false,
                    const uint8_t interpolation = 2,
                    const bool add_noise = true,
                    const double bg = 0.0,
                    const double sd = 0.0) {
  switch(TYPEOF(mat)) {
  // case NILSXP : return mat;
  case RAWSXP : return affrescale_T<RAWSXP>(mat, 1/scale, resize, interpolation, add_noise, bg, sd);
  case LGLSXP : return affrescale_T<LGLSXP>(mat, 1/scale, resize, interpolation, add_noise, bg, sd);
  case INTSXP : return affrescale_T<INTSXP>(mat, 1/scale, resize, interpolation, add_noise, bg, sd);
  case REALSXP : return affrescale_T<REALSXP>(mat, 1/scale, resize, interpolation, add_noise, bg, sd);
  default: Rcpp::stop("hpp_downscale: not supported type in 'mat'");
  }
}

#endif
