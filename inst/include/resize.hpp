/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2020 Yohann Demont                                              
                                                                                
  It is part of IFC package, please cite:                                       
  -IFC: An R Package for Imaging Flow Cytometry                                 
  -YEAR: 2020                                                                   
  -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,              
                      Jean-Pierre Marolleau, Loïc Gaççon,                       
                      INSERM, UPD, CHU Amiens                                   
                                                                                
                                                                                
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

#ifndef IFC_RESIZE_HPP
#define IFC_RESIZE_HPP

#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>  
Rcpp::Matrix<RTYPE> crop_T(Rcpp::Matrix<RTYPE> mat,
                           const R_len_t new_height = 0,
                           const R_len_t new_width = 0,
                           const bool front = false) {
  R_len_t img_c, img_r;
  img_c = mat.ncol();
  img_r = mat.nrow();
  R_len_t new_h = new_height < 0 ? std::max(0, img_r + new_height) : new_height;
  R_len_t new_w = new_width < 0 ? std::max(0, img_c +new_width) : new_width;  
  
  // new dimensions are larger than original dimensions, no need to crop
  if((img_c <= new_w) && (img_r <= new_h)) return mat;
  // new width is larger than original width and new_h is 0, no need to crop
  if((img_c <= new_w) && (new_h == 0)) return mat;
  // new height is larger than original height and new_w is 0, no need to crop
  if((img_r <= new_h) && (new_w == 0)) return mat;
  
  R_len_t ori_c, ori_r, crop_c, crop_r;
  // height resizing
  if((new_h > 0) && (new_h < img_r)) {
    // compute same amount of rows to remove on top-bottom
    ori_r = (img_r - new_h) >> 1;
    if(front) {
      crop_r = img_r - ori_r - 1;
      if((img_r - new_h) % 2) ori_r++;
    } else {
      crop_r = new_h + ori_r - 1;
    }
  } else { // no height resizing
    ori_r = 0;
    crop_r = img_r - 1;
  }
  // width resizing
  if((new_w > 0) && (new_w < img_c)) {
    // compute same amount of cols to remove on right-left
    ori_c = (img_c - new_w) >> 1;
    if(front) {
      crop_c = img_c - ori_c - 1;
      if((img_c - new_w) % 2) ori_c++;
    } else {
      crop_c = new_w + ori_c - 1;
    }
  } else { // no width resizing
    ori_c = 0;
    crop_c = img_c - 1;
  }
  return mat( Rcpp::Range(ori_r, crop_r) , Rcpp::Range(ori_c, crop_c));
}

// function to expand matrix without adding noise
template <int RTYPE> 
SEXP expand_no_noise_T(Rcpp::Matrix<RTYPE> mat,
                       const R_len_t new_height = 0,
                       const R_len_t new_width = 0,
                       const double bg = 0.0,
                       const bool front = false) {
  R_len_t img_r, img_c;
  img_r = mat.nrow();
  img_c = mat.ncol();
  // new dimensions are smaller than original dimensions, no need to expand
  if((img_r >= new_height) && (img_c >= new_width)) return mat;
  
  R_len_t i_row, i_col, ori_r, ori_c, i_out, fin_height, fin_width;
  // final height: either original or desired height if larger
  fin_height = img_r > new_height ? img_r : new_height;
  // final width: either original or desired width if larger
  fin_width = img_c > new_width ? img_c : new_width;
  // compute padding
  ori_r = (fin_height-img_r) >> 1;
  ori_c = (fin_width-img_c) >> 1;
  
  // create output matrix with expanded dimension and fill with bg
  Rcpp::Matrix<RTYPE> out = Rcpp::no_init(fin_height, fin_width);
  out.fill(bg);
  
  // write mat into center of output matrix
  if(front) {
    if((fin_height-img_r) % 2) ori_r++;
    if((fin_width-img_c) % 2) ori_c++;
    for(i_col = img_c - 1 + ori_c, i_out = mat.size() - 1; i_col >= ori_c; i_col--) {
      for(i_row = img_r - 1 + ori_r; i_row >= ori_r; i_row--) {
        out(i_row, i_col) = mat[i_out--];
      }
    }
  } else {
    for(i_col = 0; i_col < img_c; i_col++) {
      i_out = (i_col + ori_c) * fin_height + ori_r;
      for(i_row = 0; i_row < img_r; i_row++) {
        out[i_out++] = mat(i_row, i_col);
      }
    }
  }
  return out;
}

// function to expand matrix with new rows adding noisy padding bg
template <int RTYPE> 
SEXP expand_row_T(Rcpp::Matrix<RTYPE> mat,
                  const R_len_t new_height = 0,
                  const double bg = 0.0,
                  const double sd = 0.0,
                  const bool front = false) {
  R_len_t img_r = mat.nrow();
  // new dimensions are smaller than original dimensions, no need to expand
  if(img_r >= new_height) return mat;
  
  R_len_t img_c, ori_r, i;
  img_c = mat.ncol();
  // compute row padding
  ori_r = (new_height-img_r) >> 1;
  if(front) if((new_height-img_r) % 2) ori_r++;
  
  // create output matrix with expanded rows
  Rcpp::Matrix<RTYPE> out = Rcpp::no_init(new_height, img_c);
  
  // write top padding
  for(i = 0; i < ori_r; i++) out(i, Rcpp::_) = Rcpp::rnorm(img_c, bg, sd);
  // write mat into center of output matrix
  for(; i < (ori_r + img_r); i++) out(i, Rcpp::_) = mat(i - ori_r, Rcpp::_);
  // write bottom padding
  for(; i < new_height; i++) out(i, Rcpp::_) = Rcpp::rnorm(img_c, bg, sd); 
  return out;
}

// function to expand matrix with new columns adding noisy padding bg
template <int RTYPE> 
SEXP expand_col_T(Rcpp::Matrix<RTYPE> mat,
                  const R_len_t new_width = 0,
                  const double bg = 0.0,
                  const double sd = 0.0,
                  const bool front = false) {
  R_len_t img_c = mat.ncol();
  // new dimensions are smaller than original dimensions, no need to expand
  if(img_c >= new_width) return mat;
  
  R_len_t img_r, ori_c, i;
  img_r = mat.nrow();
  // compute row padding
  ori_c = (new_width-img_c) >> 1;
  if(front) if((new_width-img_c) % 2) ori_c++;
  
  // create output matrix with expanded rows
  Rcpp::Matrix<RTYPE> out = Rcpp::no_init(img_r, new_width);
  
  // write left padding
  for(i = 0; i < ori_c; i++) out(Rcpp::_, i) = Rcpp::rnorm(img_r, bg, sd);
  // write mat into center of output matrix
  for(; i < (ori_c + img_c); i++) out(Rcpp::_, i) = mat(Rcpp::_, i - ori_c);
  // write right padding
  for(; i < new_width; i++) out(Rcpp::_, i) = Rcpp::rnorm(img_r, bg, sd);
  
  return out;
}

// function to expand matrix with padding noisy bg
template <int RTYPE> 
SEXP expand_w_noise_T(Rcpp::Matrix<RTYPE> mat,
                      const R_len_t new_height = 0,
                      const R_len_t new_width = 0,
                      const double bg = 0.0,
                      const double sd = 0.0,
                      const bool front = false) {
  Rcpp::Matrix<RTYPE> M0 = expand_col_T(mat, new_width, bg, sd, front);
  return expand_row_T(M0, new_height, bg, sd, front);
}

template <int RTYPE> 
SEXP resize_T(Rcpp::Matrix<RTYPE> mat,
              const R_len_t new_height = 0,
              const R_len_t new_width = 0,
              const bool add_noise = true,
              const double bg = 0.0,
              const double sd = 0.0,
              const bool front = false) {
  Rcpp::Matrix<RTYPE> crop = crop_T(mat, new_height, new_width, front);
  Rcpp::Matrix<RTYPE> out;
  if(add_noise) {
    out = expand_w_noise_T(crop, new_height, new_width, bg, sd, front);
  } else {
    out = expand_no_noise_T(crop, new_height, new_width, bg, front);
  }
  if(mat.hasAttribute("mask")) out.attr("mask") = mat.attr("mask");
  return out;
}

//' @title Header for Matrix Cropping
//' @name cpp_crop
//' @description
//' Crops mat according to new_height and new_width parameters.
//' @param mat Matrix (Raw, Integer, Logical or Numeric).
//' @param new_height a R_len_t, giving the new height of returned mat. Default is 0 for no change. Negative values will remove rows from mat.
//' @param new_width a R_len_t, giving the new width of returned mat. Default is 0 for no change. Negative values will remove cols from mat.
//' @param front a bool, whether to apply cropping from front. Default is false.
//' @return a cropped matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP hpp_crop (SEXP mat,
               const R_len_t new_height = 0,
               const R_len_t new_width = 0, 
               const bool front = false) {
  switch(TYPEOF(mat)) {
  // case NILSXP : return mat;
  case RAWSXP : return crop_T<RAWSXP>(mat, new_height, new_width, front);
  case LGLSXP : return crop_T<LGLSXP>(mat, new_height, new_width, front);
  case INTSXP : return crop_T<INTSXP>(mat, new_height, new_width, front);;
  case REALSXP : return crop_T<REALSXP>(mat, new_height, new_width, front);
  default: Rcpp::stop("hpp_crop: not supported type in 'mat'");
  }
}

//' @title Header for Matrix Resizing
//' @name cpp_resize
//' @description
//' Resizes mat according to new_height and new_width parameters.
//' @param mat Matrix (Raw, Integer, Logical or Numeric).
//' @param new_height a R_len_t, giving the new height of returned mat. Default is 0 for no change. Negative values will remove rows from mat.
//' @param new_width a R_len_t, giving the new width of returned mat. Default is 0 for no change. Negative values will remove cols from mat.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions. Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.
//' @param front a bool, whether to apply cropping from front. Default is false.
//' @return a resized matrix with padding background if new_height or new_width is larger than original mat dimensions.
//' @keywords internal
////' @export
// [[Rcpp::export]]
SEXP hpp_resize (SEXP mat,
                 const R_len_t new_height = 0, 
                 const R_len_t new_width = 0,
                 const bool add_noise = true, 
                 const double bg = 0.0,
                 const double sd = 0.0,
                 const bool front = false) {
  switch(TYPEOF(mat)) {
  // case NILSXP : return mat;
  case RAWSXP : return resize_T<RAWSXP>(mat, new_height, new_width, add_noise, bg, sd, front);
  case LGLSXP : return resize_T<LGLSXP>(mat, new_height, new_width, add_noise, bg, sd, front);
  case INTSXP : return resize_T<INTSXP>(mat, new_height, new_width, add_noise, bg, sd, front);
  case REALSXP : return resize_T<REALSXP>(mat, new_height, new_width, add_noise, bg, sd, front);
  default: Rcpp::stop("hpp_resize: not supported type in 'mat'");
  }
}

#endif
