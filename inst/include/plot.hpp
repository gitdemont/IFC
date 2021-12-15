/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2020 Yohann Demont                                              
                                                                                
  It is part of IFC package, please cite:                                       
  -IFC: An R Package for Imaging Flow Cytometry                                 
  -YEAR: 2020                                                                   
  -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,              
                      Jean-Pierre Marolleau, Loï?c Gaççon,                       
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

#ifndef IFC_PLOT_HPP
#define IFC_PLOT_HPP

#include <Rcpp.h>
using namespace Rcpp;

// compute distance from point with coords x, y to line
// passing through points x1, y1 and x2, y2
double hpp_dist(const R_len_t x1, const R_len_t y1,
                const R_len_t x2, const R_len_t y2,
                const double x, const double y) {
  R_len_t dx = x2-x1;
  R_len_t dy = y2-y1;
  return std::abs(dx*(y1-y)-dy*(x1-x)) / std::sqrt(dx*dx+dy*dy);
}

// draw line from point x1, y1 to point x2, y2
// in a LogicalMatrix M
Rcpp::LogicalMatrix hpp_line(const R_len_t x1, const R_len_t y1,
                             const R_len_t x2, const R_len_t y2,
                             const Rcpp::LogicalMatrix M) {
  Rcpp::LogicalMatrix out = clone(M);
  R_len_t y = y1, x = x1;
  double offsetx = 0.5, offsety = 0.5;
  if(x2 < x) offsetx = -offsetx;
  if(y2 < y) offsety = -offsety;
  out(y, x) = true;
  while(!((y == y2) && (x == x2))) {
    Rcpp::NumericVector d = Rcpp::NumericVector::create(hpp_dist(x1, y1, x2, y2, x + offsetx, y + offsety),
                                                        hpp_dist(x1, y1, x2, y2, x + offsetx, y),
                                                        hpp_dist(x1, y1, x2, y2, x, y + offsety));
    switch(Rcpp::which_min(d))
    {
    case 0 :
      x = x + 2 * offsetx;
      y = y + 2 * offsety;
      break;
    case 1 :
      x = x + 2 * offsetx;
      break;
    case 2 :
      y = y + 2 * offsety;
      break;
    }
    out(y, x) = true;
  }
  return out;
}

// determines triangle vertices coordinates 
Rcpp::IntegerMatrix hpp_tr_vertices(const uint8_t size = 3) {
  Rcpp::IntegerMatrix out(3, 2);
  int half = size >> 1;
  int l = std::ceil((0.5 + 0.8666) * half);
  out(0, 1) = half;
  out(1, 0) = l;
  out(2, 0) = l;
  out(2, 1) = size > 1 ? size - 1 : 0;
  return out;
}

// create triangle shape logical matrix
Rcpp::LogicalMatrix hpp_triangle(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  Rcpp::IntegerMatrix V = hpp_tr_vertices(size);
  out = hpp_line(V(0,1),V(0,0),V(2,1),V(2,0), out);
  out = hpp_line(V(0,1),V(0,0),V(1,1),V(1,0), out);
  out = hpp_line(V(1,1),V(1,0),V(2,1),V(2,0), out);
  return out;
}

// create filled triangle shape logical matrix
Rcpp::LogicalMatrix hpp_triangle_filled(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out = hpp_triangle(size);
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    Rcpp::LogicalVector V = out(_,i_col);
    R_len_t beg = size, end = 0;
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      if(V[i_row]) {
        beg = i_row;
        break;
      }
    }
    for(R_len_t i_row = size - 1; i_row >= 0; i_row--) {
      if(V[i_row]) {
        end = i_row;
        break;
      }
    }
    for(R_len_t i_row = beg; i_row < end; i_row++) V[i_row] = true;
    out(_,i_col) = V;
  }
  return out;
}

// create filled circle shape logical matrix
Rcpp::LogicalMatrix hpp_circle_filled(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  double half = size % 2 ? size / 2 : size / 2 - 0.5;
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    double foo = i_col - half;
    foo = foo < 0 ? foo + 0.3 : foo - 0.3;
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      double bar = i_row - half;
      bar = bar < 0 ? bar + 0.3 : bar - 0.3;
      out(i_row, i_col) = std::sqrt(foo * foo + bar * bar) <= half;
    }
  }
  return out;
}

// create circle shape logical matrix
Rcpp::LogicalMatrix hpp_circle(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  double half = size % 2 ? size / 2 : size / 2 - 0.5;
  double half_1 = half - 1;
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    double foo = i_col - half;
    foo = foo < 0 ? foo + 0.3 : foo - 0.3;
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      double bar = i_row - half;
      bar = bar < 0 ? bar + 0.3 : bar - 0.3;
      double val = std::sqrt(foo * foo + bar * bar);
      out(i_row, i_col) = (val <= half) && (val > half_1);
    }
  }
  return out;
}

// create square shape logical matrix
Rcpp::LogicalMatrix hpp_square(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    out(0, i_col) = true;
    out(size - 1, i_col) = true; 
  }
  for(R_len_t i_row = 1; i_row < size - 1; i_row++) {
    out(i_row, 0) = true;
    out(i_row, size - 1) = true; 
  }
  return out;
}

// create filled square shape logical matrix
Rcpp::LogicalMatrix hpp_square_filled(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  out.fill(true);
  return out;
}

// create filled diamond shape logical matrix
Rcpp::LogicalMatrix hpp_diamond_filled(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  double half = size >> 1;
  R_len_t i = 0;
  if(size % 2) {
    for(R_len_t i_col = -half; i_col <= half; i_col++) {
      for(R_len_t i_row = -half; i_row <= half; i_row++) {
        if(i >= size * size) Rcpp::stop("Not allowed");
        out[i++] = (abs(i_col) + abs(i_row)) <= half;
      }
    }
    
  } else {
    for(R_len_t i_col = -half; i_col <= half; i_col++) {
      if(i_col == 0) i_col++;
      for(R_len_t i_row = -half; i_row <= half; i_row++) {
        if(i_row == 0) i_row++;
        if(i >= size * size) Rcpp::stop("Not allowed");
        out[i++] = (abs(i_col) + abs(i_row)) <= half;
      }
    }
  }
  return out;
}

// create diamond shape logical matrix
Rcpp::LogicalMatrix hpp_diamond(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  double half = size >> 1;
  double half_1 = half - 1;
  R_len_t i = 0;
  if(size % 2) {
    for(R_len_t i_col = -half; i_col <= half; i_col++) {
      for(R_len_t i_row = -half; i_row <= half; i_row++) {
        if(i >= size * size) Rcpp::stop("Not allowed");
        double bar = (abs(i_col) + abs(i_row));
        out[i++] = (bar <= half) && (bar > half_1);
      }
    }
    
  } else {
    for(R_len_t i_col = -half; i_col <= half; i_col++) {
      if(i_col == 0) i_col++;
      for(R_len_t i_row = -half; i_row <= half; i_row++) {
        if(i_row == 0) i_row++;
        if(i >= size * size) Rcpp::stop("Not allowed");
        double bar = (abs(i_col) + abs(i_row));
        out[i++] = (bar <= half) && (bar > half_1);
      }
    }
  }
  return out;
}

// create plus shape logical matrix
Rcpp::LogicalMatrix hpp_plus(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  double half = size % 2 ? size / 2 : size / 2 - 0.5;
  out(half, Rcpp::_) = Rcpp::rep(true, size);
  out(Rcpp::_, half) = Rcpp::rep(true, size);
  return out;
}

// create cross shape logical matrix
Rcpp::LogicalMatrix hpp_cross(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      out(i_row, i_col) = i_row == i_col || i_row == (size - 1 - i_col);
    }
  }
  return out;
}

// create gausian kernel for blurring
Rcpp::NumericMatrix hpp_gaussian(const uint8_t size = 3,
                                 const uint8_t sigma = 3) {
  Rcpp::NumericMatrix out(size, size);
  if(size == 0) return out;
  double half = size % 2 ? size / 2 : size / 2 - 0.5;
  double d = sigma * std::sqrt(2 * M_PI);
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    double foo = i_col - half;
    for(R_len_t i_row = 0; i_row < size; i_row++) {
      out(i_row, i_col) = std::exp(-std::sqrt(foo * foo + pow(i_row - half, 2.0))/(2*std::pow(sigma, 2.0))) / d;
    }
  }
  return out;
}

// check and get image array dimensions
// img is expected to be a non null IntegerVector with dimension attribute [nrow, ncol, 4]
Rcpp::IntegerVector get_dim(Rcpp::IntegerVector img) {
  Rcpp::Nullable<Rcpp::IntegerVector> x_ = img.attr("dim");
  if(x_.isNotNull()) {
    Rcpp::IntegerVector x(x_.get());
    if(x.size() == 3) {
      if(x[2] != 4) Rcpp::stop("'img' should be a 3D array of with rgba colors in 3rd dimension");
      return x;
    } else {
      Rcpp::stop("'img' should be a 3D array");
    }
  } else {
    Rcpp::stop("'img' should be a 3D array");
  }
}

//' @title Draw Shape to Image
//' @name cpp_drawat
//' @description low-level function to add shape on image
//' @param img an IntegerVector. A non null array of dimensions [nrow, ncol, 4].
//' @param coords an IntegerMatrix whose rows are points to draw and with:\cr
//' - 1st column being img col coordinate in px,\cr
//' - 2nd column being img row coordinate in px,\cr
//' - 3rd column determining if point should be drawn or not.
//' @param mask a LogicalMatrix where every true value will be added to the image.
//' @param color, a length 4 IntegerMatrix specifying rgba, from 0 to 255. Only first row will be used.
//' @details shape according to 'mask' will be drawn on 'img' centered at coordinates img[coord[, 1], coord[, 0]] if coord[, 2] is true.\cr
//' Every pixels being part of the shape will be filled with 'color'.
//' If only one 'color' is provided, this 'color' will be used for each points.
//' If more than one 'color' is provided, then if number of colors (ncol) equals the number of points 'color' will be used as is for each single point.
//' Otherwise, 'color' will be considered as a color-gradient and density will be computed.
//' /!\ please note that IFC:::densCols() is faster to compute color based on density for n < 20000 points, so it's worth using it when number of points are lower.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerVector hpp_drawat(const Rcpp::IntegerVector img,
                               const Rcpp::IntegerMatrix coords = Rcpp::IntegerMatrix(1,2),
                               const Rcpp::LogicalMatrix mask = Rcpp::LogicalMatrix(1),
                               const Rcpp::IntegerMatrix color = Rcpp::IntegerMatrix(4,1)) {
  R_len_t msk_c = mask.ncol() >> 1;
  R_len_t msk_r = mask.nrow() >> 1; 
  R_len_t col_r = color.nrow();
  R_len_t col_c = color.ncol();
  if(col_r != 4) Rcpp::stop("hpp_dens: bad density color specification");
  Rcpp::IntegerVector V = get_dim(img);
  R_len_t width  = V[1];
  R_len_t height = V[0];
  Rcpp::IntegerVector out = clone(img);
  if(color.size() == 4) { // only one-color points
    for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
      R_len_t i_row = coords(i_pt, 1);
      R_len_t i_col = coords(i_pt, 0);
      if(!(coords(i_pt, 2))) continue;
      for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col <= i_col + msk_c; f_col++) {
        for(R_len_t f_row = i_row - msk_r; f_row <= i_row + msk_r; f_row++) {
          if(mask[i_msk++] &&
             (f_col >= 0) &&
             (f_col < width) &&
             (f_row >= 0) &&
             (f_row < height)) {
            for(R_len_t i_k = 0; i_k < 4; i_k++) {
              out[i_k * height * width + f_col * height + f_row] = color[i_k];
            }
          }
        }
      }
    }
  } else {
    if(col_c == coords.nrow()) { // colors are provided for each points
      for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
        R_len_t i_row = coords(i_pt, 1);
        R_len_t i_col = coords(i_pt, 0);
        if(!(coords(i_pt, 2))) continue;
        for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col <= i_col + msk_c; f_col++) {
          for(R_len_t f_row = i_row - msk_r; f_row <= i_row + msk_r; f_row++) {
            if(mask[i_msk++] &&
               (f_col >= 0) &&
               (f_col < width) &&
               (f_row >= 0) &&
               (f_row < height)) {
              for(R_len_t i_k = 0; i_k < 4; i_k++) {
                out[i_k * height * width + f_col * height + f_row] = color(i_k,i_pt);
              }
            }
          }
        }
      }
    } else { // colors are provided as a gradient, we compute density
      Rcpp::IntegerMatrix grd(height, width);
      Rcpp::NumericMatrix den(height, width);
      Rcpp::NumericMatrix blur = hpp_gaussian(9,3);
      R_len_t blr_c = blur.ncol() >> 1;
      R_len_t blr_r = blur.nrow() >> 1;
      double den_mx = 0;
      for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
        R_len_t i_row = coords(i_pt, 1);
        R_len_t i_col = coords(i_pt, 0);
        for(R_len_t f_col = i_col - blr_c, i_blr = 0; f_col <= i_col + blr_c; f_col++) {
          for(R_len_t f_row = i_row - blr_r; f_row <= i_row + blr_r; f_row++) {
            if((f_col >= 0) &&
               (f_col < width) &&
               (f_row >= 0) &&
               (f_row < height)) {
              den(f_row, f_col) = den(f_row, f_col) + blur[i_blr++];
              if(den(f_row, f_col) > den_mx) den_mx = den(f_row, f_col);
            }
          }
        }
      }
      for(R_len_t i = 0; i < grd.size(); i++) grd[i] = (col_c - 0.001) * std::asinh(den[i] / den_mx) / std::asinh(1);
      for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
        R_len_t i_row = coords(i_pt, 1);
        R_len_t i_col = coords(i_pt, 0);
        if((!coords(i_pt, 2)) ||
           (i_col < 0) ||
           (i_col >= width) ||
           (i_row < 0) ||
           (i_row >= height)) continue;
        for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col <= i_col + msk_c; f_col++) {
          for(R_len_t f_row = i_row - msk_r; f_row <= i_row + msk_r; f_row++) {
            if(mask[i_msk++] &&
               (f_col >= 0) &&
               (f_col < width) &&
               (f_row >= 0) &&
               (f_row < height)) {
              for(R_len_t i_k = 0; i_k < 4; i_k++) {
                out[i_k * height * width + f_col * height + f_row] = color(i_k,grd(i_row, i_col));
              }
            }
          }
        }
      }
    }
  }
  return out;
}

//' @title Raster Image
//' @name cpp_raster
//' @description low-level function to create plot raster
//' @param width a uint16_t determining the returned image width.
//' @param height a uint16_t determining the returned image height.
//' @param obj a List containing drawing information:\cr
//' - pch, an integer specifying a symbol to draw. Handled are 0,1,2,3,4,5,15,17,18,20. Otherwise only a pixel will be drawn.\cr
//' - size, an integer specifying the size in pixel of the shape, from 1 to 255.\cr
//' - color an IntegerMatrix (rgba) of the color used to draw the shape.\cr
//' - coords, an IntegerMatrix whose rows are points to draw and with:\cr
//' -* 1st column being img col coordinate in px,\cr
//' -* 2nd column being img row coordinate in px,\cr
//' -* 3rd column determining if point should be drawn or not.
//' @details shape according to 'pch' will be drawn on 'img' centered at coordinates img[coord[, 1], coord[, 0]] if coord[, 2] is true.\cr
//' Every pixels being part of the shape will be filled with 'color'.
//' If only one 'color' is provided, this 'color' will be used for each points.
//' If more than one 'color' is provided, then if number of colors (ncol) equals the number of points 'color' will be used as is for each single point.
//' Otherwise, 'color' will be considered as a color-gradient and density will be computed.
//' /!\ please note that IFC:::densCols() is faster to compute color based on density for n < 20000 points, so it's worth using it when number of points are lower.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerVector hpp_raster(const uint16_t width,
                               const uint16_t height,
                               const Rcpp::List obj) {
  Rcpp::IntegerVector img(width * height * 4);
  img.attr("dim") = Rcpp::Dimension(height, width, 4);
  for(R_len_t i_obj = 0; i_obj < obj.size(); i_obj++) {
    Rcpp::List L = obj[i_obj];
    R_len_t pch = L["pch"];
    R_len_t size = L["size"];
    Rcpp::IntegerMatrix color = L["col"];
    R_len_t col_r = color.nrow();
    R_len_t col_c = color.ncol();
    if(col_r != 4) Rcpp::stop("hpp_dens: bad density color specification");
    for(R_len_t i = 0; i < color.size(); i++) if((color[i] < 0) || (color[i] > 255)) Rcpp::stop("cpp_dens: bad color specification, out-of-range [0-255]");
    Rcpp::IntegerMatrix coords = L["coords"];
    Rcpp::LogicalMatrix mask;
    switch(pch) {
    case 0 :
      mask = hpp_square(size); // size reference
      break;
    case 1 :
      mask = hpp_circle(size);
      break;
    case 2 :
      mask = hpp_triangle(size + 2);
      break;
    case 3 :
      mask = hpp_plus(size + 2);
      break;
    case 4 :
      mask = hpp_cross(size);
      break;
    case 5 :
      mask = hpp_diamond(size + 2);
      break;
    case 15 :
      mask = hpp_square_filled(size);
      break;
    case 17 :
      mask = hpp_triangle_filled(size + 2);
      break;
    case 18 :
      mask = hpp_diamond_filled(size);
      break;
    case 20 :
      mask = hpp_circle_filled(size - 2);
      break;
    default :
      mask = hpp_square_filled(1);
    break;
    }
    R_len_t msk_c = mask.ncol() >> 1;
    R_len_t msk_r = mask.nrow() >> 1;
    
    if(color.size() == 4) { // only one-color points
      for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
        R_len_t i_row = coords(i_pt, 1);
        R_len_t i_col = coords(i_pt, 0);
        if(!(coords(i_pt, 2))) continue;
        for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col <= i_col + msk_c; f_col++) {
          for(R_len_t f_row = i_row - msk_r; f_row <= i_row + msk_r; f_row++) {
            if(mask[i_msk++] &&
               (f_col >= 0) &&
               (f_col < width) &&
               (f_row >= 0) &&
               (f_row < height)) {
              for(R_len_t i_k = 0; i_k < 4; i_k++) {
                img[i_k * height * width + f_col * height + f_row] = color[i_k];
              }
            }
          }
        }
      }
    } else {
      if(col_c == coords.nrow()) { // colors are provided for each points
        for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
          R_len_t i_row = coords(i_pt, 1);
          R_len_t i_col = coords(i_pt, 0);
          if(!(coords(i_pt, 2))) continue;
          for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col <= i_col + msk_c; f_col++) {
            for(R_len_t f_row = i_row - msk_r; f_row <= i_row + msk_r; f_row++) {
              if(mask[i_msk++] &&
                 (f_col >= 0) &&
                 (f_col < width) &&
                 (f_row >= 0) &&
                 (f_row < height)) {
                for(R_len_t i_k = 0; i_k < 4; i_k++) {
                  img[i_k * height * width + f_col * height + f_row] = color(i_k,i_pt);
                }
              }
            }
          }
        }
      } else { // colors are provided as a gradient, we compute density
        Rcpp::IntegerMatrix grd(height, width);
        Rcpp::NumericMatrix den(height, width);
        Rcpp::NumericMatrix blur = hpp_gaussian(9,3);
        R_len_t blr_c = blur.ncol() >> 1;
        R_len_t blr_r = blur.nrow() >> 1;
        double den_mx = 0;
        for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
          R_len_t i_row = coords(i_pt, 1);
          R_len_t i_col = coords(i_pt, 0);
          for(R_len_t f_col = i_col - blr_c, i_blr = 0; f_col <= i_col + blr_c; f_col++) {
            for(R_len_t f_row = i_row - blr_r; f_row <= i_row + blr_r; f_row++) {
              if((f_col >= 0) &&
                 (f_col < width) &&
                 (f_row >= 0) &&
                 (f_row < height)) {
                den(f_row, f_col) = den(f_row, f_col) + blur[i_blr++];
                if(den(f_row, f_col) > den_mx) den_mx = den(f_row, f_col);
              }
            }
          }
        }
        for(R_len_t i = 0; i < grd.size(); i++) grd[i] = (col_c - 0.001) * std::asinh(den[i] / den_mx) / std::asinh(1);
        for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
          R_len_t i_row = coords(i_pt, 1);
          R_len_t i_col = coords(i_pt, 0);
          if((!coords(i_pt, 2)) ||
             (i_col < 0) ||
             (i_col >= width) ||
             (i_row < 0) ||
             (i_row >= height)) continue;
          for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col <= i_col + msk_c; f_col++) {
            for(R_len_t f_row = i_row - msk_r; f_row <= i_row + msk_r; f_row++) {
              if(mask[i_msk++] &&
                 (f_col >= 0) &&
                 (f_col < width) &&
                 (f_row >= 0) &&
                 (f_row < height)) {
                for(R_len_t i_k = 0; i_k < 4; i_k++) {
                  img[i_k * height * width + f_col * height + f_row] = color(i_k,grd(i_row, i_col));
                }
              }
            }
          }
        }
      }
    }
  }
  return img;
}

#endif
