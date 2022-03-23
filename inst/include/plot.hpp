/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2022 Yohann Demont                                              
                                                                                
  It is part of IFC package, please cite:                                       
  -IFC: An R Package for Imaging Flow Cytometry                                 
  -YEAR: 2020                                                                   
  -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,              
                      Jean-Pierre Marolleau, Loïc Garçon,                       
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
////' @export
// [[Rcpp::export]]
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
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix hpp_tr_vertices(const uint8_t size = 3) {
  Rcpp::IntegerMatrix out(3, 2);
  if(size == 0) return out;
  int half = size >> 1;
  int l = std::round((0.5 + 0.8666) * half); // TODO triangles could be less larger
  l = std::min(l, size - 1);
  // out(0, 0) = 0;
  out(0, 1) = half;
  out(1, 0) = l;
  // out(1, 1) = 0;
  out(2, 0) = l;
  out(2, 1) = size > 1 ? size - 1 : 0;
  return out;
}

// create triangle shape logical matrix
////' @export
// [[Rcpp::export]]
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
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_triangle_filled(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out = hpp_triangle(size);
  if(size == 0) return out;
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

// create triangle inside square shape
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_fourteen(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  R_len_t size_1 = size - 1;
  for(R_len_t i_col = 0; i_col < size; i_col++) {
    out(0, i_col) = true;
    out(size_1, i_col) = true; 
  }
  for(R_len_t i_row = 1; i_row < size_1; i_row++) {
    out(i_row, 0) = true;
    out(i_row, size_1) = true; 
  }
  R_len_t half = size >> 1;
  out = hpp_line(half, 0, 0, size_1, out);
  out = hpp_line(half, 0, size_1, size_1, out);
  return out;
}

// create filled circle shape logical matrix
////' @export
// [[Rcpp::export]]
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
////' @export
// [[Rcpp::export]]
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
////' @export
// [[Rcpp::export]]
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
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_square_filled(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  out.fill(true);
  return out;
}

// create filled diamond shape logical matrix
////' @export
// [[Rcpp::export]]
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
////' @export
// [[Rcpp::export]]
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
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_plus(const uint8_t size = 3) {
  Rcpp::LogicalMatrix out(size, size);
  if(size == 0) return out;
  double half = size % 2 ? size / 2 : size / 2 - 0.5;
  out(half, Rcpp::_) = Rcpp::rep(true, size);
  out(Rcpp::_, half) = Rcpp::rep(true, size);
  return out;
}

// create cross shape logical matrix
////' @export
// [[Rcpp::export]]
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

// reverse shape
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_shape_rev(const Rcpp::LogicalMatrix M) {
  Rcpp::LogicalMatrix out = Rcpp::no_init(M.nrow(), M.ncol());
  for(R_len_t i = 0; i < M.size(); i++) out[i] = M[M.size() - 1 - i];
  return out;
}

// combine 2 shapes to create a new one
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hpp_shape_combine(const Rcpp::LogicalMatrix M1,
                                      const Rcpp::LogicalMatrix M2) {
  Rcpp::LogicalMatrix out;
  Rcpp::LogicalMatrix foo;
  R_len_t pad = std::abs(M1.ncol() - M2.ncol()) >> 1;
  if(M1.ncol() > M2.ncol()) {
    out = Rcpp::clone(M1);
    foo = Rcpp::clone(M2);
  } else {
    out = Rcpp::clone(M2);
    foo = Rcpp::clone(M1);
  }
  for(R_len_t i = 0, i_col = pad; i_col < out.ncol() - pad; i_col++) {
    for(R_len_t i_row = pad; i_row < out.nrow() - pad; i_row++, i++) {
      out(i_row, i_col) = out(i_row, i_col) || foo[i];
    }
  }
  return out;
}

// create gausian kernel for blurring
Rcpp::NumericMatrix hpp_gaussian(const uint8_t size = 3,
                                 const double sigma = 3.0) {
  if((size == 0) || (size >= 35)) Rcpp::stop("hpp_gaussian: 'size' argument is not possible for blurring");
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


//' @title Image to Native Raster Conversion
//' @name cpp_as_nativeRaster
//' @description Converts 3D image array to nativeRaster
//' @param x an IntegerVecter /!\ It should be coercible to 3D array [height, width, rgba]
//' @return a nativeRaster IntegerMatrix
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix hpp_as_nativeRaster(const Rcpp::IntegerVector x) {
  Rcpp::IntegerVector V = get_dim(x);
  R_len_t h = V[0], w = V[1];
  R_len_t D = h * w;
  Rcpp::IntegerMatrix out = Rcpp::no_init(h, w);
  for(R_len_t i_row = 0, i = 0; i_row < h; i_row++) {
    for(R_len_t i_col = 0; i_col < w; i_col++, i++) {
      R_len_t d = i_col * h + i_row;
      out[i] = 
        x[d] | // x[0 * D + d]
        (x[1 * D + d] <<  8) |
        (x[2 * D + d] << 16) |
        (x[3 * D + d] << 24);
    }
  }
  out.attr("class") = "nativeRaster";
  return out;
}

//' @title Coordinates to Pixels
//' @name cpp_coord_to_px
//' @description low-level function to compute pixels coordinates
//' @param x NumericVector of x-coordinates of the points.
//' @param y NumericVector of y-coordinates of the points.
//' @param param NumericVector of parameters to scale raw points coordinates to pixels coordinates.
//' @return a 2 columns IntegerMatrix of x and y pixels coordinates.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix hpp_coord_to_px(const Rcpp::NumericVector x,
                                    const Rcpp::NumericVector y,
                                    const Rcpp::NumericVector param) {
  if(x.size() != y.size()) Rcpp::stop("cpp_coord_to_px: 'x' and 'y' should be of same size");
  Rcpp::IntegerMatrix out(x.size(), 2);
  if(param[12]) {
    for(R_len_t i = 0; i < x.size(); i++) {
      double v;
      if(x[i] < param[0]) {
        v = param[0];
      } else {
        if(x[i] > param[1]) {
          v = param[1];
        } else {
          v = x[i];
        }
      }
      out(i, 0) = ((v - param[6]) * param[4] + param[8]) / param[10];
      if(y[i] < param[2]) {
        v = param[2];
      } else {
        if(x[i] > param[3]) {
          v = param[3];
        } else {
          v = y[i];
        }
      }
      out(i, 1) = (param[9] - (v - param[7]) * param[5]) / param[11];
    }
    return out;
  } else {
    R_len_t n = 0;
    for(R_len_t i = 0; i < x.size(); i++) {
      if((x[i] >= param[0]) && (x[i] <= param[1]) && (y[i] >= param[2]) && (y[i] <= param[3])) {
        out(n, 0) = ((x[i] - param[6]) * param[4] + param[8]) / param[10];
        out(n, 1) = (param[9] - (y[i] - param[7]) * param[5]) / param[11];
        n++;
      }
    }
    return out(Rcpp::Range(0, n - 1), Rcpp::_);
  }
}

//' @title Draw Shape to Image
//' @name cpp_draw
//' @description low-level function to add shape on image
//' @param img an IntegerVector. A non null array of dimensions [nrow, ncol, 4].
//' @param coords an IntegerMatrix whose rows are points to draw and with:\cr
//' - 1st column being img col coordinate in px,\cr
//' - 2nd column being img row coordinate in px.
//' @param mask a LogicalMatrix where every true value will be added to the image.
//' @param color, a 4 rows IntegerMatrix specifying rgba, from 0 to 255.
//' @param blur_size, a uint8_t the size of the gaussian blurring kernel. Default is 9.
//' @param blur_sd, a double the sd of the gaussian blurring kernel. Default is 3.0.
//' @details shape according to 'mask' will be drawn on 'img' centered at coordinates img[coord[, 1], coord[, 0]] if coord[, 2] is true.\cr
//' Every pixels being part of the shape will be filled with 'color'.
//' If only one 'color' is provided, this 'color' will be used for each points.
//' If more than one 'color' is provided, then if number of colors (ncol) equals the number of points 'color' will be used as is for each single point.
//' Otherwise, 'color' will be considered as a color-gradient and density will be computed.
//' /!\ please note that IFC:::densCols() is faster to compute color based on density for n < 20000 points, so it's worth using it when number of points are lower.
//' @keywords internal
//' @return /!\ nothing is returned but img is modified in-place
////' @export
// [[Rcpp::export]]
void hpp_draw(Rcpp::IntegerVector img,
              const Rcpp::IntegerMatrix coords = Rcpp::IntegerMatrix(1,2),
              const Rcpp::LogicalMatrix mask = Rcpp::LogicalMatrix(1),
              const Rcpp::IntegerMatrix color = Rcpp::IntegerMatrix(4,1),
              const uint8_t blur_size = 9,
              const double blur_sd = 3.0) {
  if((mask.size() == 0) || (mask.size() >= 1225)) Rcpp::stop("hpp_draw: 'size' argument is not possible with this shape");
  R_len_t msk_c = mask.ncol() >> 1;
  R_len_t msk_r = mask.nrow() >> 1;
  R_len_t msk_c_1 = msk_c + (mask.ncol() % 2);
  R_len_t msk_r_1 = msk_r + (mask.nrow() % 2);
  R_len_t col_r = color.nrow();
  R_len_t col_c = color.ncol();
  if(col_r != 4) Rcpp::stop("hpp_draw: bad 'color' specification");
  for(R_len_t i = 0; i < color.size(); i++) if((color[i] < 0) || (color[i] > 255)) Rcpp::stop("hpp_draw: bad 'color' specification, out-of-range [0-255]");
  Rcpp::IntegerVector V = get_dim(img);
  R_len_t width  = V[1];
  R_len_t height = V[0];
  unsigned short count = 1;
  if(color.size() == 4) { // only one-color points
    Rcpp::LogicalMatrix Z(V[0], V[1]); // matrix to record points already drawn so a to skip drawing another point at same xy location
    Z.fill(true);
    for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
      if((count++ % 10000) == 0) {
        count = 1;
        Rcpp::checkUserInterrupt();
      }
      R_len_t i_row = coords(i_pt, 1);
      R_len_t i_col = coords(i_pt, 0);
      if(Z(i_row, i_col)) {
        Z(i_row, i_col) = false; // no need to draw same point at same xy location
        for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col < i_col + msk_c_1; f_col++) {
          for(R_len_t f_row = i_row - msk_r; f_row < i_row + msk_r_1; f_row++, i_msk++) {
            if(mask[i_msk] &&
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
    }
  } else {
    if(col_c == coords.nrow()) { // colors are provided for each points
      for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
        if((count++ % 10000) == 0) {
          count = 1;
          Rcpp::checkUserInterrupt();
        }
        R_len_t i_row = coords(i_pt, 1);
        R_len_t i_col = coords(i_pt, 0);
        for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col < i_col + msk_c_1; f_col++) {
          for(R_len_t f_row = i_row - msk_r; f_row < i_row + msk_r_1; f_row++, i_msk++) {
            if(mask[i_msk] &&
               (f_col >= 0) &&
               (f_col < width) &&
               (f_row >= 0) &&
               (f_row < height)) {
              for(R_len_t i_k = 0; i_k < 4; i_k++) {
                img[i_k * height * width + f_col * height + f_row] = color[i_k + 4 * i_pt];
              }
            }
          }
        }
      }
    } else { // colors are provided as a gradient, we compute density
      Rcpp::IntegerMatrix grd(height, width);
      Rcpp::NumericMatrix den(height, width);
      Rcpp::NumericMatrix blur = hpp_gaussian(blur_size, blur_sd);
      R_len_t blr_c = blur.ncol() >> 1;
      R_len_t blr_r = blur.nrow() >> 1;
      R_len_t blr_c_1 = blr_c + (blur.ncol() % 2);
      R_len_t blr_r_1 = blr_r + (blur.nrow() % 2);
      double den_mx = 0;
      double ash = std::asinh(1);
      for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
        if((count++ % 10000) == 0) {
          count = 1;
          Rcpp::checkUserInterrupt();
        }
        R_len_t i_row = coords(i_pt, 1);
        R_len_t i_col = coords(i_pt, 0);
        for(R_len_t f_col = i_col - blr_c, i_blr = 0; f_col <= i_col + blr_c_1; f_col++) {
          for(R_len_t f_row = i_row - blr_r; f_row <= i_row + blr_r_1; f_row++, i_blr++) {
            if((f_col >= 0) &&
               (f_col < width) &&
               (f_row >= 0) &&
               (f_row < height)) {
              den(f_row, f_col) = den(f_row, f_col) + blur[i_blr];
              if(den(f_row, f_col) > den_mx) den_mx = den(f_row, f_col);
            }
          }
        }
      }
      for(R_len_t i = 0; i < grd.size(); i++) grd[i] = (col_c - 0.001) * std::asinh(den[i] / den_mx) / ash;
      Rcpp::LogicalMatrix Z(V[0], V[1]); // matrix to record points already drawn so a to skip drawing another point at same xy location
      Z.fill(true);
      if(mask.size() == 1) {
        for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
          if((count++ % 10000) == 0) {
            count = 1;
            Rcpp::checkUserInterrupt();
          }
          R_len_t i_row = coords(i_pt, 1);
          R_len_t i_col = coords(i_pt, 0);
          if(Z(i_row, i_col)) {
            Z(i_row, i_col) = false;
            for(R_len_t i_k = 0; i_k < 4; i_k++) {
              img[i_k * height * width + i_col * height + i_row] = color(i_k,grd(i_row, i_col));
            }
          }
        }
      } else {
        for(R_len_t i_pt = 0; i_pt < coords.nrow(); i_pt++) {
          if((count++ % 10000) == 0) {
            count = 1;
            Rcpp::checkUserInterrupt();
          }
          R_len_t i_row = coords(i_pt, 1);
          R_len_t i_col = coords(i_pt, 0);
          if(Z(i_row, i_col)) {
            Z(i_row, i_col) = false;
            R_len_t v = 0;
            for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col < i_col + msk_c_1; f_col++) {
              for(R_len_t f_row = i_row - msk_r; f_row < i_row + msk_r_1; f_row++, i_msk++) {
                if(mask[i_msk] &&
                   (f_col >= 0) &&
                   (f_col < width) &&
                   (f_row >= 0) &&
                   (f_row < height)) {
                  if(grd(f_row, f_col) > v) v = grd(f_row, f_col);
                }
              }
            }
            for(R_len_t f_col = i_col - msk_c, i_msk = 0; f_col < i_col + msk_c_1; f_col++) {
              for(R_len_t f_row = i_row - msk_r; f_row < i_row + msk_r_1; f_row++, i_msk++) {
                if(mask[i_msk] &&
                   (f_col >= 0) &&
                   (f_col < width) &&
                   (f_row >= 0) &&
                   (f_row < height)) {
                  for(R_len_t i_k = 0; i_k < 4; i_k++) {
                    img[i_k * height * width + f_col * height + f_row] = color(i_k, v);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

//' @title Raster Image
//' @name cpp_raster
//' @description low-level function to create plot raster
//' @param width a uint16_t determining the returned image width.
//' @param height a uint16_t determining the returned image height.
//' @param obj a List containing drawing information:\cr
//' - pch, an integer specifying a symbol to draw. Handled are [0-20]. Otherwise only a pixel will be drawn.\cr
//' - size, an integer specifying the size in pixel of the shape, from 1 to 255.\cr
//' - color a 4 rows IntegerMatrix (rgba) of the color used to draw the shape.\cr
//' - coords, an IntegerMatrix whose rows are points to draw and with:\cr
//' -* 1st column being img col coordinate in px,\cr
//' -* 2nd column being img row coordinate in px.
//' - blur_size an integer controlling the size of the blurring gaussian kernel.\cr
//' - blur_sd a double controlling the sd of the blurring gaussian kernel.
//' @details shape according to 'pch' will be drawn on 'img' centered at coordinates img[coord[, 1], coord[, 0]] if coord[, 2] is true
//' and every pixels being part of the shape will be filled with 'color'.
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
    size += 3;
    Rcpp::IntegerMatrix color = L["col"];
    Rcpp::IntegerMatrix coords = L["coords"];
    Rcpp::LogicalMatrix mask;
    switch(pch) {
    case 0 :
      mask = hpp_square(size + 1);
      break;
    case 1 :
      mask = hpp_circle(size);
      break;
    case 2 :
      mask = hpp_triangle(size + 4);
      break;
    case 3 :
      mask = hpp_plus(size + 3);
      break;
    case 4 :
      mask = hpp_cross(size);
      break;
    case 5 :
      mask = hpp_diamond(size + 4);
      break;
    case 6 :
      mask = hpp_shape_rev(hpp_triangle(size + 4));
      break;
    case 7 :
      mask = hpp_shape_combine(hpp_square(size + 1), hpp_cross(size + 1));
      break;
    case 8 :
      mask = hpp_shape_combine(hpp_plus(size + 4), hpp_cross(size + 2));
      break;
    case 9 :
      mask = hpp_shape_combine(hpp_plus(size + 4), hpp_diamond(size + 4));
      break;
    case 10 :
      mask = hpp_shape_combine(hpp_circle(size), hpp_plus(size));
      break;
    case 11 :
      mask = hpp_shape_combine(hpp_triangle(size + 5), hpp_shape_rev(hpp_triangle(size + 5)));
      break;
    case 12 :
      mask = hpp_shape_combine(hpp_square(size + 1), hpp_plus(size + 1));
      break;
    case 13 :
      mask = hpp_shape_combine(hpp_circle(size), hpp_cross(size + 2));
      break;
    case 14 :
      mask = hpp_fourteen(size + 1);
      break;
    case 15 :
      mask = hpp_square_filled(size);
      break;
    case 16 :
      mask = hpp_circle_filled(size);
      break;
    case 17 :
      mask = hpp_triangle_filled(size + 2);
      break;
    case 18 :
      mask = hpp_diamond_filled(size);
      break;
    case 19 :
      mask = hpp_circle_filled(size);
      break;
    case 20 :
      mask = hpp_circle_filled(size - 2); // size reference
      break;
    /* shapes 21 to 25 are not available because they would require 2 colors (border and fill)
    case 21 :
      mask = hpp_circle(size);
      break;
    case 22 :
      mask = hpp_square(size);
      break;
    case 23 :
      mask = hpp_diamond(size + 2);
      break;
    case 24 :
      mask = hpp_triangle(size + 4);
      break;
    case 25 :
      mask = hpp_shape_rev(hpp_triangle(size + 4));
      break;*/
    default :
      mask = hpp_square_filled(1);
    break;
    }
    hpp_draw(img, coords, mask, color, L["blur_size"], L["blur_sd"]);
  }
  return img;
}

#endif
