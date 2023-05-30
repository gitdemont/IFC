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

#ifndef IFC_UTILS_HPP
#define IFC_UTILS_HPP

#include <Rcpp.h>
#include <iostream>
#include <fstream>
using namespace Rcpp;

static std::string base64_LUT = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

// check platform endian
std::string hpp_getEndian () {
  std::string out = "";
  unsigned int foo = 1;
  char *bar = (char*)&foo;
  switch(*bar) {
  case 0: {
    out = "big";
  }
    break;
  case 1: {
    out = "little";
  }
    break;
  }
  return out;
}

// combine Rcpp Vectors
template <int RTYPE> Rcpp::Vector<RTYPE>
c_vector_T( const Rcpp::Vector<RTYPE>&x, const Rcpp::Vector<RTYPE>&y) {
  Rcpp::Vector<RTYPE> out(x.size() + y.size());
  std::copy(x.begin(), x.end(), out.begin());
  std::copy(y.begin(), y.end(), out.begin() + x.size());
  return(out);
}
SEXP c_vector( SEXP x, SEXP y) {
  if((TYPEOF(x) == NILSXP) && (TYPEOF(y) == NILSXP)) return R_NilValue;
  if(TYPEOF(x) == NILSXP) return y;
  if(TYPEOF(y) == NILSXP) return x;
  if(TYPEOF(x) != TYPEOF(y)) Rcpp::stop("c_vector: 'x' and 'y' should have same RTYPE");
  switch( TYPEOF(x) ) {
  case NILSXP: return R_NilValue;
  case LGLSXP: return c_vector_T<LGLSXP>(x, y);
  case INTSXP: return c_vector_T<INTSXP>(x, y);
  case REALSXP: return c_vector_T<REALSXP>(x, y);
  case STRSXP: return c_vector_T<STRSXP>(x, y);
  case RAWSXP: return c_vector_T<RAWSXP>(x, y);
  default: Rcpp::stop("c_vector: not supported type in 'x'");
  }
}

// Ensures NumericVector is not NULL
bool nNotisNULL(const Rcpp::Nullable<Rcpp::NumericVector> x_ = R_NilValue) {
  if (x_.isNotNull()) {
    Rcpp::NumericVector x(x_.get());
    if(x.size() > 0) {
      return true;
    }
  }
  return false;
}

// Ensures IntegerVector is not NULL
bool iNotisNULL(const Rcpp::Nullable<Rcpp::IntegerVector> x_ = R_NilValue) {
  if (x_.isNotNull()) {
    Rcpp::IntegerVector x(x_.get());
    if(x.size() > 0) {
      return true;
    }
  }
  return false;
}

// Ensures LogicaVector is not NULL
bool lNotisNULL(const Rcpp::Nullable<Rcpp::LogicalVector> x_ = R_NilValue) {
  if(x_.isNotNull()) {
    Rcpp::LogicalVector x(x_.get());
    if(x.size() > 0) {
      return true;
    }
  }
  return false;
}

//' @title Multiple Pattern Fixed Matching
//' @name cpp_mpfmatch
//' @description
//' String matching with multiple pattern.
//' @param x,pattern Nullable Rcpp CharacterVector.
//' @details equivalent of as.logical(sum(unlist(lapply(pattern, grepl, x = x, fixed = TRUE)))).
//' @return a bool
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
bool hpp_mpfmatch(const Rcpp::Nullable<Rcpp::CharacterVector> x = R_NilValue,
                  const Rcpp::Nullable<Rcpp::CharacterVector> pattern = R_NilValue) {
  bool found = false;
  if(x.isNotNull()) {
    Rcpp::CharacterVector xx(x);
    if(xx.size() > 0) {
      std::string str = std::string(xx[0]);
      if(pattern.isNotNull()) {
        Rcpp::CharacterVector pat(pattern.get());
        for(R_len_t i = 0; i < pat.size(); i++) {
          if((i % 100) == 0) Rcpp::checkUserInterrupt();
          std::string p = std::string(pat[i]);
          if(str.find(p) != std::string::npos) {
            found = true;
            break;
          }
        }
      }
    }
  }
  return found;
}

//' @title Sequence of Strings Matching
//' @name cpp_seqmatch
//' @description
//' Match a sequence of strings in another string
//' @param x,y StringVector to match
//' @details smallest sequence will be searched into the largest one.
//' @return the index (starting at 1) when a match has been found. Otherwise 0.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
R_len_t hpp_seqmatch(const Rcpp::StringVector x,
                     const Rcpp::StringVector y) {
  R_len_t i = 0, j = 0, k = 0;
  if(x.size() < y.size()) {
    return hpp_seqmatch(y, x);
  } else {
    while((k < y.size()) && (i < x.size())) {
      Rcpp::checkUserInterrupt();
      j = 0;
      k = 0;
      while((i < x.size()) && (j < y.size())) {
        // Rcpp::checkUserInterrupt();
        if(x[i] == y[j]) {
          k++;
          i++;
        }
        j++;
      }
      if(k == 0) i++;
    }
  }
  return (k == y.size()) ? (1 + i - y.size()) : 0;
}

//' @title Use Rcpp to Apply Any on Matrix Rows
//' @name cpp_fast_rowAny
//' @description
//' Computes any across matrix rows
//' @param M_ a Nullable LogicalVector. /!\ But cast to LogicalMatrix
//' @return a LogicalVector.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::LogicalVector> hpp_fast_rowAny(const Rcpp::Nullable<Rcpp::LogicalVector> M_ = R_NilValue) {
  if(!lNotisNULL(M_)) return M_;
  Rcpp::LogicalVector V(M_.get());
  if(V.hasAttribute("dim")) {
    Rcpp::IntegerVector d = V.attr("dim");
    if(d.size() == 2) {
      Rcpp::LogicalMatrix M(M_.get());
      Rcpp::LogicalVector out = M(Rcpp::_, 0);
      for(R_len_t i_col = 1; i_col < M.ncol(); i_col++) out = out | M(Rcpp::_, i_col);
      return out; 
    } else {
      Rcpp::stop("hpp_fast_rowAny: input is not coercible to logical matrix");
    }
  }
  return V;
}

//' @title Use Rcpp to Apply Any on List members
//' @name cpp_fast_listAny
//' @description
//' Computes any across list members
//' @param L_ a Nullable List.
//' @return a LogicalVector.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::LogicalVector> hpp_fast_listAny(const Rcpp::Nullable<Rcpp::List> L_ = R_NilValue) {
  if(L_.isNotNull()) {
    Rcpp::List L(L_.get());
    Rcpp::Nullable<Rcpp::LogicalVector> bar = L[0];
    if(lNotisNULL(bar)) {
      Rcpp::LogicalVector out(Rcpp::clone(bar.get()));
      for(R_len_t i = 1; i < L.size(); i++) {
        Rcpp::LogicalVector foo = L[i];
        if(out.size() != foo.size()) Rcpp::stop("hpp_fast_listAny: members of 'L' should have same length");
        out = out | foo;
      }
      return out;
    } else {
      for(R_len_t i = 1; i < L.size(); i++) {
        Rcpp::Nullable<Rcpp::LogicalVector> bar = L[i];
        if(lNotisNULL(bar)) Rcpp::stop("hpp_fast_listAny: members of 'L' should have same length");
      }
    }
  }
  return R_NilValue;
}

//' @title Use Rcpp for Range
//' @name cpp_fast_range
//' @description
//' Determines range of numeric vector
//' @param x_ a Nullable NumericVector.
//' @details the behaviour is the same as R base::range(x_, na.rm = TRUE, finite = TRUE) without creating warnings
//' @return a NumericVector.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_fast_range(const Rcpp::Nullable<Rcpp::NumericVector> x_ = R_NilValue) {
  double xmax = R_NegInf, xmin = R_PosInf;
  if(nNotisNULL(x_)) {
    Rcpp::NumericVector x(x_.get());
    for(R_len_t i = 0; i < x.size(); i++) {
      if((x[i] != NA_REAL) && (x[i] != R_NaN) && (x[i] != R_NegInf) && (x[i] != R_PosInf)) {
        if(x[i] > xmax) xmax = x[i];
        if(x[i] < xmin) xmin = x[i]; 
      }
    }
  }
  if(xmin > xmax) return Rcpp::NumericVector::create(xmax, xmin);
  return Rcpp::NumericVector::create(xmin, xmax);
}

//' @title Use Rcpp for Sampling
//' @name cpp_fast_sample
//' @description
//' Create a sample of integers
//' @param n a R_len_t, max number integers to choose from.
//' @param size a R_len_t the desired size of return integers.
//' @param replace a bool determining if sampling should be done with replacement. Default is false.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerVector hpp_fast_sample(const R_len_t n = 0,
                                    const R_len_t size = 0,
                                    const bool replace = false) {
  return Rcpp::sample(n, size, replace);
}

//' @title Get Bytes Order
//' @name cpp_get_bytes_order
//' @description
//' This function expands bytes order to the whole data
//' @param obj number of objects in the data.
//' @param byt_ IntegerVector of number of bytes to take from 'ord_'.
//' @param ord_ IntegerVector bytes order.
//' @param rev bool whether to reverse order or not. Default is false.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::NumericVector> hpp_get_bytes_order (const R_len_t obj = 0,
                                                         const Rcpp::Nullable<Rcpp::IntegerVector> byt_ = R_NilValue,
                                                         const Rcpp::Nullable<Rcpp::IntegerVector> ord_ = R_NilValue,
                                                         const bool rev = false) {
  if(iNotisNULL(byt_) && iNotisNULL(ord_) && (obj > 0)) {
    Rcpp::IntegerVector byt(byt_.get());
    Rcpp::IntegerVector ord(ord_.get());
    Rcpp::IntegerVector alw = Rcpp::IntegerVector::create(1,2,3,4,5,6,7,8);
    R_len_t n = 0;
    Rcpp::List L(byt.size());
    for(R_len_t i = 0; i < byt.size(); i++) {
      if(!is_true(any(byt[i] == alw))) Rcpp::stop("hpp_get_bytes_order: 'byt' contains not allowed value:");
      Rcpp::IntegerVector v = seq_len(byt[i]);
      v = na_omit(match(ord, v));
      v = rep_len(v, byt[i]);
      if(rev) {
        L[i] = Rcpp::rev(v);
      } else {
        L[i] = v;
      }
      n += byt[i];
    }
    Rcpp::NumericVector out = Rcpp::no_init(obj * n);
    R_len_t k = 0;
    for(R_len_t g = 0; g < obj; g++) {
      for(R_len_t i = 0; i < byt.size(); i++) {
        Rcpp::IntegerVector v = L(i);
        for(R_len_t j = 0; j < v.size(); j++) {
          out[k + j] = v[j] + k;
        }
        k += byt[i];
      }
    }
    return out;
  }
  return R_NilValue;
}

// converts unsigned short to string
std::string to_string(const double x) {
  std::string out;
  std::ostringstream convert;
  convert << x;
  out = convert.str();
  return out;
}

//' @title Non Finite Values Replacement
//' @name cpp_replace_non_finite
//' @description
//' This function replaces non finite values (NA, NaN -Inf and +Inf)
//' @param V a NumericVector.
//' @param by a double used as replacement value. Default is 0.0
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::NumericVector> hpp_replace_non_finite(const Rcpp::Nullable<Rcpp::NumericVector> V_ = R_NilValue,
                                                           const double by = 0.0) {
  if(nNotisNULL(V_)) {
    Rcpp::NumericVector V(V_.get());
    Rcpp::NumericVector out = Rcpp::no_init(V.size()); // new
    Rcpp::LogicalVector a = is_infinite(V);
    Rcpp::LogicalVector b = is_nan(V);
    for(R_len_t i = 0; i < V.size(); i++) {
      out[i] = (a[i] | b[i]) ? by : V[i];
    }
    return out;
  }
  return V_;
}

// determines range of a numeric vector AND ensures that it is of finite values.
// minimal value will be clipped to -4095.0
Rcpp::NumericVector hpp_check_range(const Rcpp::NumericVector x) {
  double Min = R_PosInf, Max = R_NegInf;
  if(nNotisNULL(x)) {
    Rcpp::LogicalVector a = is_infinite(x);
    Rcpp::LogicalVector b = is_nan(x);
    for(R_len_t i = 0; i < x.size(); i++) {
      if(a[i] || b[i]) Rcpp::stop("hpp_check_range: 'x' contains non-finite values");
      if((x[i] < Min) && (x[i] > -4095.0)) Min = x[i];
      if(x[i] > Max) Max = x[i];
    }
    if(Min == R_PosInf) Min = Max;
  } else {
    Rcpp::stop("hpp_check_range: 'x' is empty");
  }
  if(Min > Max) return Rcpp::NumericVector::create(Max, Min); // New
  return Rcpp::NumericVector::create(Min, Max);
}

// Thanks to 
// http://www.libpng.org/pub/png/book/chapter10.html
// image_sample = light_out ^ gamma
// it is said that
// Once again, bear in mind that light_out and image_sample are scaled to the interval between 0 and 1;
// that is, if the sample depth is 8 bits, the file samples range between 0 and 255, so image_sample is
// obtained by dividing a given file sample by 255, in floating-point arithmetic. 
// So,
// image_sample = ymid and its range is [0,255]
// light_out = xmid and its range is [xmin,xmax]
// we have ymid / 255 = ((xmid - xmin) / (xmax - xmin)) ^ gamma
// log(ymid / 255) = log((xmid - xmin) / (xmax - xmin)) * gamma
// gamma = log(ymid / 255) / log((xmid - xmin) / (xmax - xmin))
//' @title Gamma Computation
//' @name cpp_computeGamma
//' @description
//' This function computes image gamma transformation value.
//' @param V named NumericVector of channel display properties containing 'xmin', 'xmax', 'xmid' and 'ymid'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
double hpp_computeGamma (const Rcpp::NumericVector V) {
  double V_w = V["xmax"], V_m = V["xmid"], V_y = V["ymid"];
  V_w -= V["xmin"];
  V_m -= V["xmin"];
  return (( std::log(V_y / 255) ) / (std::log(V_m / V_w)));
}

//' @title Raw to Base64 Conversion
//' @name cpp_base64_encode
//' @description
//' Converts a raw vector to base64 string.
//' @param x RawVector.
//' @return a string, representing the base64 encoding of x.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
std::string hpp_base64_encode (const Rcpp::RawVector x) {
  R_len_t i = 0, a = x.size() / 3, b = x.size() % 3;
  std::string out;
  out.reserve(((a) + (b > 0)) * 4);
  for(R_len_t idx = 0; idx < a; idx++, i += 3) { 
    uint32_t val = (x[i] << 16) + (x[i + 1] << 8) + x[i + 2];
    for(short k = 18; k >= 0; k -= 6) {
      out.push_back(base64_LUT[(val >> k) & 0x3F]);
    }
  }
  switch(b) {
  case 1: {
    uint32_t val  = x[i++] << 16;
    out.push_back(base64_LUT[(val >> 18) & 0x3F]);
    out.push_back(base64_LUT[(val >> 12) & 0x3F]);
    out.append(2,'=');
    break;
  }
  case 2: {
    uint32_t val  = (x[i] << 16) + (x[i + 1] << 8);
    out.push_back(base64_LUT[(val >> 18) & 0x3F]);
    out.push_back(base64_LUT[(val >> 12) & 0x3F]);
    out.push_back(base64_LUT[(val >>  6) & 0x3F]);
    out.push_back('=');
    break;
  }
  }
  return out;
}

//' @title BMP Writer
//' @name cpp_writeBMP
//' @description
//' Transforms 3D [0,1] image to uncompressed bmp
//' @param image, a [0,1] normalized image matrix or 3D array. If 3D array, 3rd dimension should be of length 1 or 3.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_writeBMP (const Rcpp::NumericVector image) {
  if(nNotisNULL(image)) {
    Rcpp::IntegerVector d = image.attr("dim");
    if(iNotisNULL(d)) {
      bool rgb = false;
      if(!(d.size() == 2 || d.size() == 3)) {
        Rcpp::stop("hpp_writeBMP: image should be a matrix or a 3D array");
      } else {
        if(d.size() == 3) {
          if(d[2] == 3) rgb = true;
          if(!(d[2] == 1 || rgb)) {
            Rcpp::stop("hpp_writeBMP: when 3D array is provided, 3rd dim should be 1 or 3");
          }
        }
      }
      // compute padding
      R_len_t padding = d[1] % 4;
      // compute data size from input
      R_len_t dsize = rgb ? image.size() : 3 * image.size();
      // compute extra size required for padding
      R_len_t extra = dsize + d[0] * (d[1] + padding);
      // compute final file size
      R_len_t fsize = extra + 54;
      // declare out;
      Rcpp::RawVector out = Rcpp::no_init(fsize);
      
      ////////////////// Bitmap file header (=14 Bytes)
      // 'B' 'M' header
      out[ 0] = 0x42, out[ 1] = 0x4d; 
      // final file size
      out[ 2] = fsize & 0xff;
      out[ 3] = (fsize >>  8) & 0xff;
      out[ 4] = (fsize >> 16) & 0xff;
      out[ 5] = (fsize >> 24) & 0xff;
      // app reserved
      out[ 6] = 0x00, out[ 7] = 0x00;
      out[ 8] = 0x00, out[ 9] = 0x00;
      // offset of data (14 + 40 = 54 Bytes)
      out[10] = 0x36, out[11] = 0x00, out[12] = 0x00, out[13] = 0x00;
      
      //////////////////  DIB header (bitmap information header) (=40 Bytes)
      // size = DIB (=40 Bytes)
      out[14] = 0x28, out[15] = 0x00, out[16] = 0x00, out[17] = 0x00;
      // image width
      out[18] = d[1] & 0xff;
      out[19] = (d[1] >>  8) & 0xff;
      out[20] = (d[1] >> 16) & 0xff;
      out[21] = (d[1] >> 24) & 0xff;
      // image height
      out[22] = d[0] & 0xff;
      out[23] = (d[0] >>  8) & 0xff;
      out[24] = (d[0] >> 16) & 0xff;
      out[25] = (d[0] >> 24) & 0xff;
      // color plane (=1)
      out[26] = 0x01, out[27] = 0x00;
      // bits per pixel (8 b + 8 g + 8 r = 24), no alpha
      out[28] = 0x18, out[29] = 0x00;
      // compression used (=0), no compression
      out[30] = 0x00, out[31] = 0x00, out[32] = 0x00, out[33] = 0x00;
      // image data size, including extra padding
      out[34] = extra & 0xff;
      out[35] = (extra >>  8) & 0xff;
      out[36] = (extra >> 16) & 0xff;
      out[37] = (extra >> 24) & 0xff;
      // horizontal print resolution 96 DPI x 39.3701 inches per metre yields 3779.53
      out[38] = 0xc3, out[39] = 0x0e, out[40] = 0x00, out[41] = 0x00;
      // vertical print resolution 96 DPI x 39.3701 inches per metre yields 3779.53
      out[42] = 0xc3, out[43] = 0x0e, out[44] = 0x00, out[45] = 0x00;
      // number of color palette (=0)
      out[46] = 0x00, out[47] = 0x00, out[48] = 0x00, out[49] = 0x00;
      // important colors (=0), 0 means all colors are important 
      out[50] = 0x00, out[51] = 0x00, out[52] = 0x00, out[53] = 0x00;
      
      //////////////////  copy image to out
      R_len_t n = 54;
      if(rgb) {
        for(R_len_t i = d[0] - 1; i >= 0; i--) { // image in BMP is height inverted
          for(R_len_t j = 0; j < d[1]; j++) {
            for(R_len_t k = d[2] - 1; k >= 0; k--) { // image in BMP is rgb inverted
              out[n++] = 255 * image[k * d[1] * d[0] + j * d[0] + i];
            }
          }
          // apply padding
          for(R_len_t p = 0; p < padding; p++) {
            out[n++] = 0;
          }
        }
      } else {
        for(R_len_t i = d[0] - 1; i >= 0; i--) { // image in BMP is height inverted
          for(R_len_t j = 0; j < d[1]; j++) {
            for(R_len_t k = 0; k < 3; k++) { // write 3 times the same value
              out[n++] = 255 * image[j * d[0] + i];
            }
          }
          // apply padding
          for(R_len_t p = 0; p < padding; p++) {
            out[n++] = 0;
          }
        }
      }
      return out;
    } else {
      Rcpp::stop("hpp_writeBMP: image should be a matrix or a 3D array");
    }
  }
  return R_NilValue;
}

//' @title File Raw Chunk Extraction
//' @name hpp_readchunk
//' @description
//' Reads binary chunk from file.
//' @param fname string, path to file.
//' @param offset std::size_t, where to start reading from the beginning.
//' @param nbytes uint32_t, number of bytes.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @return a RawVector
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_readchunk (const std::string fname, 
                               const std::size_t offset,
                               const uint32_t nbytes,
                               const bool verbose = false) {
  if(verbose) {
    Rcout << fname << std::endl;
    Rcout << "Extracting " << nbytes << " Bytes @offset:" << offset << std::endl;
  }
  FILE *fi = std::fopen(fname.c_str(), "rb");
  if(fi == nullptr) Rcpp::stop("hpp_readchunk: Unable to open file");
  std::fseek(fi, 0, SEEK_END);
  std::size_t filesize = std::ftell(fi);
  if((filesize - offset) < nbytes) {
    std::fclose(fi);
    Rcpp::Rcerr << "hpp_readchunk: larger nbytes [" << nbytes << "] to read than remaining filesize - offset ["<< (filesize - offset) << "]\n" << fname  << std::endl;
    Rcpp::stop("hpp_readchunk: can't read so much bytes");
  }
  if((filesize - offset) > filesize) {
    std::fclose(fi);
    Rcpp::Rcerr << "hpp_readchunk: @offset:" << offset << " points to outside of\n" << fname  << std::endl;
    Rcpp::stop("hpp_readchunk: offset is higher than file size");
  }
  Rcpp::RawVector out = Rcpp::no_init_vector(nbytes);
  std::fseek(fi, offset, SEEK_SET);
  auto s = std::fread(&out[0], sizeof(char), nbytes, fi);
  std::fclose(fi);
  if(s != nbytes) Rcpp::stop("hpp_readchunk: bad read");
  return out;
}

#endif
