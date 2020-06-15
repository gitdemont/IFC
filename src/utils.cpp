/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2020 Yohann Demont                                              
                                                                                
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

#include <Rcpp.h>
using namespace Rcpp;

// TODO: template to test if Rcpp Vector is not NULL 
// template <typename T> bool NotisNULL(Rcpp::Nullable<T> x_ = R_NilValue) {
//   if (x_.isNotNull()) {
//     T x(x_.get());
//     if(x.size() > 0) {
//       return true;
//     }
//   }
//   return false;
// }
static std::string base64_LUT = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

// Ensures NumericVector is not NULL
bool nNotisNULL(Nullable<NumericVector> x_ = R_NilValue) {
  if (x_.isNotNull()) {
    NumericVector x(x_.get());
    if(x.size() > 0) {
      return true;
    }
  }
  return false;
}

// Ensures IntegerVector is not NULL
bool iNotisNULL(Nullable<IntegerVector> x_ = R_NilValue) {
  if (x_.isNotNull()) {
    IntegerVector x(x_.get());
    if(x.size() > 0) {
      return true;
    }
  }
  return false;
}

// converts unsigned short to string
std::string to_string(uint16_t x) {
  std::string out;
  std::ostringstream convert;
  convert << x;
  out = convert.str();
  return out;
}

// determines range of a numeric vector AND ensures that it is of finite values.
// minimal value will be clipped to -4095.0
NumericVector cpp_check_range(const NumericVector x) {
  double Min = R_PosInf, Max = R_NegInf;
  if(nNotisNULL(x)) {
    for(R_len_t i = 0; i < x.size(); i++) {
      if(!Rcpp::traits::is_finite<REALSXP>(x[i])) Rcpp::stop("cpp_check_range: 'x' contains non-finite values");
      if((x[i] < Min) && (x[i] > -4095.0)) Min = x[i];
      if(x[i] > Max) Max = x[i];
    }
  } else {
    Rcpp::stop("cpp_check_range: 'x' is empty");
  }
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
// [[Rcpp::export]]
double cpp_computeGamma (const NumericVector V) {
  double V_w = V["xmax"], V_m = V["xmid"], V_y = V["ymid"];
  V_w -= V["xmin"];
  V_m -= V["xmin"];
  return (( log(V_y / 255) ) / (log(V_m / V_w)));
}

//' @title Raw to Base64 Conversion
//' @name cpp_base64_encode
//' @description
//' Converts a raw vector to base64 string.
//' @param x RawVector.
//' @return a string, representing the base64 encoding of x.
//' @keywords internal
////' @export
// [[Rcpp::export]]
std::string cpp_base64_encode ( const RawVector x) {
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
// [[Rcpp::export]]
RawVector cpp_writeBMP ( const NumericVector image) {
  if(nNotisNULL(image)) {
    IntegerVector d = image.attr("dim");
    if(iNotisNULL(d)) {
      bool rgb = false;
      if(!(d.size() == 2 || d.size() == 3)) {
        Rcerr << "cpp_writeBMP: image should be a matrix or a 3D array" << std::endl;
        Rcpp::stop("cpp_writeBMP: image should be a matrix or a 3D array");
      } else {
        if(d.size() == 3) {
          if(d[2] == 3) rgb = true;
          if(!(d[2] == 1 || rgb)) {
            Rcerr << "cpp_writeBMP: when 3D array is provided, 3rd dim should be 1 or 3" << std::endl;
            Rcpp::stop("cpp_writeBMP: when 3D array is provided, 3rd dim should be 1 or 3");
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
      RawVector out = no_init(fsize);
      
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
      R_len_t i, j, k, p;
      if(rgb) {
        for(i = d[0] -1; i >=0; i--) { // image in BMP is height inverted
          for(j = 0; j < d[1]; j++) {
            for(k = d[2] - 1; k >= 0; k--) { // image in BMP is rgb inverted
              out[n++] = 255 * image[k * d[1] * d[0] + j * d[0] + i];
            }
          }
          // apply padding
          for(p = 0; p < padding; p++) {
            out[n++] = 0;
          }
        }
      } else {
        for(i = d[0] -1; i >=0; i--) { // image in BMP is height inverted
          for(j = 0; j < d[1]; j++) {
            for(k = 0; k < 3; k++) { // write 3 times the same value
              out[n++] = 255 * image[j*d[0] + i];
            }
          }
          // apply padding
          for(p = 0; p < padding; p++) {
            out[n++] = 0;
          }
        }
      }
      return out;
    } else {
      Rcerr << "cpp_writeBMP: image should be a matrix or a 3D array" << std::endl;
      Rcpp::stop("cpp_writeBMP: image should be a matrix or a 3D array");
    }
  }
  return R_NilValue;
}
