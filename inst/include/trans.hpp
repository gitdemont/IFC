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

#ifndef IFC_TRANS_HPP
#define IFC_TRANS_HPP

#include <Rcpp.h>
#include "utils.hpp"
using namespace Rcpp;

//' @title Smooth LinLog Transformation with Rcpp
//' @name cpp_smoothLinLog
//' @description
//' Takes a numeric vector and return its transformation:
//' - to linear, if abs(x) < hyper.
//' - to log, if abs(x) > hyper.
//' @param x NumericVector.
//' @param hyper double, value where transition between Lin/Log is applied.
//' @param base double, base of Log scale.
//' @param lin_comp double, value that is used to smooth transition between Lin/Log.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_smoothLinLog (const Rcpp::NumericVector x,
                                      const double hyper = 1000.0,
                                      const double base = 10.0,
                                      const double lin_comp = 2.302585) {
  if(nNotisNULL(x)) {
    R_len_t L = x.size();
    double K = std::log(base) / lin_comp;
    Rcpp::NumericVector B = K + K * Rcpp::log((Rcpp::abs(x))/hyper);
    Rcpp::NumericVector xx = Rcpp::no_init_vector(L);
    for(R_len_t i = 0;i < L; i++) {
      if( std::fabs(x(i)) <= hyper ) {
        xx(i) = x(i) * K / hyper;
      } else {
        xx(i) = B(i);
        if(x(i)<0) xx(i) *= -1;
      }
    }
    return xx;
  } else {
    return x;
  }
}

//' @title Inverse Smooth LinLog Transformation with Rcpp
//' @name cpp_inv_smoothLinLog
//' @description
//' Takes a numeric vector and return its transformation:
//' - to linear, if abs(x) < log(base) / lin_comp.
//' - to exp, if abs(x) > log(base) / lin_comp.
//' @param x NumericVector.
//' @param hyper double, value where transition between Lin/Log is applied.
//' @param base double, base of Log scale.
//' @param lin_comp double, value that is used to smooth transition between Lin/Log.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_inv_smoothLinLog (const Rcpp::NumericVector x,
                                          const double hyper = 1000.0, 
                                          const double base = 10.0, 
                                          const double lin_comp = 2.302585) {
  if(nNotisNULL(x)) {
    R_len_t L = x.size();
    double K = std::log(base) / lin_comp;
    Rcpp::NumericVector B = hyper * Rcpp::exp((Rcpp::abs(x))/K-1);
    Rcpp::NumericVector xx = Rcpp::no_init_vector(L);
    for(R_len_t i = 0;i < L; i++) {
      if( std::fabs(x(i)) > K ) {
        xx(i) = B(i);
        if(x(i)<0) xx(i) *= -1;
      } else {
        xx(i) = x(i) * hyper / K;
      }
    }
    return xx;
  } else {
    return x;
  }
}

//' @title Uint32 to Raw Conversion
//' @name cpp_uint32_to_raw
//' @description
//' Converts unsigned 32bits integer to raw
//' @param x uint32_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_uint32_to_raw(const uint32_t x) {
  Rcpp::RawVector out(4);
  out[3] = (x >> 24) & 0xff;
  out[2] = (x >> 16) & 0xff;
  out[1] = (x >>  8) & 0xff;
  out[0] = (x      ) & 0xff;
  return out;
}

//' @title Int32 to Uint32 32bits Conversion
//' @name cpp_int32_to_uint32
//' @description
//' Converts 32bits integer from signed to unsigned
//' @param x int32_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
uint32_t hpp_int32_to_uint32 (const int32_t x) {
  uint32_t out = x;
  return out;
}

//' @title Uint32 to Int32 32bits Conversion
//' @name cpp_uint32_to_int32
//' @description
//' Converts 32bits integer from unsigned to signed
//' @param x uint32_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
int32_t hpp_uint32_to_int32 (const uint32_t x) {
  int32_t out = x;
  return out;
}

//' @title Offset to Raw Conversion
//' @name cpp_offset_to_raw
//' @description
//' Converts offset to raw
//' @param x double.
//' @param swap bool, whether to swap or not.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_offset_to_raw (const double x, const bool swap = false) {
  uint32_t a = x / 4294967296;
  uint32_t b = x - (a * 4294967296);
  if(swap) {
    Rcpp::RawVector aa = rev(hpp_uint32_to_raw(a));
    Rcpp::RawVector bb = rev(hpp_uint32_to_raw(b));
    return c_vector(bb, aa);
  }
  return c_vector(hpp_uint32_to_raw(b), hpp_uint32_to_raw(a));
}

//' @title Offset to Raw Conversion
//' @name cpp_raw_to_offset
//' @description
//' Converts raw to offset
//' @param x RawVector.
//' @param swap bool, whether to swap or not.
//' @keywords internal
////' @export
// [[Rcpp::export]]
double hpp_raw_to_offset (const Rcpp::RawVector x, const bool swap = false) {
  if(x.size() != 8) Rcpp::stop("cpp_raw_to_offset: 'x' should be a raw vector of length 8");
  if(swap) {
    return hpp_int32_to_uint32((x[3] & 0xff) + 
                               ((x[2] & 0xff) <<  8) +
                               ((x[1] & 0xff) << 16) +
                               ((x[0] & 0xff) << 24)) +
    4294967296 * hpp_int32_to_uint32((x[7] & 0xff) + 
                               ((x[6] & 0xff) <<  8) +
                               ((x[5] & 0xff) << 16) +
                               ((x[4] & 0xff) << 24));
  }
  return hpp_int32_to_uint32((x[0] & 0xff) + 
                             ((x[1] & 0xff) <<  8) +
                             ((x[2] & 0xff) << 16) +
                             ((x[3] & 0xff) << 24)) +
  4294967296 * hpp_int32_to_uint32((x[4] & 0xff) + 
                             ((x[5] & 0xff) <<  8) +
                             ((x[6] & 0xff) << 16) +
                             ((x[7] & 0xff) << 24));
}

//' @title Int64 to Uint64 64bits Conversion
//' @name cpp_int64_to_uint64
//' @description
//' Converts 64bits integer from signed to unsigned
//' @param x int64_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
uint64_t hpp_int64_to_uint64 (const int64_t x) {
  uint64_t out = x;
  return out;
}

//' @title Uint64 to Int64 64bits Conversion
//' @name cpp_uint64_to_int64
//' @description
//' Converts 64bits integer from unsigned to signed
//' @param x uint64_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
int64_t hpp_uint64_to_int64 (const uint64_t x) {
  int64_t out = x;
  return out;
}

//' @title Vectorize Int32 to Uint32 32bits Conversion
//' @name cpp_v_int32_to_uint32
//' @description
//' Converts 32bits vector of integers from unsigned to signed
//' @param V a NumericVector
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::NumericVector> hpp_v_int32_to_uint32 (Rcpp::Nullable<Rcpp::NumericVector> V = R_NilValue) {
  if(nNotisNULL(V)) {
    Rcpp::NumericVector out(V.get());
    for(R_len_t i = 0; i < out.size(); i++) out(i) = hpp_int32_to_uint32(out(i));
    return out;
  }
  return V;
}

//' @title Vectorize Int64 to Uint64 64bits Conversion
//' @name cpp_v_int64_to_uint64
//' @description
//' Converts 64bits vector of integers from unsigned to signed
//' @param V a NumericVector
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::NumericVector> hpp_v_int64_to_uint64 (Rcpp::Nullable<Rcpp::NumericVector> V = R_NilValue) {
  if(nNotisNULL(V)) {
    Rcpp::NumericVector out(V.get());
    for(R_len_t i = 0; i < out.size(); i++) out(i) = hpp_int64_to_uint64(out(i));
    return out;
  }
  return V;
}

#endif
