/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2024 Yohann Demont                                              
                                                                                
  It is part of IFC package, please cite:                                       
  -IFC: An R Package for Imaging Flow Cytometry                                 
  -YEAR: 2024                                                                   
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

#ifndef IFC_BASE64_HPP
#define IFC_BASE64_HPP

#include <Rcpp.h>
using namespace Rcpp;

// from Dirk Eddelbuettel answer on stackoverflow
// https://stackoverflow.com/a/73822223
std::string r2c (Rcpp::RawVector x) {
  const char* c = reinterpret_cast<char*>(x.begin());
  return std::string(c);
}

//' @title Raw to Base64 Conversion
//' @name cpp_base64_encode
//' @description
//' Converts a raw vector to base64 string.
//' @param x RawVector.
//' @param url a bool, whether to convert for url. Default is false.
//' @return a string, representing the base64 encoding of x.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
std::string hpp_base64_encode (const Rcpp::RawVector x, const bool url = false) {
  R_len_t i = 0, j = 0, a = x.size() / 3, b = x.size() % 3;
  Rcpp::RawVector out = Rcpp::no_init_vector(1 + ((a) + (b > 0)) * 4); //+ 1 to add 00 at the end
  const std::string LUT = url ? "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_" : "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  for(R_len_t idx = 0; idx < a; idx++, i += 3) { 
    uint32_t val = (x[i] << 16) + (x[i + 1] << 8) + x[i + 2];
    out[j++] = LUT[(val >> 18) & 0x3F];
    out[j++] = LUT[(val >> 12) & 0x3F];
    out[j++] = LUT[(val >>  6) & 0x3F];
    out[j++] = LUT[ val        & 0x3F];
  }
  switch(b) {
  case 1: {
    uint32_t val = x[i] << 16;
    out[j++] = LUT[(val >> 18) & 0x3F];
    out[j++] = LUT[(val >> 12) & 0x3F];
    out[j++] = 0x3D;
    out[j++] = 0x3D;
    break;
  }
  case 2: {
    uint32_t val = (x[i] << 16) + (x[i + 1] << 8);
    out[j++] = LUT[(val >> 18) & 0x3F];
    out[j++] = LUT[(val >> 12) & 0x3F];
    out[j++] = LUT[(val >>  6) & 0x3F];
    out[j++] = 0x3D;
    break;
  }
  }
  out[j] = 0x00; // add final 00 to terminate string in r2c conversion
  return r2c(out);
}

// helpers for base64 decoding
char convb64(const char x, char *abort) {
  if(65 <= x && x <= 90) return x - 65;
  if(97 <= x && x <= 122) return x - 71;
  if(48 <= x && x <= 57) return x + 4;
  if(x == 43) return 62;
  if(x == 47) return 63;
  if(*abort == 0) *abort = x;
  return 0;
}

char convb64_url(const char x, char *abort) {
  if(65 <= x && x <= 90) return x - 65;
  if(97 <= x && x <= 122) return x - 71;
  if(48 <= x && x <= 57) return x + 4;
  if(x == 45) return 62;
  if(x == 95) return 63;
  if(*abort == 0) *abort = x;
  return 0;
}

//' @title Base64 to Raw Conversion
//' @name cpp_base64_decode
//' @description
//' Converts a base64 string to raw vector.
//' @param x a string.
//' @param url a bool, whether to convert for url. Default is false.
//' @return a RawVector, representing the decoding of base64 string.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_base64_decode(std::string x, const bool url = false) {
  if((x.length() % 4) != 0) Rcpp::stop("hpp_base64_decode: not a base64 encoded string");
  Rcpp::RawVector out = Rcpp::no_init_vector(3 * x.length() / 4);
  if(x.length() == 0) return out;
  std::string::iterator it = x.begin();
  R_len_t i = 0;
  char abort = 0;
  if(url) {
    while(!abort && (it != x.end())) {
      out[i] = convb64_url(*(it++), &abort) << 2; 
      char b = convb64_url(*(it++), &abort);
      char c = convb64_url(*(it++), &abort);
      out[i++] += ((b & 0xF0) >> 4);
      out[i++] = ((b & 0x0F) << 4) + ((c & 0x3C) >> 2);
      out[i++] = ((c & 0x03) << 6) + convb64_url(*(it++), &abort);
    }
  } else {
    while(!abort && (it != x.end())) {
      out[i] = convb64(*(it++), &abort) << 2; 
      char b = convb64(*(it++), &abort);
      char c = convb64(*(it++), &abort);
      out[i++] += ((b & 0xF0) >> 4);
      out[i++] = ((b & 0x0F) << 4) + ((c & 0x3C) >> 2);
      out[i++] = ((c & 0x03) << 6) + convb64(*(it++), &abort);
    }
  }
  // inspect and trim last bytes
  std::string msg = "";
  if(it != x.end()) msg.append("\n-premature ending");
  Rcpp::RawVector V = Rcpp::no_init_vector(4);
  R_len_t pos = 1;
  bool bad = false;
  char foo = 0;
  for(short k = 0; k < V.size(); k++) {
    V[k] = *(--it);
    if(V[k] == abort) pos = k + 2;
    if(!bad) bad = (V[k] != 65 && V[k] != 61) && !(url ? convb64_url(V[k], &foo) : convb64(V[k], &foo));
  }
  if(!((V[3] != 61 && V[2] != 61) &&
     ((V[1] == 61 && V[0] == 61) || (V[1] != 61)))) msg.append("\n-invalid padding");
  if(bad) msg.append("\n-invalid base64 character");
  if(msg != "") Rcpp::warning(msg);
  i -= pos;
  if(i < 0) return 0;
  return out[Rcpp::Range(0, i)];
}

#endif
