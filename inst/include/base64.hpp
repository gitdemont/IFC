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
  if(x.length() == 0) return 0;
  Rcpp::RawVector out = Rcpp::no_init_vector(3 * x.length() / 4);
  
  // create LUT for conversion
  Rcpp::RawVector LUT_RAW(256);
  LUT_RAW.fill(64);
  std::string LUT = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
  std::string::iterator it_lut = LUT.begin();
  while(it_lut != LUT.end()) { 
    unsigned char v = *(it_lut++);
    LUT_RAW[v] = 65 <= v && v <= 90 ? v - 65 : 97 <= v && v <= 122 ? v - 71 : 48 <= v && v <= 57 ? v + 4 : 64;
  }
  if(url) {
    LUT_RAW[45] = 62; // -
    LUT_RAW[95] = 63; // _
  } else {
    LUT_RAW[43] = 62; // +
    LUT_RAW[47] = 63; // /
  }
  
  // convert from base64
  std::string::iterator it = x.begin();
  R_len_t i = 0, j = 0;
  short pos = 1;
  while((pos == 1) && (it != x.end())) {
    j += 4;
    unsigned char a = LUT_RAW[*(it++)];
    unsigned char b = LUT_RAW[*(it++)];
    if(a >= 64 || b >= 64) pos = 4;
    unsigned char c = LUT_RAW[*(it++)];
    if(c >= 64) pos = 3;
    unsigned char d = LUT_RAW[*(it++)];
    if(d >= 64) pos = 2;
    out[i++] = (a << 2) + ((b & 0xF0) >> 4);
    out[i++] = ((b & 0x0F) << 4) + ((c & 0x3C) >> 2);
    out[i++] = ((c & 0x03) << 6) + d;
  }
  
  // inspect and trim last bytes
  std::string msg = "";
  if(it != x.end()) msg.append("\n-premature ending");
  i -= pos;
  if(i <= 0) return 0;
  Rcpp::RawVector V = Rcpp::no_init_vector(4);
  short bad = 0;
  for(short k = 0; k < V.size(); k++) { 
    V[k] = (unsigned char) *(--it);
    if(!bad && V[k] != 61 && LUT_RAW[V[k]] >= 64) {bad = k + 1; break;}
  }
  if(!((V[3] != 61 && V[2] != 61) &&
     ((V[1] == 61 && V[0] == 61) || (V[1] != 61)))) msg.append("\n-invalid padding");
  if(bad) {
    std::stringstream s;
    s << std::uppercase << "0x" << std::hex << (unsigned short) V[bad - 1];
    msg.append("\n-invalid base64 character[");
    msg.append(s.str());
    msg.append("] @");
    msg.append(std::to_string(j - (bad - 1)));
  }
  if(msg != "") Rcpp::warning(msg);
  if(pos == 1) return out;
  return out[Rcpp::Range(0, i)];
}

#endif
