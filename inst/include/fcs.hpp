/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2024 Yohann Demont                                              
                                                                                
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

#ifndef IFC_FCS_HPP
#define IFC_FCS_HPP

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include "utils.hpp"
using namespace Rcpp;

//' @title FCS data extraction
//' @name cpp_readFCS
//' @description
//' This function reads data from FCS 3.2 files
//' @param fname string, path to file.
//' @param offset std::size_t, offset position of data start
//' @param events uint32_t, number of events to read.
//' @param b_typ Rcpp::IntegerVector, types of values to read. Allowed are 0 for"A", 1 for "F", 2 for "D", and 3 is "I".
//' @param b_siz Rcpp::IntegerVector, number of bytes to extract for each type. Allowed are 0 for 8 1 for 16, 2 for 32 and 3 for 64 bits.\cr
//' Note that whatever the input is, when 'b_typ' is 1 (float) 'b_siz' will be set to 32 and when 'b_typ' is 2 (double) 'b_siz' will be set to 64.
//' @param b_msk Rcpp::IntegerVector, bits to masks when 'b_typ' is 3 (integer). Default is R_NilValue.\cr
//' When not NULL, it should be of same length as 'b_siz' and contain only [0-64] values.
//' @param swap bool, whether to swap bytes or not. Default is false.
//' @source FCS 3.2 specifications. See, J. Spidlen et al. Data File Standard for Flow Cytometry, Version FCS 3.2. Cytometry A 99 100–102(2021) \doi{10.1002/cyto.a.24225}
//' @return a numeric vector of extracted values.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector hpp_readFCS (const std::string fname,
                                 const std::size_t offset,
                                 const uint32_t events,
                                 const Rcpp::IntegerVector b_typ,
                                 const Rcpp::IntegerVector b_siz,
                                 const Rcpp::Nullable<Rcpp::IntegerVector> b_msk = R_NilValue,
                                 const bool swap = false) {
  if(b_typ.size() != b_siz.size()) Rcpp::stop("hpp_readFCS: 'b_typ' and 'b_siz' should be of same length");
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    fi.seekg(0, std::ios::end);
    std::size_t filesize = fi.tellg();
    
    std::size_t count = 0;
    Rcpp::IntegerVector V = Rcpp::no_init_vector(b_typ.size());
    for(R_len_t i = 0; i < b_typ.size(); i++) {
      switch(b_typ[i]) {
      case 0 : { 
        if(b_siz[i] <= 0) {
          Rcpp::Rcerr << "hpp_readFCS: size[" << b_siz[i] << "] not handled for b_typ[A]" << std::endl;
          Rcpp::stop("hpp_readFCS: invalid value in 'b_siz'");
        }
        count += b_siz[i]; V[i] = b_typ[i];
        }
        break;
      case 1 : { count += 4; V[i] = b_typ[i]; }
        break;
      case 2 : { count += 8; V[i] = b_typ[i];  }
        break;
      case 3 : {
        if(b_siz[i] < 0 || b_siz[i] > 3) {
          Rcpp::Rcerr << "hpp_readFCS: size[" << b_siz[i] << "] not handled for b_typ[I]" << std::endl;
          Rcpp::stop("hpp_readFCS: invalid value in 'b_siz'");
        }
        count += std::pow(2.0, b_siz[i]);
        V[i] = b_typ[i] + b_siz[i];
        }
        break;
      default: { Rcpp::stop("hpp_readFCS: 'b_typ' not handled"); }
      }
    }
    
    if(filesize < (events * count + offset)) {
      Rcpp::Rcerr << "hpp_readFCS: can't read " << events * count << " bytes @offset[" << offset << "] with file size[" << filesize << "]" << std::endl;
      Rcpp::stop("hpp_readFCS: trying to read out of file");
    }
    
    Rcpp::NumericVector out = Rcpp::no_init_vector(events * b_typ.size());
    uint8_t  i08;
    uint16_t i16;
    uint32_t i32;
    uint64_t i64;
    float    n32;
    double   n64;
    R_len_t k = 0;
    
    fi.seekg(offset, std::ios::beg);
    if(b_msk.isNotNull()) {
      Rcpp::IntegerVector masks(b_msk.get());
      if(masks.size() != 0) {
        if(b_typ.size() != masks.size()) {
          Rcpp::Rcerr << "hpp_readFCS: when provided 'b_msk' and 'b_siz' should be of same length" << std::endl;
          Rcpp::stop("hpp_readFCS: 'b_msk' and 'b_siz' should be of same size");
        }
        std::vector<uint64_t> m(masks.size());
        for(R_len_t i = 0; i < masks.size(); i++) {
          if((masks[i] < 0) || (masks[i] > 64)) {
            Rcpp::Rcerr << "hpp_readFCS: can't handle bits[" << masks[i] << "] for 'b_msk' should be [0-64]" << std::endl;
            Rcpp::stop("hpp_readFCS: b_msk should be [0-64]");
          }
          m[i] = std::pow(2.0, masks[i]);
          m[i]--;
        }
        
        if(swap) {
          for(uint32_t o = 0; o < events; o++) {
            for(R_len_t i = 0; i < b_typ.size(); i++, k++) {
              switch(V[i]) {
              case 0: {
                //Rcpp::stop("hpp_readFCS: type[A] is not supported in FCS >= 3.2");
                std::vector<char> buf(b_siz[i]);
                fi.read(buf.data(), b_siz[i]);
                std::string s(buf.begin(), buf.end());
                out[k] = std::stod(s.c_str());
              }
                break;
              case 1: { fi.read((char *)&n32,sizeof(n32)); out[k] = bytes_swap(n32); }
                break;
              case 2: { fi.read((char *)&n64,sizeof(n64)); out[k] = bytes_swap(n64); }
                break;
              case 3: { fi.read((char *)&i08,sizeof(i08)); out[k] = bytes_swap(i08) & m[i]; }
                break;
              case 4: { fi.read((char *)&i16,sizeof(i16)); out[k] = bytes_swap(i16) & m[i]; }
                break;
              case 5: { fi.read((char *)&i32,sizeof(i32)); out[k] = bytes_swap(i32) & m[i]; }
                break;
              case 6: { fi.read((char *)&i64,sizeof(i64)); out[k] = bytes_swap(i64) & m[i]; }
                break;
              default: { Rcpp::stop("hpp_readFCS: 'b_typ' and 'b_siz' not handled"); }
              break;
              }
            }
          } 
        } else {
          for(uint32_t o = 0; o < events; o++) {
            for(R_len_t i = 0; i < b_typ.size(); i++, k++) {
              switch(V[i]) {
              case 0: { 
                //Rcpp::stop("hpp_readFCS: type[A] is not supported in FCS >= 3.2");
                std::vector<char> buf(b_siz[i]);
                fi.read(buf.data(), b_siz[i]);
                std::string s(buf.begin(), buf.end());
                out[k] = std::stod(s.c_str());
              }
                break;
              case 1: { fi.read((char *)&n32,sizeof(n32)); out[k] = n32; }
                break;
              case 2: { fi.read((char *)&n64,sizeof(n64)); out[k] = n64; }
                break;
              case 3: { fi.read((char *)&i08,sizeof(i08)); out[k] = i08 & m[i]; }
                break;
              case 4: { fi.read((char *)&i16,sizeof(i16)); out[k] = i16 & m[i]; }
                break;
              case 5: { fi.read((char *)&i32,sizeof(i32)); out[k] = i32 & m[i]; }
                break;
              case 6: { fi.read((char *)&i64,sizeof(i64)); out[k] = i64 & m[i]; }
                break;
              default: { Rcpp::stop("hpp_readFCS: 'b_typ' and 'b_siz' not handled"); }
              break;
              }
            }
          }
        }
        return out;
      }
    }
    if(swap) {
      for(uint32_t o = 0; o < events; o++) {
        for(R_len_t i = 0; i < b_typ.size(); i++, k++) {
          switch(V[i]) {
          case 0: { 
            //Rcpp::stop("hpp_readFCS: type[A] is not supported in FCS >= 3.2");
            std::vector<char> buf(b_siz[i]);
            fi.read(buf.data(), b_siz[i]);
            std::string s(buf.begin(), buf.end());
            out[k] = std::stod(s.c_str());
          }
            break;
          case 1: { fi.read((char *)&n32,sizeof(n32)); out[k] = bytes_swap(n32); }
            break;
          case 2: { fi.read((char *)&n64,sizeof(n64)); out[k] = bytes_swap(n64); }
            break;
          case 3: { fi.read((char *)&i08,sizeof(i08)); out[k] = bytes_swap(i08); }
            break;
          case 4: { fi.read((char *)&i16,sizeof(i16)); out[k] = bytes_swap(i16); }
            break;
          case 5: { fi.read((char *)&i32,sizeof(i32)); out[k] = bytes_swap(i32); }
            break;
          case 6: { fi.read((char *)&i64,sizeof(i64)); out[k] = bytes_swap(i64); }
            break;
          default: { Rcpp::stop("hpp_readFCS: 'b_typ' and 'b_siz' not handled"); }
          break;
          }
        }
      } 
    } else {
      for(uint32_t o = 0; o < events; o++) {
        for(R_len_t i = 0; i < b_typ.size(); i++, k++) {
          switch(V[i]) {
          case 0: {
            //Rcpp::stop("hpp_readFCS: type[A] is not supported in FCS >= 3.2");
            std::vector<char> buf(b_siz[i]);
            fi.read(buf.data(), b_siz[i]);
            std::string s(buf.begin(), buf.end());
            out[k] = std::stod(s.c_str());
          }
            break;
          case 1: { fi.read((char *)&n32,sizeof(n32)); out[k] = n32; }
            break;
          case 2: { fi.read((char *)&n64,sizeof(n64)); out[k] = n64; }
            break;
          case 3: { fi.read((char *)&i08,sizeof(i08)); out[k] = i08; }
            break;
          case 4: { fi.read((char *)&i16,sizeof(i16)); out[k] = i16; }
            break;
          case 5: { fi.read((char *)&i32,sizeof(i32)); out[k] = i32; }
            break;
          case 6: { fi.read((char *)&i64,sizeof(i64)); out[k] = i64; }
            break;
          default: { Rcpp::stop("hpp_readFCS: 'b_typ' and 'b_siz' not handled"); }
          break;
          }
        }
      }
    }
    return out;
  } else {
    Rcpp::stop("hpp_readFCS: Unable to open file");
  }
  return 0;
}

#endif
