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

static uint64_t pow_v[8]  = {1,256,65536,16777216,4294967296,1099511627776,281474976710656,72057594037927936};

uint64_t bytes_conv(const std::vector<char> v) {
  uint64_t x = 0;
  for(std::size_t i = 0; i < std::min(v.size(), sizeof(x)); i++) x += (v[i] & 0xff) * pow_v[i];
  return x;
}

uint64_t bytes_swap_conv(const std::vector<char> v) {
  uint64_t x = 0;
  for(std::size_t k = 0, i = std::min(v.size(), sizeof(x)); i > 0; k++) x += (v[k] & 0xff) * pow_v[--i];
  return x;
}

//' @title FCS data extraction
//' @name cpp_readFCS
//' @description
//' This function reads data from FCS 3.2 files
//' @param fname string, path to file.
//' @param offset std::size_t, offset position of data start
//' @param events uint32_t, number of events to read.
//' @param b_typ Rcpp::IntegerVector, types of values to read. Allowed are 0 for "A", 1 for "F", 2 for "D", and 3 is "I".
//' @param b_siz Rcpp::IntegerVector, number of bytes to extract for each type. Values should be higher than 0.\cr
//' Note that whatever the input is, when 'b_typ' is 1 (float) 'b_siz' will be set to 4 and when 'b_typ' is 2 (double) 'b_siz' will be set to 8.\cr
//' Note also that when 'b_typ' is 3 (integer), the only allowed values for 'b_siz' are [1-8].\cr
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
  if(b_typ.size() != b_siz.size()) {
    Rcpp::Rcerr << "hpp_readFCS: 'b_typ'[" << b_typ.size() << "] and 'b_siz'["<< b_siz.size() <<"] should be of same length" << std::endl;
    Rcpp::stop("hpp_readFCS: 'b_typ' and 'b_siz' should be of same length");
  }
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    fi.seekg(0, std::ios::end);
    std::size_t filesize = fi.tellg();
    
    std::size_t count = 0;
    for(R_len_t i = 0; i < b_typ.size(); i++) {
      count += b_siz[i];
      switch(b_typ[i]) {
      case 0 :
        if(b_siz[i] <= 0) {
          Rcpp::Rcerr << "hpp_readFCS: size[" << b_siz[i] << "] not handled for b_typ[A]" << std::endl;
          Rcpp::stop("hpp_readFCS: invalid value in 'b_siz'");
        }
        break;
      case 1 :
        if(b_siz[i] != 4) {
          Rcpp::Rcerr << "hpp_readFCS: size[" << b_siz[i] << "] not handled for b_typ[F]" << std::endl;
          Rcpp::stop("hpp_readFCS: invalid value in 'b_siz'");
        }
        break;
      case 2 :
        if(b_siz[i] != 8) {
          Rcpp::Rcerr << "hpp_readFCS: size[" << b_siz[i] << "] not handled for b_typ[D]" << std::endl;
          Rcpp::stop("hpp_readFCS: invalid value in 'b_siz'");
        } 
        break;
      case 3 :
        if((b_siz[i] < 1) || (b_siz[i] > 8)) {
          Rcpp::Rcerr << "hpp_readFCS: size[" << b_siz[i] << "] not handled for b_typ[I]" << std::endl;
          Rcpp::stop("hpp_readFCS: invalid value in 'b_siz'");
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
    float    n32;
    double   n64;
    R_len_t k = 0;
    
    fi.seekg(offset, std::ios::beg);
    if(b_msk.isNotNull()) {
      Rcpp::IntegerVector masks(b_msk.get());
      if(masks.size() != 0) {
        if(b_siz.size() != masks.size()) {
          Rcpp::Rcerr << "hpp_readFCS: when provided 'b_msk'[" << masks.size() << "] and 'b_siz'["<< b_siz.size() <<"] should be of same length" << std::endl;
          Rcpp::stop("hpp_readFCS: 'b_msk' and 'b_siz' should be of same length");
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
              switch(b_typ[i]) {
              case 0: {
                //Rcpp::stop("hpp_readFCS: type[A] is not supported in FCS >= 3.2");
                std::vector<char> buf_str(b_siz[i]);
                fi.read(buf_str.data(), b_siz[i]);
                std::string s(buf_str.begin(), buf_str.end());
                out[k] = std::stod(s.c_str());
              }
                break;
              case 1: { fi.read((char *)&n32,sizeof(n32)); out[k] = bytes_swap(n32); }
                break;
              case 2: { fi.read((char *)&n64,sizeof(n64)); out[k] = bytes_swap(n64); }
                break;
              case 3: { 
                std::vector<char>buf_int(b_siz[i]);
                fi.read(buf_int.data(),b_siz[i]); out[k] = bytes_swap_conv(buf_int) & m[i];
              }
                break;
              }
            }
          } 
        } else {
          for(uint32_t o = 0; o < events; o++) {
            for(R_len_t i = 0; i < b_typ.size(); i++, k++) {
              switch(b_typ[i]) {
              case 0: {
                //Rcpp::stop("hpp_readFCS: type[A] is not supported in FCS >= 3.2");
                std::vector<char> buf_str(b_siz[i]);
                fi.read(buf_str.data(), b_siz[i]);
                std::string s(buf_str.begin(), buf_str.end());
                out[k] = std::stod(s.c_str());
              }
                break;
              case 1: { fi.read((char *)&n32,sizeof(n32)); out[k] = n32; }
                break;
              case 2: { fi.read((char *)&n64,sizeof(n64)); out[k] = n64; }
                break;
              case 3: {
                std::vector<char>buf_int(b_siz[i]);
                fi.read(buf_int.data(),b_siz[i]); out[k] = bytes_conv(buf_int) & m[i];
              }
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
          switch(b_typ[i]) {
          case 0: {
            //Rcpp::stop("hpp_readFCS: type[A] is not supported in FCS >= 3.2");
            std::vector<char> buf_str(b_siz[i]);
            fi.read(buf_str.data(), b_siz[i]);
            std::string s(buf_str.begin(), buf_str.end());
            out[k] = std::stod(s.c_str());
          }
            break;
          case 1: { fi.read((char *)&n32,sizeof(n32)); out[k] = bytes_swap(n32); }
            break;
          case 2: { fi.read((char *)&n64,sizeof(n64)); out[k] = bytes_swap(n64); }
            break;
          case 3: {
            std::vector<char>buf_int(b_siz[i]);
            fi.read(buf_int.data(),b_siz[i]); out[k] = bytes_swap_conv(buf_int);
          }
            break;
          }
        }
      } 
    } else {
      for(uint32_t o = 0; o < events; o++) {
        for(R_len_t i = 0; i < b_typ.size(); i++, k++) {
          switch(b_typ[i]) {
          case 0: {
            //Rcpp::stop("hpp_readFCS: type[A] is not supported in FCS >= 3.2");
            std::vector<char> buf_str(b_siz[i]);
            fi.read(buf_str.data(), b_siz[i]);
            std::string s(buf_str.begin(), buf_str.end());
            out[k] = std::stod(s.c_str());
          }
            break;
          case 1: { fi.read((char *)&n32,sizeof(n32)); out[k] = n32; }
            break;
          case 2: { fi.read((char *)&n64,sizeof(n64)); out[k] = n64; }
            break;
          case 3: {
            std::vector<char>buf_int(b_siz[i]);
            fi.read(buf_int.data(),b_siz[i]); out[k] = bytes_conv(buf_int);
          }
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
