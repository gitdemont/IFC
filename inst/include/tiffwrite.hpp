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

/*
--------------------------------------------------------------------------------
|             functions described hereunder are experimental                   |
|              inputs and outputs may change in the future                     |
--------------------------------------------------------------------------------
*/

#ifndef IFC_TIFFWRITE_HPP
#define IFC_TIFFWRITE_HPP

#include <Rcpp.h>
#include "utils.hpp"
#include "import.hpp"
using namespace Rcpp;

static int tsizes[13] = {0,1,1,2,4,4,1,1,2,4,4,4,8};
static uint8_t wsizes[9] = {0,1,1,2,2,4,4,4,8};

// [[Rcpp::export(rng = false)]]
Rcpp::RawVector uint16_to_raw (const uint16_t x) {
  return Rcpp::RawVector::create((x      ) & 0xff, (x >>  8) & 0xff);
}

// [[Rcpp::export(rng = false)]]
Rcpp::RawVector uint32_to_raw (const uint32_t x) {
  return Rcpp::RawVector::create((x      ) & 0xff, (x >>  8) & 0xff, (x >> 16) & 0xff, (x >> 24) & 0xff);
}

// template to cast scalar value and copy it to vector
template <typename T> void toBYTE_T (T src,
                                     Rcpp::RawVector &dst,
                                     std::size_t at = 0,
                                     const bool swap = false) {
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(&src);
  if(swap) {
    std::reverse_copy(ptr, ptr + sizeof(src), dst.begin() + at);
  } else {
    std::copy(ptr, ptr + sizeof(src), dst.begin() + at);
  }
}

// template to cast vector to RawVector
template <typename T>
Rcpp::RawVector cast_vector_T (SEXP x,
                               const uint8_t what = 1,
                               const bool swap = false) {
  uint8_t siz = wsizes[what];
  switch( TYPEOF(x) ) {
    case RAWSXP: {
      return Rcpp::as<Rcpp::RawVector>(x);
    }
    case INTSXP: {
      std::vector<int> xx = as<std::vector<int>>(Rcpp::as<Rcpp::IntegerVector>(x));
      std::vector<T> foo(xx.begin(), xx.end());
      Rcpp::RawVector out = Rcpp::no_init_vector(foo.size() * siz);
      for(std::size_t i = 0, j = 0; i < foo.size(); i++, j += siz) toBYTE_T(foo[i], out, j, swap);
      return out;
    }
    case REALSXP: {
      std::vector<double> xx = as<std::vector<double>>(Rcpp::as<Rcpp::NumericVector>(x));
      std::vector<T> foo(xx.begin(), xx.end());
      Rcpp::RawVector out = Rcpp::no_init_vector(foo.size() * siz);
      for(std::size_t i = 0, j = 0; i < foo.size(); i++, j += siz) toBYTE_T(foo[i], out, j, swap);
      return out;
    }
    default: Rcpp::stop("cast_vector: 'x' not supported SEXPTYPE[%s]", Rcpp::type2name(x));
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
Rcpp::RawVector cast_vector ( SEXP x,
                              const uint8_t what = 1,
                              const bool swap = false) {
  switch( what ) {
  case 1: return cast_vector_T<int8_t>(x, what, swap);
  case 2: return cast_vector_T<uint8_t>(x, what, swap);
  case 3: return cast_vector_T<int16_t>(x, what, swap);
  case 4: return cast_vector_T<uint16_t>(x, what, swap);
  case 5: return cast_vector_T<int32_t>(x, what, swap);
  case 6: return cast_vector_T<uint32_t>(x, what, swap);
  case 7: return cast_vector_T<float_t>(x, what, swap);
  case 8: return cast_vector_T<double_t>(x, what, swap);
  }
  Rcpp::stop("cast_vector: bad 'what' value[%u]", what);
  return 0;
}

// template to cast image to RawVector
template <typename T>
Rcpp::RawVector cast_image_T (SEXP x,
                              const uint8_t what = 1,
                              const bool swap = false) {
  uint8_t siz = wsizes[what];
  Rcpp::RawVector def(1);
  def.attr("dims") = Rcpp::IntegerVector::create(1, 1, 1, 1, 1);
  def.attr("what") = Rcpp::IntegerVector::create(1);
  switch( TYPEOF(x) ) {
  case RAWSXP: {
    Rcpp::RawVector img = Rcpp::as<Rcpp::RawVector>(x);
    Rcpp::IntegerVector d = img.hasAttribute("dim") ? img.attr("dim") : Rcpp::IntegerVector::create(img.size(), 1, 1, 1);
    while(d.size() < 4) d.push_back(1);
    if(img.size() == 0) return def;
    Rcpp::RawVector out = Rcpp::no_init_vector(img.size());
    while(d.size() < 3) d.push_back(1);
    for(R_len_t h = 0, j = 0; h < d[0]; h++) {
      for(R_len_t w = 0; w < d[1]; w++) {
        for(R_len_t c = 0; c < d[2]; c++) {
          for(R_len_t f = 0; f < d[3]; f++) {
            out[j++] = img[h + w * d[0] + c * d[1] * d[0] + f * d[2] * d[1] * d[0]];
          }
        }
      }
    }
    out.attr("dims") = Rcpp::IntegerVector::create(siz, d[3], d[2], d[1], d[0]);
    out.attr("what") = Rcpp::IntegerVector::create(what);
    return out;
  }
  case INTSXP: {
    Rcpp::IntegerVector img = Rcpp::as<Rcpp::IntegerVector>(x);
    if(img.size() == 0) return def;
    Rcpp::IntegerVector d = img.hasAttribute("dim") ? img.attr("dim") : Rcpp::IntegerVector::create(img.size(), 1, 1, 1);
    while(d.size() < 4) d.push_back(1);
    Rcpp::RawVector out = Rcpp::no_init_vector(img.size() * siz);
    std::vector<int> xx = as<std::vector<int>>(img);
    std::vector<T> foo(xx.begin(), xx.end());
    for(R_len_t h = 0, j = 0; h < d[0]; h++) {
      for(R_len_t w = 0; w < d[1]; w++) {
        for(R_len_t c = 0; c < d[2]; c++) {
          for(R_len_t f = 0; f < d[3]; f++, j += siz) {
            toBYTE_T(foo[h + w * d[0] + c * d[1] * d[0] + f * d[2] * d[1] * d[0]], out, j, swap);
          }
        }
      }
    }
    out.attr("dims") = Rcpp::IntegerVector::create(siz, d[3], d[2], d[1], d[0]);
    out.attr("what") = Rcpp::IntegerVector::create(what);
    return out;
  }
  case REALSXP: {
    Rcpp::NumericVector img = Rcpp::as<Rcpp::NumericVector>(x);
    if(img.size() == 0) return def;
    Rcpp::IntegerVector d = img.hasAttribute("dim") ? img.attr("dim") : Rcpp::IntegerVector::create(img.size(), 1, 1, 1);
    while(d.size() < 4) d.push_back(1);
    Rcpp::RawVector out = Rcpp::no_init_vector(img.size() * siz);
    std::vector<double> xx = as<std::vector<double>>(img);
    std::vector<T> foo(xx.begin(), xx.end());
    for(R_len_t h = 0, j = 0; h < d[0]; h++) {
      for(R_len_t w = 0; w < d[1]; w++) {
        for(R_len_t c = 0; c < d[2]; c++) {
          for(R_len_t f = 0; f < d[3]; f++, j += siz) {
            toBYTE_T(foo[h + w * d[0] + c * d[1] * d[0] + f * d[2] * d[1] * d[0]], out, j, swap);
          }
        }
      }
    }
    out.attr("dims") = Rcpp::IntegerVector::create(siz, d[3], d[2], d[1], d[0]);
    out.attr("what") = Rcpp::IntegerVector::create(what);
    return out;
  }
  }
  return def;
}

//' @title Cast Image
//' @name cpp_cast_image
//' @description
//' Casts image from RAW, INT or REAL SEXP vector
//' @param x SEXP, the image to cast
//' @param what uint8_t, type to use for casting. Default is \code{1}.
//' Allowed are 1=int8_t, 2=uint8_t, 3=int16_t, 4=uint16_t, 5=int32_t, 6=uint32_t, 7=float_t, 8=double_t.
//' @param swap bool, whether single scalar values of \code{x} should be swap. Default is \code{false}.
//' @return a Rcpp::RawVector
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_cast_image ( SEXP x,
                                 const uint8_t what = 1,
                                 const bool swap = false) {
  switch( what ) {
  case 1: return cast_image_T<int8_t>(x, what, swap);
  case 2: return cast_image_T<uint8_t>(x, what, swap);
  case 3: return cast_image_T<int16_t>(x, what, swap);
  case 4: return cast_image_T<uint16_t>(x, what, swap);
  case 5: return cast_image_T<int32_t>(x, what, swap);
  case 6: return cast_image_T<uint32_t>(x, what, swap);
  case 7: return cast_image_T<float_t>(x, what, swap);
  case 8: return cast_image_T<double_t>(x, what, swap);
  }
  Rcpp::stop("hpp_cast_image: bad 'what' value[%u]", what);
  return 0;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List sub_list ( const Rcpp::List L,
                      const Rcpp::Nullable<Rcpp::IntegerVector> x_ = R_NilValue) {
  if(x_.isNotNull()) {
    Rcpp::IntegerVector x(x_.get());
    std::vector<int> xx;
    for(R_len_t i = 0; i < x.size(); i++) if((x[i] != NA_INTEGER) && (x[i] >= 0) && (x[i] < L.size())) xx.push_back(x[i]);
    if(xx.size() > 0) {
      Rcpp::IntegerVector xxx(xx.begin(), xx.end());
      return L[xxx];
    }
  }
  return L;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List hpp_tags_clean ( const Rcpp::List IFD, const R_len_t pos = 0, const bool verbose = false) {
  std::vector<double> tokeep;
  std::vector<uint16_t> tags;
  if(verbose) {
    for(R_len_t i = 0; i < IFD.size(); i++) {
      Rcpp::List ifd = IFD[i];
      if(!(ifd.containsElementNamed("tag") && ifd.containsElementNamed("typ") && ifd.containsElementNamed("map"))) {
        Rprintf("'ifd'@%f should contain 'tag', 'typ' and 'map' elements\n", i - pos);
        continue; 
      }
      bool isok = true;
      int32_t itag = ifd["tag"];
      if((itag <= 0) || (itag > 65535)) {
        Rprintf("'ifd'@%u 'tag'[=%i] should be [1-65535]\n", i - pos, itag);
        isok = false;
      }
      int32_t ityp = ifd["typ"];
      if((ityp <= 0) || (ityp > 12)) {
        Rprintf("'ifd'@%u 'tag'[=%i] typ'[%i] should be [1-12]\n", i - pos, itag, ityp);
        isok = false;
      }
      if(SEXPsize(ifd["map"]) == 0) {
        Rprintf("'ifd'@%u 'tag'[=%i] 'map' should not be a 0-length vector\n", i - pos, itag);
        isok = false;
      }
      if(((ityp == 5) || (ityp == 10)) && (SEXPsize(ifd["map"]) % 2)) {
        Rprintf("'ifd'@%u 'tag'[=%i] 'map' should be even for 'typ'[=%i]\n", i - pos, itag, ityp);
        isok = false;
      }
      if(isok && (std::find(tags.begin(), tags.end(), itag) == tags.end())) {
        tokeep.push_back(i);
        tags.push_back(itag); 
      }
    }
  } else {
    for(int i = 0; i < IFD.size(); i++) {
      Rcpp::List ifd = IFD[i];
      if(!(ifd.containsElementNamed("tag") && ifd.containsElementNamed("typ") && ifd.containsElementNamed("map"))) continue; 
      bool isok = true;
      int32_t itag = ifd["tag"];
      if((itag <= 0) || (itag > 65535)) isok = false;
      int32_t ityp = ifd["typ"];
      if((ityp <= 0) || (ityp > 12)) isok = false; 
      if(SEXPsize(ifd["map"]) == 0) isok = false;
      if(((ityp == 5) || (ityp == 10)) && (SEXPsize(ifd["map"]) % 2)) isok = false;
      if(isok && (std::find(tags.begin(), tags.end(), itag) == tags.end())) {
        tokeep.push_back(i);
        tags.push_back(itag); 
      }
    } 
  }
  if(tokeep.size() > 0) {
    if(tokeep.size() > 65535) Rcpp::stop("hpp_tags_clean: IFD can't have more than 65535 entries");
    Rcpp::List LL = sub_list(IFD, Rcpp::wrap(tokeep));
    Rcpp::IntegerVector idx = seq(0, LL.size() - 1);
    std::sort(idx.begin(), idx.end(), [&](int i, int j) { return tags[i] < tags[j]; });
    return LL[idx];
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
std::size_t hpp_strsxp_bytes (const Rcpp::Nullable<Rcpp::CharacterVector> x_ = R_NilValue) {
  std::size_t count = 0;
  if(x_.isNotNull()) {
    Rcpp::CharacterVector x(x_.get());
    for(R_len_t i = 0; i < x.size(); i++) count += (Rcpp::as<std::string>(x[i])).size();
    if(x.size()) count++;
  }
  return count;
}

// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_string_to_raw (const Rcpp::Nullable<Rcpp::CharacterVector> x_ = R_NilValue,
                                   const R_len_t extract_max = 4) {
  R_len_t count = 0;
  Rcpp::RawVector out(0);
  if(x_.isNotNull()) {
    Rcpp::CharacterVector x(x_.get());
    R_len_t m = out.size() < extract_max ? extract_max : out.size();
    for(R_len_t i = 0; i < x.size(); i++) {
      std::string y = Rcpp::as<std::string>(x[i]);
      R_len_t k = y.size();
      count += k;
      R_len_t n = k < m ? k : m;
      Rcpp::RawVector foo = Rcpp::no_init(n);
      std::copy(y.begin(), y.begin() + n, foo.begin());
      out = hpp_c(out, foo);
      m = m - n;
      if(m <= 0) break;
    }
  }
  return hpp_c(out, Rcpp::RawVector(1));
}

// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_pad_raw (const Rcpp::Nullable<Rcpp::RawVector> x_ = R_NilValue,
                             const bool front = false,
                             const R_len_t final = 4) {
  if(x_.isNotNull()) {
    Rcpp::RawVector x(x_.get());
    R_len_t s = x.size();
    if(s >= final) return x;
    Rcpp::RawVector foo(final - s);
    if(front) return hpp_c(foo, x);
    return hpp_c(x, foo);
  }
  return final;
}

// template to create tag extra content
template <int RTYPE>
Rcpp::RawVector hpp_tag_extcnt_T ( const Rcpp::Vector<RTYPE>&map,
                                   const uint16_t typ,
                                   const uint16_t tag,
                                   const bool swap = false) {
  if((typ > 12) || (typ < 1)) Rcpp::stop("hpp_tag_extcnt: 'typ' value[%u] not allowed", typ);
  if(tsizes[typ] * map.size() > 4294967295) Rcpp::stop("hpp_tag_extcnt: 'map' is too big");
  uint32_t bytes = tsizes[typ] * map.size();
  if(typ == 2) bytes = hpp_strsxp_bytes(as<Rcpp::CharacterVector>(map));
  if(bytes <= 0) return 0;
  switch(typ) {
  case 1: return as<Rcpp::RawVector>(map);
  case 2: return hpp_string_to_raw(as<Rcpp::CharacterVector>(map), bytes);
  case 3: return cast_vector(map, 4, swap);
  case 4: return cast_vector(map, 6, swap);
  case 5: return cast_vector(map, 6, swap);
  case 6: return as<Rcpp::RawVector>(map);
  case 7: return as<Rcpp::RawVector>(map);
  case 8: return cast_vector(map, 3, swap);
  case 9: return cast_vector(map, 5, swap);
  case 10: return cast_vector(map, 5, swap);
  case 11: return cast_vector(map, 7, swap);
  case 12: return cast_vector(map, 8, swap);
  }
  return 0;
}

// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_tag_extcnt ( SEXP map,
                                 const uint16_t typ,
                                 const uint16_t tag,
                                 const bool swap = false) {
  switch( TYPEOF(map) ) {
  case NILSXP: return hpp_tag_extcnt_T<RAWSXP>(Rcpp::RawVector(0), typ, tag, swap);
  case INTSXP: return hpp_tag_extcnt_T<INTSXP>(map, typ, tag, swap);
  case REALSXP: return hpp_tag_extcnt_T<REALSXP>(map, typ, tag, swap);
  case STRSXP: return hpp_tag_extcnt_T<STRSXP>(map, typ, tag, swap);
  case RAWSXP: return hpp_tag_extcnt_T<RAWSXP>(map, typ, tag, swap);
  default: Rcpp::stop("hpp_tag_extcnt: 'map' not supported SEXPTYPE[%s]", Rcpp::type2name(map));
  }
}

// template to create tag min content
template <int RTYPE>
Rcpp::RawVector hpp_tag_mincnt_T ( const Rcpp::Vector<RTYPE>&map,
                                   const uint16_t typ,
                                   const uint16_t tag,
                                   const uint32_t pos = 0,
                                   const bool swap = false) {
  if((typ > 12) || (typ < 1)) Rcpp::stop("hpp_tag_mincnt: 'typ' value[%u] not allowed", typ);
  if(tsizes[typ] * map.size() > 4294967295) Rcpp::stop("hpp_tag_mincnt: 'map' is too big");
  R_len_t bytes = tsizes[typ] * map.size();
  if(typ == 2) bytes = hpp_strsxp_bytes(as<Rcpp::CharacterVector>(map));
  bool multi = (typ == 5) || (typ == 10);
  Rcpp::RawVector out = Rcpp::no_init(12),
    typ_r = uint16_to_raw(typ),
    tag_r = uint16_to_raw(tag),
    cnt_r = uint32_to_raw(bytes / tsizes[typ] / (multi + 1)),
    val_r;
  if(bytes > 4) {
    val_r = uint32_to_raw(pos);
  } else {
    val_r = hpp_pad_raw(hpp_tag_extcnt(map, typ, tag, false), swap && (tsizes[typ] != 1));
  }
  if(swap) {
    out[ 0] = tag_r[1];
    out[ 1] = tag_r[0];
    out[ 2] = typ_r[1];
    out[ 3] = typ_r[0];
    out[ 4] = cnt_r[3];
    out[ 5] = cnt_r[2];
    out[ 6] = cnt_r[1];
    out[ 7] = cnt_r[0];
    if((bytes <= 4) && (tsizes[typ] == 1)) {
      out[ 8] = val_r[0];
      out[ 9] = val_r[1];
      out[10] = val_r[2];
      out[11] = val_r[3]; 
    } else {
      out[ 8] = val_r[3];
      out[ 9] = val_r[2];
      out[10] = val_r[1];
      out[11] = val_r[0];
    }
  } else {
    out[ 0] = tag_r[0];
    out[ 1] = tag_r[1];
    out[ 2] = typ_r[0];
    out[ 3] = typ_r[1];
    out[ 4] = cnt_r[0];
    out[ 5] = cnt_r[1];
    out[ 6] = cnt_r[2];
    out[ 7] = cnt_r[3];
    out[ 8] = val_r[0];
    out[ 9] = val_r[1];
    out[10] = val_r[2];
    out[11] = val_r[3]; 
  }
  out.attr("bytes") = bytes;
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_tag_mincnt ( SEXP map,
                                 const uint16_t typ,
                                 const uint16_t tag,
                                 const uint32_t pos = 0,
                                 const bool swap = false) {
  switch( TYPEOF(map) ) {
  case NILSXP: return hpp_tag_mincnt_T<RAWSXP>(Rcpp::RawVector(0), typ, tag, pos, swap);
  case INTSXP: return hpp_tag_mincnt_T<INTSXP>(map, typ, tag, pos, swap);
  case REALSXP: return hpp_tag_mincnt_T<REALSXP>(map, typ, tag, pos, swap);
  case STRSXP: return hpp_tag_mincnt_T<STRSXP>(map, typ, tag, pos, swap);
  case RAWSXP: return hpp_tag_mincnt_T<RAWSXP>(map, typ, tag, pos, swap);
  default: Rcpp::stop("hpp_tag_mincnt: 'map' not supported SEXPTYPE[%s]", Rcpp::type2name(map));
  }
}

//' @title IFD Tag Writer
//' @name cpp_writeIFD
//' @description
//' Writes TIFF IFD (Image Field Directory).
//' @param img RawVector, an encoded image. It should contain 'dims', 'what' and 'comp' attributes.
//' @param tags List, extra tags to be included to IFD. Expecting a list whose sub-elements are list containing:\cr
//' -'tag' uint16_t, tag number, it should be [1-65535],\cr
//' -'typ' uint16_t typ number, it should be [1-12],\cr
//' -'map' SEXP vector of values to write, it should not be empty and should be even for 'typ' 5 and 10.
//' @param offset uint32_t, position of the IFD beginning. Default is \code{0}.
//' @param endianness std::string, "little" or "big".
//' @param rgb bool, whether to write channels as rgb. Default is \code{false}.
//' @param last bool, whether IFD is the last one. Default is \code{false}.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is \code{false}.
//' @return a RawVector
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_writeIFD ( const Rcpp::RawVector img,
                               const Rcpp::List tags,
                               const uint32_t offset = 0,
                               const std::string endianness = "little",
                               const bool rgb = false,
                               const bool last = false,
                               const bool verbose = false) {
  if(!(img.hasAttribute("dims") && img.hasAttribute("what") && img.hasAttribute("comp"))) Rcpp::stop("hpp_writeIFD: bad 'img' format");
  Rcpp::IntegerVector d = img.attr("dims");
  if(d.size() != 5) Rcpp::stop("hpp_writeIFD: 'img' illegal attr(img, \"dims\")");
  uint8_t what = img.attr("what");
  if(what <= 0 || what >= 9) Rcpp::stop("hpp_writeIFD: 'img' illegal attr(img, \"what\")");
  uint16_t comp = img.attr("comp");
  bool swap = endianness != hpp_getEndian();
  Rcpp::CharacterVector ver = hpp_getversion("IFC");
  if(ver.size()) ver = hpp_c(Rcpp::CharacterVector::create(" "), ver);
  ver = hpp_c(Rcpp::CharacterVector::create("IFC"), ver);
  
  Rcpp::List default_ifd = Rcpp::List::create(
    Rcpp::List::create(_["tag"] = 256, _["typ"] = d[2] > 65535 ? 4 : 3, _["map"] = d[3]),
    Rcpp::List::create(_["tag"] = 257, _["typ"] = d[3] > 65535 ? 4 : 3, _["map"] = d[4]),
    Rcpp::List::create(_["tag"] = 258, _["typ"] = 3, _["map"] = rep(8 * d[0], d[2])),
    Rcpp::List::create(_["tag"] = 259, _["typ"] = 3, _["map"] = comp),
    Rcpp::List::create(_["tag"] = 262, _["typ"] = 3, _["map"] = (rgb && (d[2] == 3)) ? 2 : 1),
    Rcpp::List::create(_["tag"] = 273, _["typ"] = 4, _["map"] = 0),
    Rcpp::List::create(_["tag"] = 277, _["typ"] = 3, _["map"] = (rgb && (d[2] == 3)) ? 3 : d[2]),
    Rcpp::List::create(_["tag"] = 278, _["typ"] = d[4] > 65535 ? 4 : 3, _["map"] = d[4]),
    Rcpp::List::create(_["tag"] = 279, _["typ"] = 4, _["map"] = 0),
    Rcpp::List::create(_["tag"] = 284, _["typ"] = 3, _["map"] = 1),
    Rcpp::List::create(_["tag"] = 305, _["typ"] = 2, _["map"] = Rcpp::collapse(ver)),
    Rcpp::List::create(_["tag"] = 339, _["typ"] = 3, _["map"] = rep(what > 6 ? 3 : what % 2 ? 2 : 1, d[2]))
  );
  if(!((rgb && (d[2] == 3)) || (d[2] <= 1))) {
    default_ifd = hpp_c(default_ifd, Rcpp::List::create(Rcpp::List::create(_["tag"] = 338, _["typ"] = 3, _["map"] = rep(0, d[2] - 1)))); 
  }
  Rcpp::List cleaned_ifd = hpp_tags_clean(hpp_c(default_ifd, tags), default_ifd.size() - 1, verbose);
  
  Rcpp::RawVector off(4);
  Rcpp::RawVector ext_content(0);
  Rcpp::RawVector min_content = Rcpp::no_init_vector(12 * cleaned_ifd.size());
  std::size_t npos = (offset + 12 * cleaned_ifd.size() + 2 + 4) % 4294967296;
  
  R_len_t k = -1;
  for(R_len_t i = 0; i < cleaned_ifd.size(); i++) {
    Rcpp::List ifd = cleaned_ifd[i];
    uint16_t tag = ifd["tag"];
    if(tag == 273) k = i;
    if(tag == 279) {
      off = uint32_to_raw(npos);
      npos = (npos + img.size()) % 4294967296;
      ext_content = hpp_c(ext_content, img);
      Rcpp::RawVector v = hpp_tag_mincnt(Rcpp::NumericVector::create(img.size()), 4, 279, npos, swap);
      std::copy(v.begin(), v.end(), min_content.begin() + i * 12);
      continue;
    }
    Rcpp::RawVector v = hpp_tag_mincnt(ifd["map"], ifd["typ"], ifd["tag"], npos, swap);
    uint32_t bytes = v.attr("bytes");
    std::copy(v.begin(), v.end(), min_content.begin() + i * 12);
    if(bytes > 4) {
      npos = (npos + bytes) % 4294967296;
      ext_content = hpp_c(ext_content, hpp_tag_extcnt(ifd["map"], ifd["typ"], ifd["tag"], swap));
    }
  }
  if(last) npos = 0;
  Rcpp::RawVector EI = uint16_to_raw(cleaned_ifd.size());
  Rcpp::RawVector PI = uint32_to_raw(npos);
  if(swap) {
    std::reverse(off.begin(), off.end());
    std::reverse(EI.begin(), EI.end());
    std::reverse(PI.begin(), PI.end());
  }
  if(k >= 0) std::copy(off.begin(), off.end(), min_content.begin() + k * 12 + 8);
  Rcpp::RawVector out = hpp_c(hpp_c(EI, min_content), hpp_c(PI, ext_content));
  out.attr("offset") = npos;
  return out;
}

#endif
