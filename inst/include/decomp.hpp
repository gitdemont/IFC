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

#ifndef IFC_DECOMP_HPP
#define IFC_DECOMP_HPP

#include <Rcpp.h>
#include "utils.hpp"
#include "trans.hpp"
using namespace Rcpp;

//' @title Raw to 32 bits Integer Conversion
//' @name cpp_raw_to_int32
//' @description
//' Converts raw vector to 32 bits integer vector
//' @param x RawVector.
//' @param bits uint8_t, bits depth. Default is 16.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector hpp_raw_to_int32 (const Rcpp::RawVector x,
                                      const uint8_t bits = 16) {
  if(!((bits == 8) || (bits == 16) || (bits == 24) || (bits == 32))) Rcpp::stop("hpp_raw_to_int32: 'bits' should be either 8, 16, 24 or 32");
  int b = bits / 8;
  if(x.size() % b) Rcpp::stop("hpp_raw_to_int32: 'x' size is not a multiple of 'bits'");
  Rcpp::IntegerVector out = Rcpp::no_init_vector(x.size() / b);
  uint8_t bb = b - 1;
  for(R_len_t i = 0, k = 0; i < out.size(); i++) {
    uint32_t r = 0;
    for(uint8_t j = 0; j <= bb; j++) r += x[k++] << (8 * j);
    out[i] = r;
  }
  out.attr("bits") = bits;
  return out;
}

//' @title Force Integer Sign
//' @name cpp_sign_int
//' @description
//' Force <32 bits integers to be signed
//' @param x IntegerVector.
//' @param bits uint8_t, bits depth. Default is 16.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector hpp_sign_int (const Rcpp::IntegerVector x,
                                  const uint8_t bits = 16) {
  if(!((bits == 8) || (bits == 16) || (bits == 24) || (bits == 32))) Rcpp::stop("hpp_sign_int: 'bits' should be either 8, 16, 24 or 32");
  if(bits == 32) return x;
  Rcpp::IntegerVector out = Rcpp::no_init_vector(x.size());
  uint32_t MAXR = std::pow(2.0, bits - 1) - 1;
  uint32_t MAXV = std::pow(2.0, bits);
  for(R_len_t i = 0; i < out.size(); i++) out[i] = hpp_int32_to_uint32(x[i]) > MAXR ? hpp_int32_to_uint32(x[i]) - MAXV : x[i];
  if(x.hasAttribute("dim")) out.attr("dim") = x.attr("dim");
  return out;
}

//' @title NONE Decompression
//' @name cpp_none_Decomp
//' @description
//' Operates none decompression of compressed image.
//' @param raw_chnk, a RawVector of compressed image.
//' @param imgWidth uint32_t, Width of the decompressed image. Default is 1.
//' @param imgHeight uint32_t, Height of the decompressed image. Default is 1.
//' @param nb_channels uint32_t, number of channels of the decompressed image. Default is 1.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List hpp_none_Decomp (const Rcpp::RawVector raw_chnk, 
                            const uint32_t imgWidth = 1, 
                            const uint32_t imgHeight = 1, 
                            const uint32_t nb_channels = 1,
                            const bool verbose = false) {
  R_len_t nbytes = raw_chnk.size();
  if((nb_channels * imgWidth * imgHeight * nbytes) != 0) {
    Rcpp::List out(nb_channels);
    uint32_t tile_width = imgWidth / nb_channels;
    Rcpp::IntegerVector V = hpp_raw_to_int32(raw_chnk, 8 * nbytes / (imgWidth * imgHeight));
    V.attr("dim") = Rcpp::Dimension(imgWidth, imgHeight);
    Rcpp::IntegerMatrix img = Rcpp::transpose(Rcpp::as<Rcpp::IntegerMatrix>(V));
    for(uint32_t i = 0; i < nb_channels; i++) {
      out[i] = img(Rcpp::_, Rcpp::Range(tile_width * i, tile_width * (i + 1) - 1));
      if(V.hasAttribute("bits")) out[i] = hpp_sign_int(out[i], V.attr("bits"));
    }
    if(V.hasAttribute("bits")) out.attr("bits") = V.attr("bits");
    return out;
  } else {
    Rcpp::stop("hpp_none_Decomp: raw_chnk, imgWidth, imgHeight and nb_channels should be >0");    
  }
  return R_NilValue;
}

//' @title RLE Decompression
//' @name cpp_rle_Decomp
//' @description
//' Operates RLE decompression of compressed image.
//' @param raw_chnk, a RawVector of compressed image.
//' @param imgWidth uint32_t, Width of the decompressed image. Default is 1.
//' @param imgHeight uint32_t, Height of the decompressed image. Default is 1.
//' @param nb_channels uint32_t, number of channels of the decompressed image. Default is 1.
//' @param removal uint8_t, object removal method. Default is no removal. Otherwise, if\cr
//' -1, for clipped removal: height OR width clipped pixels will be set to -1.\cr
//' -2, height clipped removal: height clipped pixels will be set to -1.\cr
//' -3, width clipped removal: width clipped pixels will be set to -1.\cr
//' -4, only keep background: background pixels will be set to 1 and all others to -1.\cr
//' -5, only keep foreground: foreground pixels will be set to 1 and all others to -1
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}\cr
//' /verb{
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'     this list of conditions and the following disclaimer in the documentation
//'     and/or other materials provided with the distribution.
//'  
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' }
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List hpp_rle_Decomp (const RawVector raw_chnk,
                           const uint32_t imgWidth = 1,
                           const uint32_t imgHeight = 1,
                           const uint32_t nb_channels = 1,
                           const uint8_t removal = 0,
                           const bool verbose = false) {
  R_len_t nbytes = raw_chnk.size();
  if((nb_channels * imgWidth * imgHeight * nbytes) != 0) {
    Rcpp::List out(nb_channels);
    uint32_t tile_width = imgWidth / nb_channels;
    Rcpp::IntegerMatrix img = Rcpp::no_init(imgWidth,imgHeight);
    R_len_t L = imgWidth * imgHeight, runLength = 0, j = 0;
    
    switch(removal) {
    case 1: { // clipped removal
      for(R_len_t k = 0; k < nbytes; k++) {
      int value = raw_chnk[k++];
      if(value > 1) value = -1;
      R_len_t off = runLength;
      runLength = off + (raw_chnk[k] & 0xff) + 1;
      if (runLength > L) Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
      for(j = off; j < runLength; j++) img[j] = value;
    }
      break;
    }
    case 2: { // height clipped removal
      for(R_len_t k = 0; k < nbytes; k++) {
      int value = raw_chnk[k++];
      if(value == 2) value = -1;
      R_len_t off = runLength;
      runLength = off + (raw_chnk[k] & 0xff) + 1;
      if (runLength > L) Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
      for(j = off; j < runLength; j++) img[j] = value;
    }
      break;
    }
    case 3: { // width clipped removal
      for(R_len_t k = 0; k < nbytes; k++) {
      int value = raw_chnk[k++];
      if(value == 3) value = -1;
      R_len_t off = runLength;
      runLength = off + (raw_chnk[k] & 0xff) + 1;
      if (runLength > L) Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
      for(j = off; j < runLength; j++) img[j] = value;
    }
      break;
    }
    case 4: { // only keep background
      for(R_len_t k = 0; k < nbytes; k++) {
      int value = raw_chnk[k++];
      value = (value == 0) ? 1:-1;
      R_len_t off = runLength;
      runLength = off + (raw_chnk[k] & 0xff) + 1;
      if (runLength > L) Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
      for(j = off; j < runLength; j++) img[j] = value;
    }
      break;
    }
    case 5: { // only keep non clipped foreground
      for(R_len_t k = 0; k < nbytes; k++) {
      int value = raw_chnk[k++];
      if(value != 1) value = -1;
      R_len_t off = runLength;
      runLength = off + (raw_chnk[k] & 0xff) + 1;
      if (runLength > L) Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
      for(j = off; j < runLength; j++) img[j] = value;
    }
      break;
    }
    default: { // no removal 
      for(R_len_t k = 0; k < nbytes; k++) {
      int value = raw_chnk[k++];
      R_len_t off = runLength;
      runLength = off + (raw_chnk[k] & 0xff) + 1;
      if (runLength > L) Rcpp::stop("hpp_rle_Decomp: Buffer overrun");
      for(j = off; j < runLength; j++) img[j] = value;
    }
      break;
    }
    }
    Rcpp::IntegerMatrix timg = Rcpp::transpose(img);
    for(uint32_t i = 0; i < nb_channels; i++) {
      out[i] = timg(Rcpp::_, Rcpp::Range(tile_width * i, tile_width * (i + 1) - 1));
    }
    return out;
  } else {
    Rcpp::stop("hpp_rle_Decomp: raw_chnk, imgWidth, imgHeight and nb_channels should be >0");    
  }
  return R_NilValue;
}

//' @title GRAY Decompression
//' @name cpp_gray_Decomp
//' @description
//' Operates GrayScale decompression of compressed image.
//' @param raw_chnk, a RawVector of compressed image.
//' @param imgWidth uint32_t, Width of the decompressed image. Default is 1.
//' @param imgHeight uint32_t, Height of the decompressed image. Default is 1.
//' @param nb_channels uint32_t, number of channels of the decompressed image. Default is 1.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}\cr
//' /verb{
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'    this list of conditions and the following disclaimer in the documentation
//'    and/or other materials provided with the distribution.
//' 
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' }
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List hpp_gray_Decomp (const Rcpp::RawVector raw_chnk, 
                            const uint32_t imgWidth = 1, 
                            const uint32_t imgHeight = 1, 
                            const uint32_t nb_channels = 1,
                            const bool verbose = false) {
  R_len_t nbytes = raw_chnk.size();
  if((nb_channels * imgWidth * imgHeight * nbytes) != 0) {
    Rcpp::List out(nb_channels);
    uint32_t tile_width = imgWidth / nb_channels;
    Rcpp::IntegerVector lastRow(imgWidth + 1);
    Rcpp::IntegerMatrix img = Rcpp::no_init(imgHeight, imgWidth + 1);
    for(uint32_t y = 0 ; y < imgHeight ; y++) img(y, 0) = 0;
    bool odd = false;
    
    R_len_t k = 0;
    for(uint32_t y = 0 ; y < imgHeight ; y++) {
      for(uint32_t x = 1 ; x <= imgWidth ; x++) {
        int value = 0;
        short shift = 0, nibble = -1;
        while((nibble & 0x8)) {
          if(odd) {
            nibble = raw_chnk[k++] >> 4;
          } else {
            if(k >= nbytes) Rcpp::stop("hpp_gray_Decomp: Buffer overrun");
            nibble = raw_chnk[k] & 0xf;
          }
          odd = !odd;
          value += (nibble & 0x7) << shift;
          shift += 3;
        }
        if(nibble & 0x4) value |= - (1 << shift);
        lastRow[x] += value;
        img(y,x) = img(y,x - 1) + lastRow[x];
      }
    }
    if(k != nbytes - odd) Rcpp::stop("hpp_gray_Decomp: Bad decompression");
    for(uint32_t i = 0; i < nb_channels; i++) {
      out[i] = img(Rcpp::_, Rcpp::Range(1 + tile_width * i, tile_width * (i + 1)));
    }
    return out;
  } else {
    Rcpp::stop("hpp_gray_Decomp: raw_chnk, imgWidth, imgHeight and nb_channels should be >0");    
  }
  return R_NilValue;
}

//' @title IFC_object Decompression
//' @name cpp_decomp
//' @description
//' Operates decompression of compressed image stored in TIFF file.
//' @param fname string, path to file.
//' @param offset std::size_t, position of the beginning of compressed image.
//' @param nbytes uint32_t, number of bytes of compressed image.
//' @param imgWidth uint32_t, Width of the decompressed image. Default is 1.
//' @param imgHeight uint32_t, Height of the decompressed image. Default is 1.
//' @param nb_channels uint32_t, number of channels of the decompressed image. Default is 1.
//' @param removal uint8_t, object removal method. Only apply for 30818 compression. Default is 0.\cr
//' -1, for clipped removal: height OR width clipped pixels will be set to -1.\cr
//' -2, height clipped removal: height clipped pixels will be set to -1.\cr
//' -3, width clipped removal: width clipped pixels will be set to -1.\cr
//' -4, only keep background: background pixels will be set to 1 and all others to 0.\cr
//' -5, only keep foreground: foreground pixels will be set to 1 and all others to 0.
//' @param compression uint32_t, compression algorithm used. Default is 30818.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}\cr
//' /verb{
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'    this list of conditions and the following disclaimer in the documentation
//'    and/or other materials provided with the distribution.
//' 
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' }
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List hpp_decomp (const std::string fname, 
                       const std::size_t offset, 
                       const uint32_t nbytes, 
                       const uint32_t imgWidth = 1, 
                       const uint32_t imgHeight = 1, 
                       const uint32_t nb_channels = 1,
                       const uint8_t removal = 0,
                       const uint32_t compression = 30818,
                       const bool verbose = false) {
  Rcpp::RawVector raw_chnk = hpp_readchunk(fname, offset, nbytes, verbose);
  switch(compression) {
  case 1: return hpp_none_Decomp(raw_chnk, imgWidth, imgHeight, nb_channels, verbose);
  case 30817: return hpp_gray_Decomp(raw_chnk, imgWidth, imgHeight, nb_channels, verbose);
  case 30818: return hpp_rle_Decomp(raw_chnk, imgWidth, imgHeight, nb_channels, removal, verbose);
  }
  Rcpp::Rcerr << "hpp_decomp: can't deal with compression format: " << compression << std::endl;
  Rcpp::stop("hpp_decomp: can't deal with compression format");   
  return R_NilValue;
}

//' @title GRAY Decompression to RAW
//' @name cpp_gray_rawDecomp
//' @description
//' Operates GrayScale decompression to raw of compressed image.
//' @param raw_chnk, a RawVector of compressed image..
//' @param imgWidth uint32_t, Width of the decompressed image. Default is 1.
//' @param imgHeight uint32_t, Height of the decompressed image. Default is 1.
//' @param swap bool, whether to swap bytes or not. Default is false.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}\cr
//' /verb{
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'    this list of conditions and the following disclaimer in the documentation
//'    and/or other materials provided with the distribution.
//' 
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' }
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_gray_rawDecomp (const Rcpp::RawVector raw_chnk,
                                    const uint32_t imgWidth = 1,
                                    const uint32_t imgHeight = 1,
                                    const bool swap = false,
                                    const bool verbose = false) {
  R_len_t nbytes = raw_chnk.size();
  if((imgWidth * imgHeight * nbytes) != 0) {
    Rcpp::RawVector out = Rcpp::no_init_vector(imgWidth * imgHeight * 2);
    Rcpp::IntegerVector lastRow(imgWidth + 1);
    Rcpp::IntegerMatrix img = Rcpp::no_init(imgHeight, imgWidth + 1);
    for(uint32_t y = 0 ; y < imgHeight ; y++) img(y, 0) = 0;
    bool odd = false;
    
    R_len_t k = 0;
    R_len_t i = 0;
    if(swap) {
      for(uint32_t y = 0; y < imgHeight ; y++) {
        for(uint32_t x = 1 ; x <= imgWidth ; x++) {
          int value = 0;
          short shift = 0, nibble = -1;
          while((nibble & 0x8)) {
            if(odd) {
              nibble = raw_chnk[k++] >> 4;
            } else {
              if(k >= nbytes) Rcpp::stop("hpp_gray_rawDecomp: Buffer overrun");
              nibble = raw_chnk[k] & 0xf;
            }
            odd = !odd;
            value += (nibble & 0x7) << shift;
            shift += 3;
          }
          if(nibble & 0x4) value |= - (1 << shift);
          lastRow[x] += value;
          img(y,x) = img(y,x - 1) + lastRow[x];
          if((i + 1) >= out.size()) Rcpp::stop("hpp_gray_rawDecomp: wrong size");
          // TODO: check swap
          out[i++] = (img(y,x) >> 24) & 0xff;
          out[i++] = (img(y,x) >> 16) & 0xff;
        }
      }
    } else {
      for(uint32_t y = 0; y < imgHeight ; y++) {
        for(uint32_t x = 1 ; x <= imgWidth ; x++) {
          int value = 0;
          short shift = 0, nibble = -1;
          while((nibble & 0x8)) {
            if(odd) {
              nibble = raw_chnk[k++] >> 4;
            } else {
              if(k >= nbytes) Rcpp::stop("hpp_gray_rawDecomp: Buffer overrun");
              nibble = raw_chnk[k] & 0xf;
            }
            odd = !odd;
            value += (nibble & 0x7) << shift;
            shift += 3;
          }
          if(nibble & 0x4) value |= - (1 << shift);
          lastRow[x] += value;
          img(y,x) = img(y,x - 1) + lastRow[x];
          if((i + 1) >= out.size()) Rcpp::stop("hpp_gray_rawDecomp: wrong size");
          out[i++] = (img(y,x)     ) & 0xff;
          out[i++] = (img(y,x) >> 8) & 0xff;
        }
      }
    }
    if(k != nbytes - odd) Rcpp::stop("hpp_gray_rawDecomp: Bad decompression");
    return out;
  } else {
    Rcpp::stop("hpp_gray_rawDecomp: raw_chnk, imgWidth and imgHeight should be >0");    
  }
  return R_NilValue;
}

//' @title RLE Decompression to RAW
//' @name cpp_rle_rawDecomp
//' @description
//' Operates RLE decompression to raw of compressed image.
//' @param raw_chnk, a RawVector of compressed image..
//' @param imgWidth uint32_t, Width of the decompressed image. Default is 1.
//' @param imgHeight uint32_t, Height of the decompressed image. Default is 1.
//' @param swap bool, whether to swap bytes or not. Default is false.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}\cr
//' /verb{
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'     this list of conditions and the following disclaimer in the documentation
//'     and/or other materials provided with the distribution.
//'  
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' }
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_rle_rawDecomp (const Rcpp::RawVector raw_chnk,
                                   const uint32_t imgWidth = 1,
                                   const uint32_t imgHeight = 1,
                                   const bool swap = false,
                                   const bool verbose = false) {
  R_len_t nbytes = raw_chnk.size();
  R_len_t L = imgWidth * imgHeight, runLength = 0, j = 0;
  if(L * nbytes != 0) {
    Rcpp::RawVector out = Rcpp::no_init_vector(L * 2);
    if(swap) {
      for(R_len_t k = 0; k < nbytes; k++) {
        int value = raw_chnk[k++];
        R_len_t off = runLength;
        runLength = off + (raw_chnk[k] & 0xff) + 1;
        if(runLength > L) Rcpp::stop("hpp_rle_rawDecomp: Buffer overrun");
        for(j = off; j < runLength; j++) {
          if((j*2 + 1) >= out.size()) Rcpp::stop("hpp_rle_rawDecomp: wrong size");
          out[j*2] = 0x00;
          out[j*2 + 1] = value & 0xff;
        }
      }
    } else {
      for(R_len_t k = 0; k < nbytes; k++) {
        int value = raw_chnk[k++];
        R_len_t off = runLength;
        runLength = off + (raw_chnk[k] & 0xff) + 1;
        if(runLength > L) Rcpp::stop("hpp_rle_rawDecomp: Buffer overrun");
        for(j = off; j < runLength; j++) {
          if((j*2 + 1) >= out.size()) Rcpp::stop("hpp_rle_rawDecomp: wrong size");
          out[j*2] = value & 0xff;
          out[j*2 + 1] = 0x00;
        }
      }
    }
    return out;
  } else {
    Rcpp::stop("hpp_rle_rawDecomp: raw_chnk, imgWidth and imgHeight should be >0");
  }
  return R_NilValue;
}

//' @title IFC_object Decompression to RAW
//' @name cpp_rawdecomp
//' @description
//' Operates decompression to raw of compressed image stored in TIFF file.
//' @param fname string, path to file.
//' @param offset std::size_t, position of the beginning of compressed image.
//' @param nbytes uint32_t, number of bytes of compressed image.
//' @param imgWidth uint32_t, Width of the decompressed image. Default is 1.
//' @param imgHeight uint32_t, Height of the decompressed image. Default is 1.
//' @param compression uint32_t, compression algorithm used. Default is 30818.
//' @param swap bool, whether to swap bytes or not. Default is false.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/ome/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{https://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}\cr
//' /verb{
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'    this list of conditions and the following disclaimer in the documentation
//'    and/or other materials provided with the distribution.
//' 
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' }
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector hpp_rawdecomp (const std::string fname, 
                               const std::size_t offset, 
                               const uint32_t nbytes, 
                               const uint32_t imgWidth = 1, 
                               const uint32_t imgHeight = 1, 
                               const uint32_t compression = 30818,
                               const bool swap = false,
                               const bool verbose = false) {
  Rcpp::RawVector raw_chnk = hpp_readchunk(fname, offset, nbytes, verbose);
  switch(compression) {
  case 1: return raw_chnk;
  case 30817: return hpp_gray_rawDecomp(raw_chnk, imgWidth, imgHeight, swap, verbose);
  case 30818: return hpp_rle_rawDecomp(raw_chnk, imgWidth, imgHeight, swap, verbose);
  }
  Rcpp::Rcerr << "hpp_rawdecomp: can't deal with compression format: " << compression << std::endl;
  Rcpp::stop("hpp_rawdecomp: can't deal with compression format");   
  return R_NilValue;
}

#endif
