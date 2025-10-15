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

#define STRICT_R_HEADERS
#include <Rcpp.h>
#include "../inst/include/align.hpp"
#include "../inst/include/assert.hpp"
#include "../inst/include/base64.hpp"
#include "../inst/include/gate.hpp"
#include "../inst/include/utils.hpp"
#include "../inst/include/import.hpp"
#include "../inst/include/tiff.hpp"
#include "../inst/include/color.hpp"
#include "../inst/include/trans.hpp"
#include "../inst/include/scan.hpp"
#include "../inst/include/decomp.hpp"
#include "../inst/include/extract.hpp"
#include "../inst/include/resize.hpp"
#include "../inst/include/group.hpp"
#include "../inst/include/plot.hpp"
#include "../inst/include/fcs.hpp"
#include "../inst/include/tiffwrite.hpp"
using namespace Rcpp;

//' @title Get Current Compilation Bits Depth
//' @name cpp_getBits
//' @description
//' Retrieve size of std::size_t.
//' @return an unsigned integer
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
unsigned int cpp_getBits () {
  std::size_t foo = 0;
  return sizeof(foo);
}

// FROM align
//' @title Spatial Offsets Image Correction
//' @name cpp_align
//' @description
//' This function uses bilinear interpolation to apply spatial offset correction on image
//' @param mat, a NumericMatrix.
//' @param dx, a double x spatial offset. It has to be within ]-1,+1[. Default is NA_REAL for no change.
//' @param dy, a double y spatial offset. It has to be within ]-1,+1[. Default is NA_REAL for no change.
//' @details It is intended to be applied on raw images matrices from .rif files so has to generate spatial offset corrected image matrices.\cr
//' See William E. Ortyn et al. Sensitivity Measurement and Compensation in Spectral Imaging. Cytometry A 69 852-862 (2006).
//' \doi{10.1002/cyto.a.20306}
//' @return a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_align(const Rcpp::NumericMatrix mat,
                              const double dx = NA_REAL,
                              const double dy = NA_REAL) {
  return hpp_align(mat, dx, dy);
}
// END align

// FROM assert
//' @title Input Parameters Assertive Tool
//' @name cpp_assert
//' @description
//' Ensures that x respects several parameters
//' @param len IntegerVector, of allowed length for x. Default is NULL, for not checking this parameter.
//' @param cla CharacterVector, of allowed classes of x. Default is NULL, for not checking this parameter.
//' @param typ CharacterVector, of allowed types of x. Default is NULL, for not checking this parameter.
//' @param Robject of allowed values for x (will be passed to cpp_allowed). Default is NULL, for not checking this parameter.
//' @param fun CharacterVector, function to execute when mandatory parameters are not met. Default is "stop". Allowed are "stop","warning","message","return".
//' fun is placed in cpp_assert() in order to check it is correct before being used in assert() function.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalVector cpp_assert(const RObject x,
                               const Rcpp::Nullable<Rcpp::IntegerVector> len = R_NilValue,
                               const Rcpp::Nullable<Rcpp::CharacterVector> cla = R_NilValue,
                               const Rcpp::Nullable<Rcpp::CharacterVector> typ = R_NilValue,
                               const RObject alw = R_NilValue,
                               const Rcpp::CharacterVector fun = "stop") {
  return hpp_assert(x, len, cla, typ, alw, fun);
}
// END assert

// FROM gate
//' @title Ellipse Boundaries to Coordinates
//' @name cpp_ell_coord
//' @description
//' Transforms ellipse boundaries to usefull coordinates.
//' @param bound_x NumericVector, x-boundaries of the ellipse.
//' @param bound_y NumericVector, y-boundaries of the ellipse.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_ell_coord (const Rcpp::NumericVector bound_x,
                                   const Rcpp::NumericVector bound_y) {
  return hpp_ell_coord (bound_x, bound_y);
}

//' @title Point in Gate
//' @name cpp_pnt_in_gate
//' @description
//' This function checks if points lie in a polygon or ellipse.
//' @param pnts NumericMatrix, a 2-columns matrix with (x and y) coordinates of the points of interest.
//' @param gate NumericMatrix, a 2-columns matrix defining polygon vertices or ellipse boundaries.
//' @param algorithm int, used for computation. Default is 1.\cr
//' 1: Trigonometry.\cr
//' 2: Special case = axes-aligned rectangle.\cr
//' 3: Special case = axes-aligned ellipse.\cr
//' 4: Special case = axes-aligned rectangle not including upper limits.
//' @param epsilon double, epsilon threshold value. Default is 0.000000000001
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::LogicalVector cpp_pnt_in_gate (const Rcpp::NumericMatrix pnts,
                                     const Rcpp::NumericMatrix gate,
                                     const int algorithm = 1,
                                     const double epsilon = 0.000000000001) {
  return hpp_pnt_in_gate (pnts, gate, algorithm, epsilon);
}
// END gate

// FROM utils
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
bool cpp_mpfmatch(const Rcpp::Nullable<Rcpp::CharacterVector> x = R_NilValue,
                  const Rcpp::Nullable<Rcpp::CharacterVector> pattern = R_NilValue) {
  return hpp_mpfmatch(x, pattern);
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
R_len_t cpp_seqmatch(const Rcpp::StringVector x,
                     const Rcpp::StringVector y) {
  return hpp_seqmatch(x,y);
}
//' @title Use Rcpp to Apply Any on Matrix Rows
//' @name cpp_fast_rowAny
//' @description
//' Computes any across matrix rows.
//' @param M_ a Nullable LogicalVector. /!\ But cast to LogicalMatrix
//' @return a LogicalVector.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::LogicalVector> cpp_fast_rowAny(const Rcpp::Nullable<Rcpp::LogicalVector> M_ = R_NilValue) {
  return hpp_fast_rowAny(M_);
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
Rcpp::Nullable<Rcpp::LogicalVector> cpp_fast_listAny(const Rcpp::Nullable<Rcpp::List> L_ = R_NilValue) {
  return hpp_fast_listAny(L_);
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
Rcpp::NumericVector cpp_fast_range(const Rcpp::Nullable<Rcpp::NumericVector> x_ = R_NilValue) {
  return hpp_fast_range(x_);
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
Rcpp::IntegerVector cpp_fast_sample(const R_len_t n = 0,
                                    const R_len_t size = 0,
                                    const bool replace = false) {
  return hpp_fast_sample(n, size, replace);
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
Rcpp::Nullable<Rcpp::NumericVector> cpp_get_bytes_order (const R_len_t obj = 0,
                                                         const Rcpp::Nullable<Rcpp::IntegerVector> byt_ = R_NilValue,
                                                         const Rcpp::Nullable<Rcpp::IntegerVector> ord_ = R_NilValue,
                                                         const bool rev = false) {
  return hpp_get_bytes_order (obj, byt_, ord_, rev);
}

//' @title Non Finite Values Replacement
//' @name cpp_replace_non_finite
//' @description
//' This function replaces non finite values (NA, NaN -Inf and +Inf)
//' @param V_ a NumericVector.
//' @param by a double used as replcaement value. Default is 0.0
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::NumericVector> cpp_replace_non_finite(const Rcpp::Nullable<Rcpp::NumericVector> V_ = R_NilValue,
                                                           const double by = 0.0) {
  return hpp_replace_non_finite (V_, by);
}

//' @title Gamma Computation
//' @name cpp_computeGamma
//' @description
//' This function computes image gamma transformation value.
//' @param V named NumericVector of channel display properties containing 'xmin', 'xmax', 'xmid' and 'ymid'.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
double cpp_computeGamma (const Rcpp::NumericVector V) {
  return hpp_computeGamma (V);
}

//' @title Matrix to Matrix Writer According to Mask with Offsets
//' @name cpp_mark
//' @description
//' Writes matrix \code{'B'} in matrix \code{'A'} according to \code{'mask'}.
//' @param A a NumericMatrix.
//' @param B a NumericMatrix.
//' @param mask a NumericMatrix.
//' @param xoff x offset in \code{'A'} to start writing \code{'B'}.
//' @param yoff x offset in \code{'A'} to start writing \code{'B'}.
//' @param invert a logical. Default is \code{false}.
//' When \code{false}, the default, values of \code{'B'} are written into \code{'A'} when \code{'mask'} is not \code{0.0}.
//' When \code{true}, values of '\code{1-B}' are written into \code{'A'} when \code{'mask'} is not \code{0.0}.
//' @details indices resulting from writing B outside of A will trigger error.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_mark (const Rcpp::NumericMatrix A,
                              const Rcpp::NumericMatrix B,
                              const Rcpp::NumericMatrix mask,
                              const R_len_t xoff = 0,
                              const R_len_t yoff = 0,
                              const bool invert = false) {
  return hpp_mark (A, B, mask, xoff, yoff, invert);
}

//' @title Matrix to Matrix Writer According to Mask with Offsets
//' @name cpp_mark2
//' @description
//' Writes matrix \code{'B'} in matrix \code{'A'} according to \code{'mask'}.
//' @param A a NumericMatrix.
//' @param B a NumericMatrix.
//' @param mask a NumericMatrix.
//' @param xoff x offset in \code{'A'} to start writing \code{'B'}.
//' @param yoff x offset in \code{'A'} to start writing \code{'B'}.
//' @param invert a logical. Default is \code{false}.
//' When \code{false}, the default, values of \code{'B'} are written into \code{'A'} when \code{'mask'} is not \code{0.0}.
//' When \code{true}, values of '\code{1-B}' are written into \code{'A'} when \code{'mask'} is not \code{0.0}.
//' @details indices resulting from writing B outside of A will be skipped.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_mark2 (const Rcpp::NumericMatrix A,
                               const Rcpp::NumericMatrix B,
                               const Rcpp::NumericMatrix mask,
                               const R_len_t xoff = 0,
                               const R_len_t yoff = 0,
                               const bool invert = false) {
  return hpp_mark2 (A, B, mask, xoff, yoff, invert);
}

//' @title BMP Writer
//' @name cpp_writeBMP
//' @description
//' Transforms 3D [0,1] image to uncompressed bmp
//' @param image, a [0,1] normalized image matrix or 3D array. If 3D array, 3rd dimension should be of length 1 or 3.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector cpp_writeBMP (const Rcpp::NumericVector image) {
 return hpp_writeBMP (image); 
}

//' @title File Raw Chunk Extraction
//' @name cpp_readchunk
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
Rcpp::RawVector cpp_readchunk (const std::string fname, 
                               const std::size_t offset,
                               const uint32_t nbytes,
                               const bool verbose = false) {
  return hpp_readchunk(fname, offset, nbytes, verbose); 
}
// END utils

// FROM import
//' @title Get Current Package Version
//' @name cpp_getversion
//' @description
//' Retrieve package version number.
//' @param pkg s std::string. Default is \code{"IFC"}
//' @return a Rcpp::CharacterVector
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::CharacterVector cpp_getversion(const std::string pkg = "IFC") {
  return hpp_getversion(pkg);
}
// END import

// FROM base64
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
std::string cpp_base64_encode (const Rcpp::RawVector x, const bool url = false) {
 return hpp_base64_encode (x, url);
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
Rcpp::RawVector cpp_base64_decode(std::string x, const bool url = false) {
  return hpp_base64_decode (x, url);
}
// END base64

// FROM tiff
//' @title TIFF Checker
//' @name cpp_checkTIFF
//' @description
//' Checks if file is a TIFF.
//' @details If file is a TIFF it returns endianness of file, 'big' or 'little.
//' Otherwise, it shows an error and returns an empty string.
//' @param fname string, path to file.
//' @source TIFF 6.0 specifications archived from web \url{https://web.archive.org/web/20211209104854/https://www.adobe.io/open/standards/TIFF.html}
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
std::string cpp_checkTIFF (const std::string fname) {
  return hpp_checkTIFF(fname);
}

//' @title IFC_offsets Computation without Id Determination
//' @name cpp_getoffsets_noid
//' @description
//' Returns offsets of the IFDs (Image Field Directory) within a TIFF file.
//' @param fname string, path to file.
//' @param obj_count R_len_t, numbers of objects present in the file. Default is 0.
//' If obj_count <= 0 then progress_bar is forced to false.
//' @param display_progress bool, whether to display a progress bar. Default is false.
//' @param pb a List of class `IFC_progress` containing a progress bar of class `txtProgressBar`, `winProgressBar` or `Progress`. Default is R_Nilvalue.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source TIFF 6.0 specifications archived from web \url{https://web.archive.org/web/20211209104854/https://www.adobe.io/open/standards/TIFF.html}
//' @return an numeric vector with offsets of IFDs found.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_getoffsets_noid(const std::string fname,
                                        const R_len_t obj_count = 0,
                                        const bool display_progress = false,
                                        const Rcpp::Nullable<Rcpp::List> pb = R_NilValue,
                                        const bool verbose = false) {
  return hpp_getoffsets_noid(fname, obj_count, display_progress, pb, verbose);
}

//' @title IFD Tags Extraction
//' @name cpp_getTAGS
//' @description
//' Returns TAGS contained within an IFD (Image Field Directory) entry.
//' @param fname string, path to file.
//' @param offset std::size_t, position of the IFD beginning.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is 'false'.
//' @param trunc_bytes uint32_t maximal number of individual scalar to extract BYTE/ASCII/SBYTE/UNDEFINED for TAGS (1, 2, 6 or 7). Default is 12.\cr
//' However, if less is found, less is returned in map.
//' Note that, if 0 is provided, it will be automatically set to 1.
//' @param force_trunc whether to force truncation for all TAGS types. Default is FALSE.\cr
//' If 'true', 'trunc_bytes' will be used for TAGS (3, 4, 5, 8, 9, 10, 11 and 12) to extract desired number of individual scalar corresponding to each types.
//' @source TIFF 6.0 specifications archived from web \url{https://web.archive.org/web/20211209104854/https://www.adobe.io/open/standards/TIFF.html}
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_getTAGS (const std::string fname, 
                        const std::size_t offset, 
                        const bool verbose = false, 
                        const uint8_t trunc_bytes = 12, 
                        const bool force_trunc = false) {
  return hpp_getTAGS (fname, offset, verbose, trunc_bytes, force_trunc); 
}

//' @title IFD Fast Tags Extraction
//' @name cpp_fastTAGS
//' @description
//' Returns TAGS contained within an IFD (Image Field Directory) entry.
//' @param fname string, path to file.
//' @param offset std::size_t, position of the IFD beginning.
//' @param swap bool, whether to swap bytes or not. Default is 'false'.
//' @source TIFF 6.0 specifications archived from web \url{https://web.archive.org/web/20211209104854/https://www.adobe.io/open/standards/TIFF.html}
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_fastTAGS (const std::string fname, 
                         const std::size_t offset, 
                         const bool swap = false) {
  return hpp_fastTAGS (fname, offset, swap); 
}

//' @title IFC_offsets Computation with Object Identification
//' @name cpp_getoffsets_wid
//' @description
//' Returns offsets of the IFD (Image Field Directory) within a TIFF file.
//' @param fname string, path to file.
//' @param obj_count R_len_t, numbers of objects present in the file. Default is 0.
//' If obj_count <= 0 then progress_bar is forced to false.
//' @param display_progress bool, whether to display a progress bar. Default is false.
//' @param pb a List of class `IFC_progress` containing a progress bar of class `txtProgressBar`, `winProgressBar` or `Progress`. Default is R_Nilvalue.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source TIFF 6.0 specifications archived from web \url{https://web.archive.org/web/20211209104854/https://www.adobe.io/open/standards/TIFF.html}
//' @return a list of numeric vectors with OBJECT_ID, TYPE and OFFSET of IFDs found.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_getoffsets_wid(const std::string fname,
                              const R_len_t obj_count = 0,
                              const bool display_progress = false,
                              const Rcpp::Nullable<Rcpp::List> pb = R_NilValue,
                              const bool verbose = false) {
  return hpp_getoffsets_wid(fname, obj_count, display_progress, pb, verbose); 
}

//' @title Checksum for RIF/CIF
//' @name cpp_checksum
//' @description
//' Computes sum of img IFDs (Image Field Directory) offsets of objects 0, 1, 2, 3 and 4.
//' @param fname string, path to file.
//' @source TIFF 6.0 specifications archived from web \url{https://web.archive.org/web/20211209104854/https://www.adobe.io/open/standards/TIFF.html}
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
std::size_t cpp_checksum(const std::string fname) {
  return hpp_checksum(fname);
}
// END tiff

// FROM color
//' @title Matrix HSV to RGB Conversion
//' @name cpp_M_HSV2RGB
//' @description
//' Converts grayscale [0,1] mat to 3D rgb array according to hsv space.
//' hue and saturation determines tint whereas v is given by each element of mat
//' @param mat NumericMatrix, [0,1].
//' @param h double, [0,1], hue. Default is 0.0
//' @param s double, [0,1], saturation. Default is 0.0
//' @return a NumericVector with 3 dimensions attribute i.e. a 3D array
//' - 1st Dim is matrix rows count,
//' - 2nd Dim is matrix cols count,
//' - 3rd Dim is RGB
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_M_HSV2RGB (const Rcpp::NumericMatrix mat,
                                   const double h = 0.0,
                                   const double s = 0.0) {
  return hpp_M_HSV2RGB (mat, h, s);
}
// END color

// FROM trans
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
Rcpp::NumericVector cpp_smoothLinLog (const Rcpp::NumericVector x,
                                      const double hyper = 1000.0,
                                      const double base = 10.0,
                                      const double lin_comp = 2.302585) {
  return hpp_smoothLinLog (x, hyper, base, lin_comp);
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
Rcpp::NumericVector cpp_inv_smoothLinLog (const Rcpp::NumericVector x,
                                          const double hyper = 1000.0,
                                          const double base = 10.0,
                                          const double lin_comp = 2.302585) {
  return hpp_inv_smoothLinLog (x, hyper, base, lin_comp);
}

//' @title Uint32 to Raw Conversion
//' @name cpp_uint32_to_raw
//' @description
//' Converts unsigned 32bits integer to raw
//' @param x uint32_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector cpp_uint32_to_raw (const uint32_t x) {
  return hpp_uint32_to_raw (x);
}

//' @title Int32 to Uint32 32bits Conversion
//' @name cpp_int32_to_uint32
//' @description
//' Converts 32bits integer from signed to unsigned
//' @param x int32_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
uint32_t cpp_int32_to_uint32 (const int32_t x) {
  return hpp_int32_to_uint32 (x);
}

//' @title Uint32 to Int32 32bits Conversion
//' @name cpp_uint32_to_int32
//' @description
//' Converts 32bits integer from unsigned to signed
//' @param x uint32_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
int32_t cpp_uint32_to_int32 (const uint32_t x) {
  return hpp_uint32_to_int32 (x);
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
Rcpp::RawVector cpp_offset_to_raw (const double x, const bool swap = false) {
  return hpp_offset_to_raw(x, swap);
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
double cpp_raw_to_offset (const Rcpp::RawVector x, const bool swap = false) {
  return hpp_raw_to_offset(x, swap);
}

//' @title Int64 to Uint64 64bits Conversion
//' @name cpp_int64_to_uint64
//' @description
//' Converts 64bits integer from signed to unsigned
//' @param x int64_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
uint64_t cpp_int64_to_uint64 (const int64_t x) {
  return hpp_int64_to_uint64 (x);
}

//' @title Uint64 to Int64 64bits Conversion
//' @name cpp_uint64_to_int64
//' @description
//' Converts 64bits integer from unsigned to signed
//' @param x uint64_t.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
int64_t cpp_uint64_to_int64 (const uint64_t x) {
  return hpp_uint64_to_int64 (x);
}

//' @title Vectorize Int32 to Uint32 32bits Conversion
//' @name cpp_v_int32_to_uint32
//' @description
//' Converts 32bits vector of integers from unsigned to signed
//' @param V a NumericVector
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::NumericVector> cpp_v_int32_to_uint32 (Rcpp::Nullable<Rcpp::NumericVector> V = R_NilValue) {
  return hpp_v_int32_to_uint32(V);
}

//' @title Vectorize Int64 to Uint64 64bits Conversion
//' @name cpp_v_int64_to_uint64
//' @description
//' Converts 64bits vector of integers from unsigned to signed
//' @param V a NumericVector
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::Nullable<Rcpp::NumericVector> cpp_v_int64_to_uint64 (Rcpp::Nullable<Rcpp::NumericVector> V = R_NilValue) {
  return hpp_v_int64_to_uint64(V);
}
// END trans

// FROM scan
//' @title File Scanner
//' @name cpp_scanFirst
//' @description
//' Scans file for 1st occurence of a target string.
//' If found, it returns the position in bytes of the target.
//' Otherwise, it returns 0.
//' @param fname string, path to file.
//' @param raw a Rcpp::RawVector, exact target to be searched for. When converted to string it should be of at least 1 character and not exceed 1024 characters.
//' @param start size_t, position where to begin search.
//' It can't be superior or equal than file size or end (when end is different from 0 and inferior than file size).
//' @param end size_t, position where to stop searching. Default is 0.
//' Search will end up at this position unless it is higher than file size.
//' In such case, search will end up when file end will be reached.
//' @param buf_size uint8_t, size of buffer used to search for target (in kilo-Bytes, will be forced to be between 2 and 1024). Default is 64.
//' @return size_t index of first target character found within target plus 1 or 0 if not found.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
std::size_t cpp_scanFirst(const std::string fname, 
                          const Rcpp::RawVector raw, 
                          const std::size_t start = 0, 
                          const std::size_t end = 0, 
                          const uint8_t buf_size = 64) {
  return hpp_scanFirst(fname, raw, start, end, buf_size); 
}
// END scan

// FROM group
//' @title Fast Factorize Vector
//' @name cpp_fast_factor
//' @description
//' Makes factor with Rcpp
//' @param x a SEXP only NILSXP, LGLSXP, INTSXP, REALSXP, STRSXP, RAWSXP are handled.
//' @param handleNA a bool specifying if NA should be returned as NA or if they should be attributed a unique integer.
//' @details returned object should not be of class `factor` so has to prevent malformed factor' error when printing the result in R.
//' For this reason, it is aimed to be wrap in an R function to handle factor class / handleNA balance.\cr
//' e.g. either:
//' - structure(cpp_fast_factor(df, TRUE), class = "factor"), or,
//' - structure(cpp_fast_factor(unclass(df), FALSE), levels = NULL)
//' @source adaptation from Kevin Ushey code \url{https://gallery.rcpp.org/articles/fast-factor-generation/}
//' @return an IntegerVector, with attributes "levels" being the non-NA unique value(s) found
//' and "lvs" the total number of unique values found (including NA).
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
SEXP cpp_fast_factor( SEXP x, 
                      const bool handleNA = true) {
  return hpp_fast_factor(x, handleNA);
}

//' @title Data Frame Merge Groups with Rcpp
//' @name cpp_group_df
//' @description
//' Computes global grouping factor from data.frame columns
//' @param df a DataFrame.
//' @return an IntegerVector with the resulting global grouping factor and a "table" attribute representing the amount of scalar in each resulting level.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector cpp_group_df(const Rcpp::DataFrame df) {
  return hpp_group_df(df);
}
// END group

// FROM plot
//' @title Coordinates to Pixels
//' @name cpp_coord_to_px
//' @description low-level function to compute pixels coordinates
//' @param x NumericVector of x-coordinates of the points.
//' @param y NumericVector of y-coordinates of the points.
//' @param param NumericVector of parameters to scale raw points coordinates to pixels coordinates.
//' @return a 2 columns NumericMatrix of x and y pixels coordinates.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_coord_to_px (const Rcpp::NumericVector x,
                                     const Rcpp::NumericVector y,
                                     const Rcpp::NumericVector param) {
  return hpp_coord_to_px(x, y, param);
}

//' @title Image to Native Raster Conversion
//' @name cpp_as_nativeRaster
//' @description Converts 3D image array to nativeRaster
//' @param x an IntegerVecter /!\ It should be coercible to 3D array [height, width, rgba]
//' @return a nativeRaster IntegerMatrix
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix cpp_as_nativeRaster(const Rcpp::IntegerVector x) {
  return hpp_as_nativeRaster(x);
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
//' @param blur_size, a R_len_t the size of the gaussian blurring kernel. Default is 9.
//' @param blur_sd, a double the sd of the gaussian blurring kernel. Default is 3.0.
//' @details shape according to 'mask' will be drawn on 'img' centered at coordinates coords[, 1], coords[, 0]
//' and every pixels being part of the shape will be filled with 'color'.
//' If only one 'color' is provided, this 'color' will be used for each points.
//' If more than one 'color' is provided, then if number of colors (ncol) equals the number of points 'color' will be used as is for each single point.
//' Otherwise, 'color' will be considered as a color-gradient and density will be computed.
//' /!\ please note that IFC:::densCols() is faster to compute color based on density for n < 20000 points, so it's worth using it when number of points are lower.
//' @keywords internal
//' @return /!\ nothing is returned but img is modified in-place
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector cpp_draw (const Rcpp::IntegerVector img,
                              const Rcpp::IntegerMatrix coords = Rcpp::IntegerMatrix(1,2),
                              const Rcpp::LogicalMatrix mask = Rcpp::LogicalMatrix(1),
                              const Rcpp::IntegerMatrix color = Rcpp::IntegerMatrix(4,1),
                              const R_len_t blur_size = 9,
                              const double blur_sd = 3.0) {
  Rcpp::IntegerVector out = Rcpp::clone(img);
  hpp_draw(out, coords, mask, color, blur_size, blur_sd);
  return out;
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
//' @param bg_ a Nullable IntegerVector that will be cast to 3D array when not NULL. Default is R_NilValue.\cr
//' When not NULL, its dimensions should be the same as required by 'width' and 'height', otherwise an error will be thrown.\cr
//' When not NULL, it will serve as a background to draw new points on top of it.
//' @details shape according to 'pch' will be drawn centered at coordinates obj$coord[, 1], obj$coord[, 0]
//' and every pixels being part of the shape will be filled with 'color'.
//' If only one 'color' is provided, this 'color' will be used for each points.
//' If more than one 'color' is provided, then if number of colors (ncol) equals the number of points 'color' will be used as is for each single point.
//' Otherwise, 'color' will be considered as a color-gradient and density will be computed.
//' /!\ please note that IFC:::densCols() is faster to compute color based on density for n < 20000 points, so it's worth using it when number of points are lower.
//' @return an IntegerVector with dimensions [height, width, 4]
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector cpp_raster (const uint16_t width,
                                const uint16_t height,
                                const Rcpp::List obj,
                                const Rcpp::Nullable <Rcpp::IntegerVector> bg_ = R_NilValue) {
  return hpp_raster(width, height, obj, bg_);
}
// END plot

// FROM resize
//' @title Matrix Cropping
//' @name cpp_crop
//' @description
//' Crops mat according to new_height and new_width parameters.
//' @param mat a numeric matrix.
//' @param new_height an unsigned integer, giving the new height of returned mat. Default is 0 for no change.
//' @param new_width an unsigned integer, giving the new width of returned mat. Default is 0 for no change.
//' @return a cropped matrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_crop (Rcpp::NumericMatrix mat,
                              const R_len_t new_height = 0,
                              const R_len_t new_width = 0) {
  return hpp_crop(mat, new_height, new_width);
}

//' @title Matrix Resizing
//' @name cpp_resize
//' @description
//' Resizes mat according to new_height and new_width parameters.
//' @param mat a numeric matrix.
//' @param new_height an unsigned integer, giving the new height of returned mat. Default is 0 for no change.
//' @param new_width an unsigned integer, giving the new width of returned mat. Default is 0 for no change.
//' @param add_noise logical, if true adds normal noise when at least one new dimension is larger than original mat dimensions 
//' Rcpp::rnorm() function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.
//' @return a resized matrix with padding background if new_height or new_width is larger than original mat dimensions.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_resize (const Rcpp::NumericMatrix mat, 
                                const R_len_t new_height = 0, 
                                const R_len_t new_width = 0,
                                const bool add_noise = true, 
                                const double bg = 0.0, const double sd = 0.0) {
  return hpp_resize(mat, new_height, new_width, add_noise, bg, sd);
}
// END resize

// FROM decomp
//' @title IFC_object Decompression
//' @name cpp_decomp
//' @description
//' Operates decompression of compressed image stored in TIFF file.
//' @param fname string, path to file.
//' @param offset std::size_t, position of the beginning of compressed image.
//' @param nbytes uint32_t, number of bytes of compressed image.
//' @param imgWidth R_len_t, Width of the decompressed image. Default is 1.
//' @param imgHeight R_len_t, Height of the decompressed image. Default is 1.
//' @param nb_channels R_len_t, number of channels of the decompressed image. Default is 1.
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
//' \verb{
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
Rcpp::List cpp_decomp (const std::string fname, 
                       const std::size_t offset, 
                       const uint32_t nbytes, 
                       const uint32_t imgWidth = 1, 
                       const uint32_t imgHeight = 1, 
                       const uint32_t nb_channels = 1,
                       const uint8_t removal = 0,
                       const uint32_t compression = 30818,
                       const bool verbose = false) {
  return hpp_decomp (fname, offset, nbytes, imgWidth, imgHeight, nb_channels, removal, compression, verbose);
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
//' \verb{
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
Rcpp::RawVector cpp_rawdecomp (const std::string fname, 
                               const std::size_t offset, 
                               const uint32_t nbytes, 
                               const uint32_t imgWidth = 1, 
                               const uint32_t imgHeight = 1, 
                               const uint32_t compression = 30818,
                               const bool swap = false,
                               const bool verbose = false) {
  return hpp_rawdecomp (fname, offset, nbytes, imgWidth, imgHeight, compression, swap, verbose);
}
// END decomp


// FROM extract
//' @title Matrix Normalization
//' @name cpp_normalize
//' @description
//' Normalizes a finite matrix to [0,1]
//' @param mat a finite NumericMatrix.
//' @param input_range a finite NumericVector, sets the range of the input intensity values. Values outside this range are clipped. Default is \code{[0.0,4095.0]}.
//' @param full_range if \code{'full_range'} is \code{true}, then \code{'input_range'} will be set to \code{[0.0,4095.0]} and \code{'gamma'} forced to \code{1.0}. Default is \code{false}.
//' @param force_range if \code{'force_range'} is \code{true}, then \code{'input_range'} will be adjusted to \code{'mat'} range in \code{[-4095.0, +inf]} and \code{'gamma'} forced to \code{1.0}. Default is \code{false}.\cr
//' Note that this parameter takes the precedence over \code{'input_range'} and \code{'full_range'}.
//' @param gamma correction. Default is \code{1.0}, for no correction.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_normalize (const Rcpp::NumericMatrix mat, 
                                   const Rcpp::NumericVector input_range = Rcpp::NumericVector::create(0.0,4095.0),
                                   const bool full_range = false,
                                   const bool force_range = false, 
                                   const double gamma = 1.0) {
  return hpp_normalize (mat, input_range, full_range, force_range, gamma);
}

//' @title Matrix Cleanser
//' @name cpp_cleanse
//' @description
//' Replaces values in matrix mat according to mask msk.
//' Depending of \code{'add_noise'} parameter, values of \code{'mat'} will be replaced with noise or not.
//' @param mat a NumericMatrix.
//' @param msk a IntegerMatrix.
//' @param add_noise bool, whether to add normal noise or not, \code{Rcpp::Rf_rnorm(bg, sd)} function is used. Default is \code{true}.
//' @param bg double, mean value of the background added if \code{'add_noise'} is \code{true}. Default is \code{0.0}.
//' @param sd double, standard deviation of the background added if \code{'add_noise'} is \code{true}. Default is \code{0.0}.
//' @return a NumericMatrix cleansed according to \code{'msk'}.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_cleanse (const Rcpp::NumericMatrix mat, 
                                 const Rcpp::IntegerMatrix msk,
                                 const bool add_noise = true, 
                                 const double bg = 0.0, const double sd = 0.0) {
  return hpp_cleanse (mat, msk, add_noise, bg, sd);
}

//' @title Equal Sized Matrix to Matrix Writer According to Mask
//' @name cpp_mask
//' @description
//' Writes matrix \code{'B'} in matrix \code{'A'} according to \code{'mask'}.
//' If \code{'mask'} is not \code{0.0} \code{'B'} is written, \code{'A'} otherwise.
//' @param A a NumericMatrix.
//' @param B a NumericMatrix.
//' @param mask a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix cpp_mask (const Rcpp::NumericMatrix A,
                              const Rcpp::NumericMatrix B,
                              const Rcpp::NumericMatrix mask) {
  return hpp_mask (A, B, mask);
}

//' @title Matrix Transformation
//' @name cpp_transform
//' @description
//' Function to normalize, colorize and add background to images.
//' @param mat NumericMatrix.
//' @param color NumericVector, whose members are h,s,v color.
//' This vector has to be named with 1st name being the name of this color.
//' @param msk IntegerMatrix.
//' @param size a length 2 IntegerVector, of final dimensions (height,width) of the image. Default is \code{0,0} for no change.
//' @param mode string, color mode export. Either \code{"rgb"}, \code{"gray"} or \code{"raw"}. Default is \code{"raw"}.
//' @param type uint16_t image object type. Default is \code{2}.
//' @param input_range a finite NumericVector, only apply when \code{'mode'} is not \code{"raw"}, sets the range of the input intensity values. Values outside this range are clipped. Default is \code{[0.0,4095.0]}.
//' @param add_noise bool, whether to add normal noise or not. Default is \code{true}.
//' @param bg double, mean value of the background added if \code{'add_noise'} is \code{true}. Default is \code{0.0}.
//' @param sd double, standard deviation of the background added if \code{'add_noise'} is \code{true}. Default is \code{0.0}.
//' @param full_range bool, only apply when \code{'mode'} is not \code{"raw"}, if \code{'full_range'} is \code{true}, then \code{'input_range'} will be set to \code{[0.0,4095.0]} and \code{'gamma'} forced to \code{1.0}. Default is \code{false}.
//' @param force_range bool, only apply when \code{'mode'} is not \code{"raw"}, if \code{'force_range'} is \code{true}, then \code{'input_range'} will be adjusted to object range in \code{[-4095.0, +inf]} and \code{'gamma'} forced to \code{1.0}. Default is \code{false}.\cr
//' Note that this parameter takes the precedence over \code{'input_range'} and \code{'full_range'}.
//' @param gamma correction. Default is \code{1.0}, for no correction.
//' @param spatialX X offset correction. Default is \code{0.0} for no change.
//' @param spatialY Y offset correction. Default is \code{0.0} for no change.
//' @details When a mask is detected, \code{'add_noise'}, \code{'full_range'} and \code{'force_range'} are set to \code{false}, \code{'bg'} and \code{'sd'} to \code{0.0}, \code{'input_range'} to \code{[0.0,3.0]} and \code{'gamma'} to \code{1.0}.\cr\cr
//' Experimental (as of v0.2.0.501): when \code{'mode'} is not \code{"raw"}, if \code{'input_range'} is within \code{]0,1[}, it will be used as \code{'probs'} argument to \link[stats]{quantile} to determine clipping range from image object, in addition, \code{'gamma'} will be forced to \code{1.0}.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::NumericVector cpp_transform(const Rcpp::NumericMatrix mat,
                                  const Rcpp::NumericVector color,
                                  const Rcpp::IntegerMatrix msk,
                                  const Rcpp::IntegerVector size = Rcpp::IntegerVector::create(0,0),
                                  const std::string mode = "raw",
                                  const uint16_t type = 2,
                                  const Rcpp::NumericVector input_range = Rcpp::NumericVector::create(0.0,4095.0),
                                  const bool add_noise = true,
                                  const double bg = 0.0,
                                  const double sd = 0.0,
                                  const bool full_range = false,
                                  const bool force_range = false,
                                  const double gamma = 1.0,
                                  const double spatialX = 0.0,
                                  const double spatialY = 0.0) {
  return hpp_transform(mat, color, msk, size, mode, type, input_range, add_noise, bg, sd, full_range, force_range, gamma, spatialX, spatialY); 
}

//' @title IFC_object Extraction
//' @name cpp_extract
//' @description
//' Extracts object from ifd
//' @param fname string, path to file
//' @param ifd List, ifd information of class IFC_ifd
//' @param colors List of colors to use.
//' @param channels DataFrame, channels information.
//' @param physicalChannel CharacterVector of indices for each channel.
//' @param xmin NumericVector of minimal values for each channel.
//' @param xmax NumericVector of maximal values for each channel.
//' @param spatialX NumericVector of X spatial offset correction for each channel.
//' @param spatialY NumericVector of Y spatial offset correction for each channel.
//' @param removal IntegerVector of removal method to be used for each channel.
//' @param add_noise LogicalVector of whether to \code{'add_noise'} for each channel.
//' @param full_range LogicalVector of whether to use \code{'full_range'} for each channel.
//' @param force_range LogicalVector of whether to use \code{'force_range'} for each channel.
//' @param gamma NumericVector of the \code{'gamma'} for each channel.
//' @param chan_to_extract IntegerVector, channels to extract.
//' @param extract_msk uint8_t, type of mask to extract:\cr
//' -\code{0}: no mask,\cr
//' -\code{1}: at least one \code{"raw"},\cr
//' -\code{2}: at least one \code{"clipped"},\cr
//' -\code{3}: at least one \code{"masked"},\cr
//' -\code{4}: at least one \code{"MC"} mask.
//' @param mode string, color mode export. Either \code{"rgb"}, \code{"gray"} or \code{"raw"}. Default is \code{"raw"}.
//' @param size a length 2 IntegerVector of final dimensions (height,width) of the image. Default is \code{0,0} for no change.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is \code{false}.
//' @details When \code{'add_noise'} is \code{false}, background will be automatically set to minimal pixel value for \code{"masked"} and \code{"MC"} \code{'removal'} method.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::List cpp_extract (const std::string fname,
                        const Rcpp::List ifd,
                        const Rcpp::List colors,
                        const Rcpp::CharacterVector physicalChannel,
                        const Rcpp::NumericVector xmin,
                        const Rcpp::NumericVector xmax,
                        const Rcpp::NumericVector spatialX,
                        const Rcpp::NumericVector spatialY,
                        const Rcpp::IntegerVector removal,
                        const Rcpp::LogicalVector add_noise,
                        const Rcpp::LogicalVector full_range,
                        const Rcpp::LogicalVector force_range,
                        const Rcpp::NumericVector gamma,
                        const Rcpp::IntegerVector chan_to_extract,
                        const uint8_t extract_msk = 0,
                        const std::string mode = "raw",
                        const Rcpp::IntegerVector size = Rcpp::IntegerVector::create(0,0),
                        const bool verbose = false) {
  return hpp_extract (fname, ifd, colors, physicalChannel, xmin, xmax, spatialX, spatialY, removal, add_noise, full_range, force_range, gamma, chan_to_extract, extract_msk, mode, size, verbose);
}
// END extract

// FROM fcs
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
Rcpp::NumericVector cpp_readFCS (const std::string fname,
                                 const std::size_t offset,
                                 const uint32_t events,
                                 const Rcpp::IntegerVector b_typ,
                                 const Rcpp::IntegerVector b_siz,
                                 const Rcpp::Nullable<Rcpp::IntegerVector> b_msk = R_NilValue,
                                 const bool swap = false) {
  return hpp_readFCS(fname, offset, events, b_typ, b_siz, b_msk, swap);
 }
// END fcs

// FROM tiffwrite
//' @title Cast Image
//' @name cpp_cast_image
//' @description
//' Casts image from RAW, INT or REAL SEXP vector
//' @param x SEXP, the image to cast
//' @param what uint8_t, type to use for casting. Default is \code{1}.
//' Allowed are 1=uint8_t, 2=int8_t, 3=uint16_t, 4=int16_t, 5=uint32_t, 6=int32_t, 7=float_t, 8=double_t.
//' @param swap bool, whether single scalar values of \code{x} should be swap. Default is \code{false}.
//' @return a Rcpp::RawVector
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector cpp_cast_image ( SEXP x,
                                 const uint8_t what = 1,
                                 const bool swap = false) {
  return hpp_cast_image(x, what, swap);
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
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::RawVector cpp_writeIFD ( const Rcpp::RawVector img,
                               const Rcpp::List tags,
                               const uint32_t offset = 0,
                               const std::string endianness = "little",
                               const bool rgb = false,
                               const bool last = false,
                               const bool verbose = false) {
  return hpp_writeIFD(img, tags, offset, endianness, rgb, last, verbose);
}
// END tiffwrite
