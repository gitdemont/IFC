/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2020 Yohann Demont                                              
                                                                                
  It is part of IFC package, please cite:                                       
  -IFC: An R Package for Imaging Flow Cytometry                                 
  -YEAR: 2020                                                                   
  -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,              
                      Jean-Pierre Marolleau, Loï�c Gaççon,                       
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
// template <typename T> bool NotisNULL(Rcpp::Nullable<T> x_ = R_NilValue); //can't deduce type of list
bool nNotisNULL(Rcpp::Nullable<Rcpp::NumericVector> x_ = R_NilValue);
bool iNotisNULL(Rcpp::Nullable<Rcpp::IntegerVector> x_ = R_NilValue);
std::string to_string(uint16_t x);
Rcpp::NumericVector cpp_check_range(const Rcpp::NumericVector x);
double cpp_computeGamma (const Rcpp::NumericVector V);
std::string cpp_base64_encode ( const Rcpp::RawVector x);
Rcpp::RawVector cpp_writeBMP ( const Rcpp::NumericVector image);

#endif
