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

#ifndef IFC_IMPORT_HPP
#define IFC_IMPORT_HPP

#include <Rcpp.h>
#include "utils.hpp"
using namespace Rcpp;

// import setPB from 'IFC' package
void hpp_setPB (SEXP pb,
                double value = 0.0,
                const std::string title = "",
                const std::string label = "") {
  Rcpp::Environment IFC = Rcpp::Environment::namespace_env("IFC");
  Rcpp::Function setPB = IFC["setPB"];
  if(SEXPsize(pb) != 0) setPB(pb, value, title, label);
}

// import basename from 'base' package
std::string hpp_basename (const std::string fname) {
  Rcpp::Environment base("package:base");
  Rcpp::Function basename = base["basename"];
  Rcpp::Nullable<Rcpp::CharacterVector> out_ = basename(fname);
  if(out_.isNotNull()) {
    Rcpp::CharacterVector out(out_.get());
    return as<std::string>(out[0]);
  }
  return "";
}

// use stats::quantile to compute input_range for matrix normalization
Rcpp::NumericVector hpp_quantile(const Rcpp::NumericMatrix mat,
                                 const double prob1 = 0.0,
                                 const double prob2 = 1.0) {
  Rcpp::Environment senv = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function sfun = senv["quantile"];
  return sfun(mat, Rcpp::Named("probs", Rcpp::NumericVector::create(prob1, prob2)));
}

//' @title Get Current Package Version
//' @name cpp_getversion
//' @description
//' Retrieve package version number.
//' @param pkg s std::string. Default is \code{"IFC"}
//' @return a Rcpp::CharacterVector
//' @keywords internal
////' @export
// [[Rcpp::export(rng = false)]]
Rcpp::CharacterVector hpp_getversion(const std::string pkg = "IFC") {
  Rcpp::Environment env = Rcpp::Environment::namespace_env(pkg);
  Rcpp::List L = env[".pkgenv"];
  Rcpp::List LL(0);
  Rcpp::CharacterVector sep = Rcpp::CharacterVector::create(".");
  if(L.containsElementNamed("version")) {
    LL = L["version"];
  } else {
    Rcpp::Environment uenv = Rcpp::Environment::namespace_env("utils");
    Rcpp::Function ufun = uenv["packageVersion"];
    LL = ufun(pkg);
  }
  if(LL.size()) {
    Rcpp::List LLL = LL[0];
    Rcpp::CharacterVector V(0);
    for(R_len_t i = 0; i < LLL.size(); i++) {
      if(V.size()) V = hpp_c(V, sep);
      V = hpp_c(V, Rcpp::CharacterVector::create(std::to_string(as<int>(LLL[i]))));
    }
    return V;
  }
  return 0;
}

#endif
