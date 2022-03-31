/*
  This file is released under the GNU General Public License, Version 3, GPL-3  
  Copyright (C) 2022 Yohann Demont                                              
                                                                                
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

#ifndef IFC_SPLIT_HPP
#define IFC_SPLIT_HPP

#include <Rcpp.h>
using namespace Rcpp;

// from https://gallery.rcpp.org/articles/fast-factor-generation
// slightly modified to handle raw and logical vectors
// it can also handle NA, NaN, NULL, -Inf, and +Inf
// depending on handleNA, NAs will be returned as NA or as an integer
template <int RTYPE>
Rcpp::IntegerVector fast_factor_T( const Rcpp::Vector<RTYPE>& x,
                                   const bool handleNA = true,
                                   const int from = 0) {
  Rcpp::Vector<RTYPE> lev = Rcpp::sort_unique(x);
  Rcpp::IntegerVector out = Rcpp::no_init(x.size());
  if(lev.size() == 1) {
    if(is_true(any(is_na(lev))) && handleNA) {
      out.fill(NA_INTEGER);
    } else {
      out.fill(1); 
    }
    if(from == LGLSXP) {
      out.attr("levels") = as<Rcpp::CharacterVector>(as<Rcpp::LogicalVector>(na_omit(lev)));
    } else {
      out.attr("levels") = as<Rcpp::CharacterVector>(Rcpp::na_omit(lev)); 
    }
  } else {
    out = Rcpp::match(x, lev);
    switch( from ) {
    case LGLSXP: {
      out.attr("levels") = as<Rcpp::CharacterVector>(as<Rcpp::LogicalVector>(na_omit(lev)));
    }
      break;
    case REALSXP: {
      Rcpp::CharacterVector ll = na_omit(as<Rcpp::CharacterVector>(lev));
      Rcpp::LogicalVector nan = is_nan(lev);
      if(is_true(any(nan))) {
        R_len_t bad = 0;
        for(R_len_t i = 0; i < ll.size(); i++) {
          if(ll[i] == "NaN") {
            bad = i + 1;
            break;
          }
        }
        Rcpp::LogicalVector nan_all = is_nan(x);
        for(R_len_t i = 0; i < out.size(); i++) if(nan_all[i]) out[i] = bad;
      }
      out.attr("levels") = ll;
      R_len_t k = lev.size();
      if(!handleNA && (lev.size() != ll.size())) for(R_len_t i = 0; i < out.size(); i++) if(out[i] == NA_INTEGER) out[i] = k;
    }
      break;
    default: {
      out.attr("levels") = as<Rcpp::CharacterVector>(Rcpp::na_omit(lev));
    }
    }
    if(from != REALSXP) {
      if(handleNA) {
        R_len_t bad = 0;
        for(R_len_t i = 0; i < lev.size(); i++) {
          if(lev[i] == NA_INTEGER) {
            bad = i + 1;
            break;
          }
        }
        if(bad > 0) for(R_len_t i = 0; i < out.size(); i++) if(out[i] == bad) out[i] = NA_INTEGER; 
      }
    }
  }
  out.attr("lvs") = lev.size();
  // no class attribution since if NA are present in x and handleNA is false,
  // printing the result will lead to 'malformed factor' error in R
  // out.attr("class") = "factor";
  return out;
}

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
//' @source adaptation from Kevin Ushey code \url{https://gallery.rcpp.org/articles/fast-factor-generation}
//' @return an IntegerVector, with attributes "levels" being the non-NA unique value(s) found
//' and "lvs" the total number of unique values found (including NA).
//' @keywords internal
////' @export
// [[Rcpp::export]]
SEXP hpp_fast_factor( SEXP x,
                      const bool handleNA = true) {
  switch( TYPEOF(x) ) {
  case NILSXP: return fast_factor_T(Rcpp::IntegerVector(0));
  case LGLSXP: return fast_factor_T(as<Rcpp::IntegerVector>(x), handleNA, LGLSXP);
  case INTSXP: return fast_factor_T<INTSXP>(x, handleNA);
  case REALSXP: return fast_factor_T<REALSXP>(x, handleNA, REALSXP);
  case STRSXP: return fast_factor_T<STRSXP>(x, handleNA);
  case RAWSXP: return fast_factor_T(as<Rcpp::IntegerVector>(x), handleNA);
  default: Rcpp::stop("hpp_fast_factor: not supported type in 'x'");
  }
}

//' @title Data Frame Groups Merging with Rcpp
//' @name cpp_group_df
//' @description
//' Computes global grouping factor from data.frame columns
//' @param df a DataFrame.
//' @return an IntegerVector with the resulting global grouping factor and a "table" attribute representing the amount of scalar in each resulting level.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::IntegerVector hpp_group_df(const Rcpp::DataFrame df) {
  R_len_t ll = 0;
  Rcpp::IntegerVector bar(df.nrows());                           // will hold global df[, ] factor
  for(R_len_t l = 0, i_col = 0; i_col < df.size(); i_col++) {
    Rcpp::IntegerVector foo = hpp_fast_factor(df[i_col], false); // hold current df[, i_col] factor
    l = foo.attr("lvs");
    if(l == 1) {
      foo.fill(0);
    } else {
      foo = foo - 1;
    }
    if(i_col == 0) {                                             // 1st loop over df columns
      bar = foo;                                                 // first copy to global factor
      ll = l;                                                    // first copy to global levels size
    } else {                                                     // concatenate factor / levels size
      bar += foo * ll;
      ll *= l;
    }
  }
  
  // create table count
  Rcpp::IntegerVector tab(ll);
  if(ll == 1) {
    tab[0] = bar.size();
  } else {
    for(Rcpp::IntegerVector::iterator i = bar.begin(); i != bar.end(); ++i) tab[*i]++;
  }
  // create output
  bar.attr("table") = tab;
  bar.attr("lvs") = ll;
  return bar;
}

#endif
