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

#ifndef IFC_CBIND_HPP
#define IFC_CBIND_HPP

#include <Rcpp.h>
#include "utils.hpp"
using namespace Rcpp;

//' @title Fast Dataframe and Matrix Column Binding
//' @name cpp_fast_cbind_DF_M
//' @description
//' Combines data.frame and matrix by columns
//' @param Df_ a Nullable DataFrame.
//' @param M_ a Nullable NumericVector. /!\ But cast to NumericMatrix.
//' @param add_id a bool determining if 1st column of returned object should be given 1 to nrow integers
//' @return a DataFrame.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::DataFrame> hpp_fast_cbind_DF_M(const Rcpp::Nullable<Rcpp::DataFrame> Df_ = R_NilValue,
                                                    const Rcpp::Nullable<Rcpp::NumericVector> M_ = R_NilValue,
                                                    const bool add_id = false) {
  if(Df_.isNotNull()) {
    Rcpp::DataFrame Df(Df_.get());
    if(nNotisNULL(M_)) {
      Rcpp::NumericMatrix M(M_.get());
      if(Df.nrows() == M.nrow()) {
        Rcpp::DataFrame out = clone(Df);
        for(R_len_t i_col = 0; i_col < M.ncol(); i_col++) out.push_back(M(Rcpp::_, i_col));
        if(add_id) {
          Rcpp::IntegerVector ids = seq_along(M(Rcpp::_, 0));
          out.push_front(ids);
        }
        return out;
      } else {
        Rcpp::stop("hpp_fast_cbind_DF_M: 'Df' and 'M' should have same number of rows");
      }
    } 
  } else {
    if(nNotisNULL(M_)) {
      Rcpp::NumericMatrix M(M_.get());
      Rcpp::DataFrame out = as<Rcpp::DataFrame>(M);
      if(add_id) {
        Rcpp::IntegerVector ids = seq_along(M(Rcpp::_, 0));
        out.push_front(ids);
      }
      return out;
    }
  }
  return Df_;
}

//' @title Fast Matrix and Dataframe Column Binding
//' @name cpp_fast_cbind_M_DF
//' @description
//' Combines matrix and data.frame by columns
//' @param M_ a Nullable NumericVector. /!\ But cast to NumericMatrix.
//' @param Df_ a Nullable DataFrame.
//' @param add_id a bool determining if 1st column of returned object should be given 1 to nrow integers
//' @return a DataFrame.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::DataFrame> hpp_fast_cbind_M_DF(const Rcpp::Nullable<Rcpp::NumericVector> M_ = R_NilValue,
                                                    const Rcpp::Nullable<Rcpp::DataFrame> Df_ = R_NilValue,
                                                    const bool add_id = false) {
  if(Df_.isNotNull()) {
    Rcpp::DataFrame Df(Df_.get());
    if(nNotisNULL(M_)) {
      Rcpp::NumericMatrix M(M_.get());
      if(Df.nrows() == M.nrow()) {
        Rcpp::DataFrame out = clone(Df);
        for(R_len_t i_col = M.ncol() - 1; i_col >= 0; i_col--) out.push_front(M(Rcpp::_, i_col));
        if(add_id) {
          Rcpp::IntegerVector ids = seq_along(M(Rcpp::_, 0));
          out.push_front(ids);
        }
        return out;
      } else {
        Rcpp::stop("hpp_fast_cbind_M_DF: 'M' and 'Df' should have same number of rows");
      }
    }
  } else {
    if(nNotisNULL(M_)) {
      Rcpp::NumericMatrix M(M_.get());
      Rcpp::DataFrame out = as<Rcpp::DataFrame>(M);
      if(add_id) {
        Rcpp::IntegerVector ids = seq_along(M(Rcpp::_, 0));
        out.push_front(ids);
      }
      return out;
    }
  }
  return Df_;
}

//' @title Dataframe and Dataframe Column Binding
//' @name cpp_fast_cbind_DF_DF
//' @description
//' Combines numeric matrix by columns
//' @param Df1_ a Nullable DataFrame.
//' @param Df2_ a Nullable DataFrame.
//' @param add_id a bool determining if 1st column of returned object should be given 1 to nrow integers
//' @return a DataFrame.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::DataFrame> hpp_fast_cbind_DF_DF(const Rcpp::Nullable<Rcpp::DataFrame> Df1_ = R_NilValue,
                                                     const Rcpp::Nullable<Rcpp::DataFrame> Df2_ = R_NilValue,
                                                     const bool add_id = false) {
  if(Df1_.isNotNull() && Df2_.isNotNull() ) {
    Rcpp::DataFrame Df1(Df1_.get());
    Rcpp::DataFrame Df2(Df2_.get());
    if(Df1.nrows() == Df2.nrows()) {
      Rcpp::DataFrame out = clone(Df1);
      for(R_len_t i_col = 0; i_col < Df2.size(); i_col++) out.push_back(Df2[i_col]);
      if(add_id) {
        Rcpp::IntegerVector ids = Rcpp::no_init(out.nrows());
        for(R_len_t i = 0; i < ids.size(); i++) ids[i] = i + 1;
        out.push_front(ids);
      }
      return out;
    } else {
      Rcpp::stop("hpp_fast_cbind_DF_DF: 'Df1' and 'Df2' should have same number of rows");
    }
  }
  if(Df1_.isNotNull()) {
    Rcpp::DataFrame out(Df1_.get());
    if(add_id) {
      Rcpp::IntegerVector ids = Rcpp::no_init(out.nrows());
      for(R_len_t i = 0; i < ids.size(); i++) ids[i] = i + 1;
      out.push_front(ids);
    }
    return out;
  }
  if(Df2_.isNotNull()) {
    Rcpp::DataFrame out(Df2_.get());
    if(add_id) {
      Rcpp::IntegerVector ids = Rcpp::no_init(out.nrows());
      for(R_len_t i = 0; i < ids.size(); i++) ids[i] = i + 1;
      out.push_front(ids);
    }
    return out;
  }
  return Df1_;
}

//' @title Matrix and Matrix Column Binding
//' @name cpp_fast_cbind_M_M
//' @description
//' Combines numeric matrix by columns
//' @param M1_ a Nullable NumericVector. /!\ But cast to NumericMatrix.
//' @param M2_ a Nullable NumericVector. /!\ But cast to NumericMatrix.
//' @param add_id a bool determining if 1st column of returned object should be given 1 to nrow integers
//' @return a NumericVector.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::NumericVector> hpp_fast_cbind_M_M(const Rcpp::Nullable<Rcpp::NumericVector> M1_ = R_NilValue,
                                                       const Rcpp::Nullable<Rcpp::NumericVector> M2_ = R_NilValue,
                                                       const bool add_id = false) {
  if(nNotisNULL(M1_) && nNotisNULL(M2_)) {
    Rcpp::NumericMatrix M1(M1_.get());
    Rcpp::NumericMatrix M2(M2_.get());
    R_len_t mat_r = M1.nrow(), mat_c = M1.ncol();
    if(mat_r != M2.nrow()) Rcpp::stop("hpp_fast_cbind_M_M: 'M1' and 'M2' should have same number of rows");
    Rcpp::NumericMatrix out = no_init(mat_r, mat_c + M2.ncol() + add_id);
    if(add_id) out(Rcpp::_, 0) = seq_along(out(Rcpp::_, 0));
    for(R_len_t i_col = 0; i_col < mat_c; i_col++) out(Rcpp::_, i_col + add_id) = M1(Rcpp::_, i_col);
    for(R_len_t i_col = 0; i_col < M2.ncol(); i_col++) out(Rcpp::_, i_col + add_id + mat_c) = M2(Rcpp::_, i_col);
    return out;
  }
  if(nNotisNULL(M2_)) {
    Rcpp::NumericMatrix out(M2_.get());
    if(add_id) out(Rcpp::_, 0) = seq_along(out(Rcpp::_, 0));
    return out;
  }
  if(nNotisNULL(M1_)) {
    Rcpp::NumericMatrix out(M1_.get());
    if(add_id) out(Rcpp::_, 0) = seq_along(out(Rcpp::_, 0));
    return out;
  }
  return M1_;
}

//' @title Fast Dataframe and List Column Binding
//' @name cpp_fast_cbind_DF_L
//' @description
//' Combines data.frame and list by columns
//' @param Df_ a Nullable DataFrame.
//' @param L_ a Nullable List.
//' @param add_id a bool determining if 1st column of returned object should be given 1 to nrow integers
//' @return a DataFrame.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::DataFrame> hpp_fast_cbind_DF_L(const Rcpp::Nullable<Rcpp::DataFrame> Df_ = R_NilValue,
                                                    const Rcpp::Nullable<Rcpp::List> L_ = R_NilValue,
                                                    const bool add_id = false) {
  if(Df_.isNotNull()) {
    Rcpp::DataFrame Df(Df_.get());
    if(L_.isNotNull()) {
      Rcpp::List L(L_.get());
      R_len_t k = -1;
      for(R_len_t i = 0; i < L.size(); i++) {
        if(k >= 0) {
          if(k != Rf_length(L[i])) Rcpp::stop("hpp_fast_cbind_DF_L: 'Df' and 'L' should have same number of rows");
        } else {
          k = Rf_length(L[i]);
        }
      }
      Rcpp::DataFrame out = clone(Df);
      for(R_len_t i_col = 0; i_col < L.size(); i_col++) out.push_back(L[i_col]);
      if(add_id) {
        Rcpp::IntegerVector ids = no_init(out.nrows());
        for(R_len_t i = 0; i < ids.size(); i++) ids[i] = i + 1; 
        out.push_front(ids);
      }
      return out;
    } else {
      if(L_.isNotNull()) {
        Rcpp::List L(L_.get());
        Rcpp::DataFrame out = as<Rcpp::DataFrame>(L);
        if(add_id) {
          Rcpp::IntegerVector ids = no_init(out.nrows());
          for(R_len_t i = 0; i < ids.size(); i++) ids[i] = i + 1; 
          out.push_front(ids);
        }
        return out;
      }
    }
  }
  return Df_;
}

//' @title Fast List and Dataframe Column Binding
//' @name cpp_fast_cbind_L_DF
//' @description
//' Combines list and data.frame by columns
//' @param L_ a Nullable List.
//' @param Df_ a Nullable DataFrame.
//' @param add_id a bool determining if 1st column of returned object should be given 1 to nrow integers
//' @return a DataFrame.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::Nullable<Rcpp::DataFrame> hpp_fast_cbind_L_Df(const Rcpp::Nullable<Rcpp::List> L_ = R_NilValue,
                                                    const Rcpp::Nullable<Rcpp::DataFrame> Df_ = R_NilValue,
                                                    const bool add_id = false) {
  if(Df_.isNotNull()) {
    Rcpp::DataFrame Df(Df_.get());
    if(L_.isNotNull()) {
      Rcpp::List L(L_.get());
      R_len_t k = -1;
      for(R_len_t i = 0; i < L.size(); i++) {
        if(k >= 0) {
          if(k != Rf_length(L[i])) Rcpp::stop("hpp_fast_cbind_L_Df: 'Df' and 'L' should have same number of rows");
        } else {
          k = Rf_length(L[i]);
        }
      }
      Rcpp::DataFrame out = clone(Df);
      for(R_len_t i_col = L.size() - 1; i_col >=0 ; i_col--) out.push_front(L[i_col]);
      if(add_id) {
        Rcpp::IntegerVector ids = no_init(out.nrows());
        for(R_len_t i = 0; i < ids.size(); i++) ids[i] = i + 1; 
        out.push_front(ids);
      }
      return out;
    } else {
      if(L_.isNotNull()) {
        Rcpp::List L(L_.get());
        Rcpp::DataFrame out = as<Rcpp::DataFrame>(L);
        if(add_id) {
          Rcpp::IntegerVector ids = no_init(out.nrows());
          for(R_len_t i = 0; i < ids.size(); i++) ids[i] = i + 1; 
          out.push_front(ids);
        }
        return out;
      }
    }
  }
  return Df_;
}

#endif
