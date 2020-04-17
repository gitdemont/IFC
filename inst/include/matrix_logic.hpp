#ifndef IFC_MATRIX_LOGIC_HPP
#define IFC_MATRIX_LOGIC_HPP

#include <Rcpp.h>

//' @title Matrix List And Logic
//' @name cpp_AND_M
//' @description
//' This function takes a list of matrices and returns the AND operation applied on these matrices.
//' @param list a list of logical matrices.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix cpp_AND_M(const Rcpp::List list) {
  R_len_t L = list.length();
  if(L < 1) Rcpp::stop("cpp_AND_M: 'list' should contain at least 1 matrix");
  Rcpp::LogicalMatrix MAT = Rcpp::clone(Rcpp::as<Rcpp::LogicalMatrix>(list[0]));
  R_len_t mat_r, mat_c; 
  mat_r = MAT.nrow();
  mat_c = MAT.ncol();
  if(L > 1) for(R_len_t i = 1; i < L; i++) {
    Rcpp::LogicalMatrix CUR_M = Rcpp::clone(Rcpp::as<Rcpp::LogicalMatrix>(list[i]));
    if(mat_r != CUR_M.nrow()) Rcpp::stop("cpp_AND_M: 'All matrices in 'list' should have same number of rows/columns");
    if(mat_c != CUR_M.ncol()) Rcpp::stop("cpp_AND_M: 'All matrices in 'list' should have same number of rows/columns");
    for(R_len_t i_row = 0; i_row < mat_r; i_row ++) {
      MAT(i_row, Rcpp::_) = CUR_M(i_row, Rcpp::_) & MAT(i_row, Rcpp::_);
    }
  }
  return MAT;
}

//' @title Matrix List Or Logic
//' @name cpp_OR_M
//' @description
//' This function takes a list of matrices and returns the OR operation applied on these matrices.
//' @param list a list of logical matrices.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix cpp_OR_M(const Rcpp::List list) {
  R_len_t L = list.length();
  if(L < 1) Rcpp::stop("cpp_OR_M: 'list' should contain at least 1 matrix");
  Rcpp::LogicalMatrix MAT = Rcpp::clone(Rcpp::as<Rcpp::LogicalMatrix>(list[0]));
  R_len_t mat_r, mat_c; 
  mat_r = MAT.nrow();
  mat_c = MAT.ncol();
  if(L > 1) for(R_len_t i = 1; i < L; i++) {
    Rcpp::LogicalMatrix CUR_M = Rcpp::clone(Rcpp::as<Rcpp::LogicalMatrix>(list[i]));
    if(mat_r != CUR_M.nrow()) Rcpp::stop("cpp_OR_M: 'All matrices in 'list' should have same number of rows/columns");
    if(mat_c != CUR_M.ncol()) Rcpp::stop("cpp_OR_M: 'All matrices in 'list' should have same number of rows/columns");
    for(R_len_t i_row = 0; i_row < mat_r; i_row ++) {
      MAT(i_row, Rcpp::_) = CUR_M(i_row, Rcpp::_) | MAT(i_row, Rcpp::_);
    }
  }
  return MAT;
}

//' @title Matrix Neg Logic
//' @name cpp_NEG_M
//' @description
//' This function takes a logical matrix and returns its negation.
//' @param mat LogicalMatrix.
//' @return a logical matrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix cpp_NEG_M(const Rcpp::LogicalMatrix mat) {
  Rcpp::LogicalMatrix OUT_M = Rcpp::no_init_matrix(mat.nrow(), mat.ncol());
  R_len_t mat_r = mat.nrow();
  for(R_len_t i_row = 0; i_row < mat_r; i_row++) {
    OUT_M(i_row, Rcpp::_) = !OUT_M(i_row, Rcpp::_);
  }
  return OUT_M;
}

#endif
