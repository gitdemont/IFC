#ifndef IFC_MATRIX_LOGIC_HPP
#define IFC_MATRIX_LOGIC_HPP

#include <Rcpp.h>
Rcpp::LogicalMatrix cpp_AND_M(const Rcpp::List list);
Rcpp::LogicalMatrix cpp_OR_M(const Rcpp::List list);
Rcpp::LogicalMatrix cpp_NEG_M(const Rcpp::LogicalMatrix mat);

#endif
