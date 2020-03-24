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
