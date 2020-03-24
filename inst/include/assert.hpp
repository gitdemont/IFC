#ifndef IFC_ASSERT_HPP
#define IFC_ASSERT_HPP

#include <Rcpp.h>
Rcpp::LogicalVector cpp_assert(RObject x,
                               Rcpp::Nullable<IntegerVector> len = R_NilValue,
                               Rcpp::Nullable<CharacterVector> cla = R_NilValue,
                               Rcpp::Nullable<CharacterVector> typ = R_NilValue,
                               RObject alw = R_NilValue,
                               Rcpp::CharacterVector fun = "stop");

#endif
