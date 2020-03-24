#ifndef IFC_TRANS_HPP
#define IFC_TRANS_HPP

#include <Rcpp.h>
Rcpp::NumericVector cpp_smoothLinLog (const Rcpp::NumericVector x, 
                                      const double hyper = 1000.0, 
                                      const double base = 10.0, 
                                      const double lin_comp = 2.302585);
Rcpp::NumericVector cpp_inv_smoothLinLog (const Rcpp::NumericVector x, 
                                          const double hyper = 1000.0, 
                                          const double base = 10.0, 
                                          const double lin_comp = 2.302585);
uint32_t cpp_int32_to_uint32 ( const int32_t x);
int32_t cpp_uint32_to_int32 ( const uint32_t x);
Rcpp::StringVector cpp_num_to_string(const Rcpp::NumericVector x, 
                                     const unsigned char precision = 16);

#endif
