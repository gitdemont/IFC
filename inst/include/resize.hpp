#ifndef IFC_RESIZE_HPP
#define IFC_RESIZE_HPP

#include <Rcpp.h>
Rcpp::NumericMatrix cpp_crop ( Rcpp::NumericMatrix mat,
                               const R_len_t new_height = 0,
                               const R_len_t new_width = 0);
  
Rcpp::NumericMatrix cpp_resize (const Rcpp::NumericMatrix mat, 
                                  const R_len_t new_height = 0, 
                                  const R_len_t new_width = 0,
                                  const bool add_noise = true, 
                                  const double bg = 0.0, const double sd = 0.0);

#endif
