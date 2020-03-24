#ifndef IFC_COLOR_HPP
#define IFC_COLOR_HPP

#include <Rcpp.h>
void cpp_HSV2RGB(double *r, double *g, double *b,
                 const double h = 0.0, const double s = 0.0, const double v = 0.0);
Rcpp::NumericVector cpp_M_HSV2RGB (const Rcpp::NumericMatrix mat, 
                                   const double h = 0.0,
                                   const double s = 0.0);

#endif
