#ifndef IFC_GATE_HPP
#define IFC_GATE_HPP

#include <Rcpp.h>
Rcpp::NumericVector cpp_ell_coord (const Rcpp::NumericVector bound_x, 
                                   const Rcpp::NumericVector bound_y);
Rcpp::LogicalVector cpp_pnt_in_gate (const Rcpp::NumericMatrix pnts,
                                     const Rcpp::NumericMatrix gate, 
                                     const int algorithm = 1, 
                                     const double epsilon = 0.000000000001);

#endif
