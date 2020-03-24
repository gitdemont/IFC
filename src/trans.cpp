#include <Rcpp.h>
#include "../inst/include/utils.hpp"
using namespace Rcpp;

//' @title Smooth LinLog Transformation with Rcpp
//' @name cpp_smoothLinLog
//' @description
//' Takes a numeric vector and return its transformation:
//' - to linear, if abs(x) < hyper.
//' - to log, if abs(x) > hyper.
//' @param x NumericVector.
//' @param hyper double, value where transition between Lin/Log is applied.
//' @param base double, base of Log scale.
//' @param lin_comp double, value that is used to smooth transition between Lin/Log.
//' @keywords internal
////' @export
// [[Rcpp::export]]
NumericVector cpp_smoothLinLog (const NumericVector x, const double hyper = 1000.0, const double base = 10.0, const double lin_comp = 2.302585) {
  if(nNotisNULL(x)) {
    R_len_t L = x.size();
    double K = std::log(base) / lin_comp;
    NumericVector B = K + K * Rcpp::log((Rcpp::abs(x))/hyper);
    NumericVector xx = no_init_vector(L);
    for(R_len_t i = 0;i < L; i++) {
      if( std::fabs(x(i)) <= hyper ) {
        xx(i) = x(i) * K / hyper;
      } else {
        xx(i) = B(i);
        if(x(i)<0) xx(i) *= -1;
      }
    }
    return xx;
  } else {
    return x;
  }
}

//' @title Inverse Smooth LinLog Transformation with Rcpp
//' @name cpp_inv_smoothLinLog
//' @description
//' Takes a numeric vector and return its transformation:
//' - to linear, if abs(x) < log(base) / lin_comp.
//' - to exp, if abs(x) > log(base) / lin_comp.
//' @param x NumericVector.
//' @param hyper double, value where transition between Lin/Log is applied.
//' @param base double, base of Log scale.
//' @param lin_comp double, value that is used to smooth transition between Lin/Log.
//' @keywords internal
////' @export
// [[Rcpp::export]]
NumericVector cpp_inv_smoothLinLog (const NumericVector x, const double hyper = 1000.0, const double base = 10.0, const double lin_comp = 2.302585) {
  if(nNotisNULL(x)) {
    R_len_t L = x.size();
    double K = std::log(base) / lin_comp;
    NumericVector B = hyper * Rcpp::exp((Rcpp::abs(x))/K-1);
    NumericVector xx = no_init_vector(L);
    for(R_len_t i = 0;i < L; i++) {
      if( std::fabs(x(i)) > K ) {
        xx(i) = B(i);
        if(x(i)<0) xx(i) *= -1;
      } else {
        xx(i) = x(i) * hyper / K;
      }
    }
    return xx;
  } else {
    return x;
  }
}

//' @title Int32 to Uint32 32bits Conversion
//' @name cpp_int32_to_uint32
//' @description
//' Converts 32bits integer from signed to unsigned
//' @param x int32_t.
//' @keywords internal
////' @export
// [[Rcpp::export]]
uint32_t cpp_int32_to_uint32 ( const int32_t x) {
  uint32_t out = x;
  return out;
}

//' @title Uint32 to Int32 32bits Conversion
//' @name cpp_uint32_to_int32
//' @description
//' Converts 32bits integer from unsigned to signed
//' @param x uint32_t.
//' @keywords internal
////' @export
// [[Rcpp::export]]
int32_t cpp_uint32_to_int32 ( const uint32_t x) {
  int32_t out = x;
  return out;
}

//' @title Numeric to String Formatting
//' @name cpp_num_to_string
//' @description
//' Formats numeric to string used for features, images, ... values conversion when exporting to xml.
//' @param x NumericVector.
//' @param precision, number of significant decimal digits to keep when x < 1. Default is 16.
//' @return a string vector.
//' @keywords internal
////' @export
// [[Rcpp::export]]
Rcpp::StringVector cpp_num_to_string(const Rcpp::NumericVector x, const unsigned char precision = 16) {
  Rcpp::StringVector out(x.size());
  std::stringstream stream;
  std::string s;
  unsigned char relative;
  for(R_len_t i =0; i < x.size(); i++) {
    if(x(i) == R_PosInf) {
      out(i) = "+Inf";
      continue;
    }
    if(x(i) == R_NegInf) {
      out(i) = "-Inf";
      continue;
    }
    if(std::isnan(x[i])) { // for NaN or NA
      out(i) = "NaN";
      continue;
    }
    if(x[i] < 1) {
      relative = 0;
    } else {
      relative = std::log10(x[i]);
      if(x[i] > pow(10, relative)) {
        relative = relative + 1;
      } else {
        relative = relative + 2;
      }
    }
    stream << std::fixed << std::setprecision(precision - relative) << x(i);
    s = stream.str();
    s.erase(s.find_last_not_of('0') + 1, s.size());
    s.erase(s.find_last_not_of('.') + 1, s.size());
    out(i) = s;
    stream.str(std::string()); // clear stream
  }
  return out;
}
