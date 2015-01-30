// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// c_plane_dist
NumericMatrix c_plane_dist(NumericVector x, NumericVector y, IntegerVector from, IntegerVector to);
RcppExport SEXP sphere_c_plane_dist(SEXP xSEXP, SEXP ySEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP );
        NumericMatrix __result = c_plane_dist(x, y, from, to);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// c_sphere_dist
NumericMatrix c_sphere_dist(NumericVector lat, NumericVector lon, IntegerVector from, IntegerVector to);
RcppExport SEXP sphere_c_sphere_dist(SEXP latSEXP, SEXP lonSEXP, SEXP fromSEXP, SEXP toSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type lat(latSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type lon(lonSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type from(fromSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type to(toSEXP );
        NumericMatrix __result = c_sphere_dist(lat, lon, from, to);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}