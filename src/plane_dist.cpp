#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix c_plane_dist(NumericVector x, NumericVector y, IntegerVector from, IntegerVector to) {
  int k, l, i, j, n = x.length();
  double d;
  NumericMatrix D(n,n);
  for(k=0; k < from.length(); k++){
    i = from[k]-1;
    for(l=0; l < to.length(); l++) {
      j = to[l]-1;
      d = sqrt( pow(x[i]-x[j]  ,2)+pow(y[i]-y[j], 2)   );
      D(i,j) = d;
    }
  }
  return D;
}