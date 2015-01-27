#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

double min(double a, double b){
  if(a < b) return a;
  return b;
}

// [[Rcpp::export]]
NumericMatrix c_sphere_dist(NumericVector lat, NumericVector lon, IntegerVector from, IntegerVector to) {
  int k, l, i, j, n = lon.length();
  double dlat, dlon, a, d;
  NumericMatrix D(n,n);
  for(k=0; k < from.length(); k++){
    i = from[k]-1;
    for(l=0; l < to.length(); l++) {
      j = to[l]-1;
      dlat = fabs(lat[i]-lat[j]);
      dlon = fabs(lon[i]-lon[j]);
      dlon = min(dlon, 2*PI-dlon); 
      a = pow(sin(dlat/2),2) + cos(lat[i])*cos(lat[j])*pow(sin(dlon/2),2);
      d = 2*asin(sqrt(a));
      //printf("[%f,%f,%f,%f]", dlat, dlon, a, d);
      D(i,j) = d;
    }
  }
  
  return D;
}