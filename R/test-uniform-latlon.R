#' Test of uniformity on a sphere
#' 
#' Independent test of uniformity in azimuth and inclination.
#' 
#' @param ai azimuth-inclination matrix (n x 2)
#' @details
#' Returns the p value of Chi2 test for both angles
#' @export

#### test uniformity on a sphere
test.unifsphere <- function(ai){
  #' we must back transform so that we can test uniformity on a unit interval
  u <- ai[,1]/(2*pi)
  v <- (cos(ai[,2])+1)/2
  #' then we test their uniformity
  bin_u <- hist(u, plot=FALSE)$counts
  bin_v <- hist(v, plot=FALSE)$counts
  pu <- test.chi2(bin_u)
  pv <- test.chi2(bin_v)
  c(p.azi=pu, p.inc=pv)
}


#' Chi2 test for counts in equal size bins
#' @export
test.chi2 <- function(x){
  k <- length(x) # bins
  E <- sum(x)/k
  x2 <- sum( (x-E)^2 )/E
  c(p=1-pchisq(x2, k-1))
}



