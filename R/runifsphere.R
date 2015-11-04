#' Sample unit sphere uniformly
#' 
#' @param n sample size
#' @param spherical Return azimuth-inclination? (Def: FALSE)
#' @details
#' Default is to return unit vectors.
#' @export

runifsphere <- function(n, spherical=FALSE){
  u <- runif(n)
  v <- runif(n)
  a <- u * 2 * pi
  i <- acos(2*v - 1)
  ai<- cbind(azi=a, inc=i)
  if(spherical) ai
  else ai2xyz(ai)
}


#' Sample unit circle uniformly
#' 
#' @param n sample size
#' @param polar Return polar coordinates? (Def: FALSE)
#' @details
#' Default is to return unit vectors.
#' @export

runifcircle <- function(n, polar=FALSE){
  a <- runif(n,  0, 2 * pi)
  if(polar) cbind(r=1, phi=a)
  else cbind(cos(a),sin(a))
}


