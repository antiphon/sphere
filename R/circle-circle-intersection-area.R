#' Circle - circle intersection area
#' 
#' @param d Distance between the circle centers
#' @param r Radius of the first circle
#' @param R Radius of the second circle
#' 
#' @details
#' 
#' @export

circlecirclearea <- function(d, r, R){
  if(d > r + R) 0
  else if(d == 0) min(c(r,R))^2*pi 
  else
  r^2 * acos( (d^2 + r^2 - R^2)/(2*d*r) ) + R^2 * acos((d^2+R^2-r^2)/(2*d*R)) - 
    0.5 * sqrt( (-d +r+R)*(d+r-R)*(d-r+R)*(d+r+R)   )
}