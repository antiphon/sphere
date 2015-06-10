#' Circle - circle intersection area
#' 
#' @param d Distance between the circle centers
#' @param r Radius of the first circle
#' @param R Radius of the second circle
#' 
#' @export

circlecirclearea <- function(d, r, R){
  m <- min(c(r,R))
  M <- max(c(r,R))
  if(d > r + R) 0 # not overlapping
  else if(d + m < M ) m^2*pi  # colocated or the other is completely inside the other
  else 
  r^2 * acos( (d^2 + r^2 - R^2)/(2*d*r) ) + R^2 * acos((d^2+R^2-r^2)/(2*d*R)) - 
    0.5 * sqrt( (-d +r+R)*(d+r-R)*(d-r+R)*(d+r+R)   )
}