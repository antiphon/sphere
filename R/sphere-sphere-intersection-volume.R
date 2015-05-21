#' Sphere - sphere intersection volume
#' 
#' @param d Distance between the sphere centers
#' @param r Radius of the first sphere
#' @param R Radius of the second sphere
#' 
#' @details
#' 
#' @export

spherespherearea <- function(d, r, R){
  if(d > r + R) 0
  else if(d == 0) min(c(r,R))^3*pi*4/3 
    else
      pi*(R+r-d)^2 * (d^2 + 2*d*r - 3*r^2 + 2*d*R + 6*r*R - 3*R^2)/(12*d)
}