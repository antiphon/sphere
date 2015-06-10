#' Sphere - sphere intersection volume
#' 
#' @param d Distance between the sphere centers
#' @param r Radius of the first sphere
#' @param R Radius of the second sphere
#' 
#' @export

spherespherearea <- function(d, r, R){
  m <- min(c(r,R))
  M <- max(c(r,R))
  if(d > r + R) 0
  else if(d + m < R) m^3*pi*4/3  # one engulfs the other
    else
      pi*(R+r-d)^2 * (d^2 + 2*d*r - 3*r^2 + 2*d*R + 6*r*R - 3*R^2)/(12*d)
}