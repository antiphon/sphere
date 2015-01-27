#' distance on the plane
#' 
#' @export
#' @useDynLib sphere

dist.plane <- function(x, y, from, to, usec=TRUE){
  n <- length(x)
  if(length(y)!=n | n<2) stop("lat lon should be given for n>2 points.") 
  if(missing(from)) from <- 1:n
  if(missing(to)) to <- 1:n
  
  D <- 
    if(usec) c_plane_dist(x, y, from, to)
  else {
    d <- function(ij) {
      sqrt((x[ij[1]]-x[ij[2]])^2+(y[ij[1]-y[ij[2]]])^2)
    }  
    D <- diag(0, n)
    D[from, to] <- sapply(to, function(j) sapply(from, function(i) d(c(i,j)) ))
    D
  }  
  
  D
}

