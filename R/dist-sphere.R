#' distance on the surface of a sphere
#' 
#' @useDynLib sphere
#' @export

dist.sphere <- function(lat, lon, R=1, from, to, usec=TRUE){
  n <- length(lat)
  if(length(lon)!=n | n<2) stop("lat lon should be given for n>2 points.") 
  if(missing(from)) from <- 1:n
  if(missing(to)) to <- 1:n
  
  D <- 
    if(usec) c_sphere_dist(lat, lon, from, to)
    else {
      d <- function(ij) {
        i <- ij[1]
        j <- ij[2]
        dlat <- abs(lat[i]-lat[j])
        dlon <- abs(lon[i]-lon[j])
        dlon <- min(dlon, 2*pi-dlon) 
        a <- sin(dlat/2)^2 + cos(lat[i])*cos(lat[j])*sin(dlon/2)^2
        2*asin(sqrt(a))
      }  
      D <- diag(0, n)
      D[from, to] <- sapply(to, function(j) sapply(from, function(i) d(c(i,j)) ))
      D
    }  

  R*D
}

