#' generate grids on sphere
#' 
#' @param N number of refinenemts 
#' @param lower return also the lower hemishpere?
#' @param ico use triangulation based on icosahedron (req. rgl)? N is then takes as refinement steps, should be <5!
#' 
#' #' Two ways:
#' 
#'   A) a decreasing radius ring of N points from equator to pole with alternating 2pi/N shift. 
#'   Concentrates on the poles, useful for statitical estimation if the pole area is where things happen.
#'   
#'   B) Triangulation by refining an icosahedron.
#'   
#' Returns latitude-longitude matrix. use ll2xyz or ll2ai for transformation.
#' @import rgl
#' @export

sphere.grid  <- function(N=2, lower=TRUE, ico=TRUE) {
  if(!ico){
    a1 <- seq(0, pi/2, length=N)
    ll <- NULL
    i <- -1
    for(a in a1[-N]){
      m <- N  
      i<-i+1
      d <- 2*pi/m
      b <- seq((i%%2)* d/2, 2*pi-(1-i%%2) * d/2, by=d) 
      ll <- rbind(ll, cbind(lat=a, lon=b-pi))
    }
    coords <- rbind(ll, cbind(pi/2, 0))
    if(lower){
      l <- coords[coords[,1]!=0, ]
      l[,1] <- -l[,1]
      coords <- rbind( coords, l)
    }
  }
  else{
    s0 <- icosahedron3d()
    for(i in 1:N) s0 <- subdivision3d(s0)
    vb <- t(s0$vb[1:3,])/s0$vb[4,]
    coords <- xyz2ll(vb)
    if(!lower) coords <- coords[coords[,1]>=0, ]
  }
  coords
}

#' spherical unit vectors in xyz
#' @export
unit.vecs <- function(...) {
  e <- sphere.grid(...)
  ll2xyz(e)
}
