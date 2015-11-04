#' Smoothed unit sphere values
#' 
#' Krige on a sphere, return a plottable orb.
#' 
#' @param latlon latitude-longitude of data
#' @param v values on data points
#' @param N refinenement of the icosahedron triangulation
#' @param s smoothing sd for gaussian smoothing
#' 
#' @import rgl
#' @export 

sphere.smooth <- function(latlon, v, N=3, s=.25){
  #' the smoothing locations
  ico <- icosahedron3d()
  for(i in 1:N) ico <- subdivision3d(ico)
  newxyz <- t(ico$vb[1:3,])/ico$vb[4,]
  newlatlon <- xyz2ll(newxyz)
  #' predict
  f <- sphere.predict(latlon, v, newlatlon, s=s)
  # assign to icosahedron vertices
  ico$vb <- t( t(ico$vb)/ico$vb[4,])
  ico$data <- if(missing(v)) latlon else cbind(latlon, v)
  ico$prediction <- f
  class(ico) <- c("spheresmooth", class(ico))
  ico
}


#' plot the 3d sphere of predictions
#' @exportMethod plot
#' @import rgl
#' @export
plot.spheresmooth <- function(obj, col=heat.colors, col_zlim, ..., data=FALSE) {
  co <- values2colors(obj$prediction, col=col, zlim=col_zlim, ...)
  cols <- co[c(obj$it)]
  plot3d(obj, col=cols, aspect=FALSE, ...)
  if(data){
    cod <- values2colors(obj$data[,3], col=col, ...)
    text3d(ll2xyz(obj$data[,1:2]), col=cod, texts=1:nrow(obj$data))
  }
}




