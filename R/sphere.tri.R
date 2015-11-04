#' Triangulation of a sphere
#'
#' For now, retriangulate icoshedron3d() from rgl.
#'  
#' @import rgl
#' @export

sphere.tri <- function(N=2) {
  g <- icosahedron3d()
  for(i in 1:N) g <- subdivision3d(g)
  g
}