#' intensity on a sphere
#' 
#' Gaussian kernel
#' 
#' @param latlon observations
#' @param s sd for gaussian kernel
#' @param nlat if we make a grid, use nlat x (2*lat)
#' @param newlatlon if given predict here, else on grid
#' @param plot should we plot (only with grid)
#' @param points add the points to the plot
#' @export

intensity.sphere <- function(latlon, s=0.25, nlat=20, newlatlon, plot=FALSE, points=FALSE, ...) {
  #' gaussian kernel
  kern <- function(r) dnorm(r, 0, s)
  #'
  #' aggregation points
  if(missing(newlatlon)){ # create a grid
    lat <- seq(-pi/2, pi/2, length=nlat)
    lon <- seq(-pi, pi, length=2*nlat)
    grid <- expand.grid(lat, lon)
    isgrid <- TRUE
  }
  else{ # by user
    lat <- newlatlon[,1]
    lon <- newlatlon[,2]
    grid <- newlatlon
    isgrid <- FALSE
  }
  #
  n<-nrow(latlon)
  m <- nrow(grid)
  
  y <- sphere.smooth(latlon=latlon, newlatlon=grid, s=s)
  
  v <- if(isgrid) t(matrix(y, byrow=F, nrow=length(lat))) else y
  
  if(plot & isgrid) {
    rgltools.setup()
    rgltools.originarrows(col="blue")
    rgltools.im2sphere(v, lat, lon, lit=F, ...)
    if(points) points3d(xyz2globe(ll2xyz(latlon)), col="green")
  }
  
  
  list(x=lon, y=lat, z=v, type=ifelse(isgrid, "grid", "user"))
}