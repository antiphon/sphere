#' plotting of 3d ellipses
library(devtools)
load_all(".")
library(rgl)
set.seed(1)
pp <- rellipsoid(1000, c(1,.5,0.4))

x <- e <- ellipsoid_OLS(pp)
#plot3d(x, aspect=F)
#s<-ellipsoid_shape(2, x$semi_axes, x$rot)
#plot3d(s, col=values2colors(s$vb[3,s$it]))
#par(mfrow=c(3,3), mar=c(0,0,0,0))

#for(i in 1:9) persp(e, theta=5+i*15, border=NA)

## overlaying two ellipses
pp2 <- rellipsoid(1000, c(1,.3,0.4))
x2 <- e2 <- ellipsoid_OLS(pp2)

colmap <- function(z,a) values2colors(z, col=function(...) heat.colors(..., alpha=a))

pm <- persp(e, colmap=function(v)colmap(v,1), border=NA)
persp(e2, add=T, pmat=pm, colmap=function(v)colmap(v,0.125), border=NA)

x <- list(x,x2)

m <- mean_ellipsoids(x)

plot(m)
