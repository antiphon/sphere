#' test smooth spehre
#' 
library(devtools)
load_all(".")
library(rgl)
set.seed(123)
latlon <- sphere.grid(N=5, lower=T)
n<-20
latlon <- cbind( runif(n, -pi/2, pi/2), runif(n, -pi, pi)   )
v <- rnorm(nrow(latlon), latlon[,1], .1)
ico <- sphere.smooth(latlon, v, N=3, s=.125)
plot(ico, data=T, lit=FALSE)


