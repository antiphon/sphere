#' summary
library(devtools)
load_all(".")

comp <- 0.6
axis <- c(1/comp, comp)
an <- runif(M<-30, 0, 2*pi)
x0 <- rnorm(M, 1, 0.1)*cbind(cos(an),sin(an))
x0 <- t(c(1/comp,comp)*t(x0))
xxx <- {a <- seq(0,2*pi, length=200); cbind(cos(a),sin(a))}
##########
el <- ellipsoid_OLS(x0)

### Dim 2:
f<-el$rot %*% c(1,0)
a<-ang <- atan2(f[2],f[1])

R <- matrix(c(cos(a), sin(a), -sin(a),cos(a)),2)
vm <- t(el$semi_axis*t(xxx))
fe <- vm%*%R

plot(x0, asp=1, cex=.1)
points(vm, pch=19, col=1:10)
points(fe, pch=1, col=1:10)
### OK

# Dim 3:
axis <- c(1/sqrt(c(comp,comp*2)), comp)
rot_euler <- c(0.5, pi/4, 0.4)
x0 <- rellipsoid(1000, axes = axis, R=EulerAngles2rotationMatrix(rot_euler))

el <- ellipsoid_OLS(x0)

print(data.frame(est=el$rot_angle, true=rot_euler))

library(rgl)
sh <- ellipsoid_shape(2, axes=el$semi_axis, R=el$rot)
sh2 <- ellipsoid_shape(2, axes=el$semi_axis, R=EulerAngles2rotationMatrix(el$rot_angle))

plot3d(x0, aspect=FALSE)
shade3d(sh, col=3, lit=F, alpha=0.4)
shade3d(sh2, col=4, lit=F, alpha=0.4)
