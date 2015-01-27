# ellipsoid fitting in 3d
library(devtools)
load_all(".")

library(rgl)

comp <- 0.9
axis <- c(1/sqrt(comp), 1/sqrt(comp), comp)
rot <- rot_ai <- c(pi/5, 0)
x0 <- rellipsoid(100, axis, rot_ai = rot)
sh <- ellipsoid_shape(2, axis, rot)
plot3d(sh ,aspect=F, lit=T, col=2, alpha=0.2)


##########
el <- ellipsoid_OLS(x0)
R <- el$rot 
for(i in 1:3){v <- R[,i] * el$semi_axis[i]; lines3d(c(0,v[1]), c(0,v[2]), c(0,v[3]), col=i) }


R%*%diag(1,3)
