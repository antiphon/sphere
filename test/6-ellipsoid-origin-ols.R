# ellipsoid fitting with center 0
library(devtools)
load_all(".")

library(rgl)

comp <- 0.9
axis <- c(1/sqrt(comp), 1/sqrt(comp), comp)
rot <- rot_ai <- c(pi/5, 0)
x0 <- rellipsoid(100, axis, noise=0.1)
sh <- ellipsoid_shape(2, axis)
plot3d(sh ,aspect=F, lit=T, col=2, alpha=0.2)


##########
el0 <- ellipsoid_OLS_origin(x0)
el <- ellipsoid_OLS(x0)

ci0<-confint(el0)
ci<-confint(el)

print(ci[11:13,])
print(ci0[11:13,])
