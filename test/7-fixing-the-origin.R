# test origin fixing
# CI ellipsiestimaatille
library(devtools)
load_all(".")

library(rgl)
library(mvtnorm)

comp <- .6
axis <- c(1/comp, comp, 1)
x0 <- rellipsoid(500, axis, noise = 0.1 )
##########
el <- ellipsoid_OLS(x0)
el0 <- ellipsoid_OLS(x0, T)
summary(el)
summary(el0)


print(confint(el))
print(confint(el0))
