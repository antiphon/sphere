#' summary
library(devtools)
load_all(".")

comp <- 0.6
#axis <- c(1/sqrt(comp), 1/sqrt(comp), comp)
#rot <- rot_ai <- c(pi/5, 0)
axis <- c(1/comp, comp)
an <- runif(M<-30, 0, 2*pi)
x0 <- rnorm(M, 1, 0.1)*cbind(cos(an),sin(an))

#plot(x0, asp=1)
##########
el <- ellipsoid_OLS(x0)
summary(el)
