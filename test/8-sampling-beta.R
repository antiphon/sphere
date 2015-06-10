# sampling from the betas
library(devtools)
library(mvtnorm)
load_all(".")


v <- rellipsoid(axes = c(1,1.01,2), 100, noise=.1)
x <- ellipsoid_OLS(v)

b<-sample_ellipse_beta(x, 1000, tol = 1e-4)

par(mfrow=c(2,4))
apply(b, 2, hist)
hist(rowSums(b^2)-1)


plot(v, asp=1)
els <- apply(b, 1, ellipsoid_from_beta, d=2)

sapply(els, plot, col="gray90")


ci<-confint(x, nsim=2500, ellipse_contrast_2d, tol=1e-4)
cio <- confint(x, nsim=2500, function(...)ellipse_contrast_2d(..., out=TRUE), tol=1e-4)
