# CI ellipsiestimaatille
library(devtools)
load_all(".")

library(rgl)
library(mvtnorm)

comp <- .7
axis <- c(1, comp)
x0 <- rellipsoid(50, axis, noise = 0.5 )
##########
el <- ellipsoid_OLS(x0)

print( sqrt(el$ols_fit$s2 ))

r <- sqrt(rowSums(x0^2))
u <- x0/r
r_pred <- predict(el, u)


plot(x0, asp=1)
plot(el, lwd=3, col=3)



# loop some time
#  rr <- list()
# for(s2 in c(0.05, 0.1)){
#   rr[[paste0("s",s2)]] <- sapply(1:30, function(v) ellipsoid_OLS(rellipsoid(500, c(axis,1), noise=sqrt(s2)))$ols_fit$s2)
# }
# print(sapply(rr, summary))


# Confidence intervals using simulation:

# compression estime
# f <- function(els){
#   R <- lapply(lapply(els, getElement, "A"), ellipse_solve_rota)
#   axes <- t(sapply(R, function(b) b$axes ))
#   angs <- sapply(lapply(R, getElement, "R"), function(R) {f <- R %*% c(1,0)
#                                                           atan2(f[2],f[1])})
#   ab <- sapply(1:length(angs), function(i) if(angs[i] < 3*pi/4 & angs[i] > pi/4) axes[i,1:2] else axes[i,2:1]  )
#   ab[1,]-ab[2,]
# }
# 
 
system.time(print(confint(el, maxiter=0, ellipse_contrast_2d)))
system.time(print(confint(el, maxiter=500, tol=1e-3, ellipse_contrast_2d)))
