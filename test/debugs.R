library(Kdirectional)

library(devtools)

library(rgl)

load_all(".")

#s <- readRDS("../../papers/directed_summaries/code/example_patterns/intensity_250_small_cut_compressed_0.8.rds")
# x<-s[[2]]
# 
# el <- fry_ellipsoids(x, border=TRUE, eps=0.1,r_adjust = 1, nangles=2)$el  
# ave <- mean_ellipsoids( el )

#load("bug.rda")
#e<-ellipsoid_OLS(x)



load("debug_mean_ellipsoid1.rds")
els <- aves1
m <- mean_ellipsoids(els)
