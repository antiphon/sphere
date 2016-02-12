#' Predict on sphere using gaussian regression (Kriging)
#' 
#' if no v is given, compute intensity.
#' @param latlon latitude-longitude matrix
#' @param v value at locations. if not given, we do intensity estimation.
#' @param newlatlon new locations
#' @param s gaussian kernel sd
#' 
#' use xyz2ll or ai2ll to transform from unitvecs or (azimuth, inclination) (physical)
#' 
#' @export
sphere.predict <- function(latlon, v, newlatlon, s=0.25) {
  int <- missing(v)
  a<-dnorm(0,0,s)
  kern <- function(r) dnorm(r, 0, s)/a
  m <- nrow(newlatlon)
  n <- nrow(latlon)
  all <- rbind(as.matrix(newlatlon), latlon)
  from <- if(int) 1:m else 1:(n+m)
  to <- if(int) 1:n+m else 1:(n+m)
  
  D <- dist.sphere(all[,1], all[,2], from=from, to=to)
  K <- kern(D[1:m, 1:n+m])
  
  # intensity or prediction
  hedge <- 0.0001
  if(!int){
    Kx <- kern(D[1:n+m, 1:n+m]) 
    while("try-error" %in% is( Ki<-try(solve(Kx), silent=TRUE) ) ){
      Kx <- Kx + diag(hedge, n)
      cat("hedging.\n")
    }
    y <- K%*%Ki%*%v
  }
  else y <- rowSums(K) 
  c(y)
}
