#' smooth using gaussian regression on a 2d plane
#' 
#' if no v is given, compute intensity.
#' 
#' @useDynLib sphere
#' @export
plane.smooth <- function(xy, v, newxy, s=0.25) {
  int <- missing(v)
  a<-dnorm(0,0,s)
  kern <- function(r) dnorm(r, 0, s)/a
  m <- nrow(newxy)
  n <- nrow(xy)
  all <- rbind(as.matrix(newxy), xy)
  from <- if(int) 1:m else 1:(n+m)
  to <- if(int) 1:n+m else 1:(n+m)
  
  D <- dist.plane(all[,1], all[,2], from=from, to=to)
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
