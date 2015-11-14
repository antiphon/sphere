#' values to colors
#' 
#' @param v the values
#' @param n number of colors
#' @param zlim limits
#' @param color function, e.g. heat.colors, gray.colors
#' 
#' @export

values2colors <- function(v, n=100, zlim, col=heat.colors, na.col="gray50", ...){
  zlim0 <- range(v, na.rm = TRUE)
  if(missing(zlim)) zlim <- zlim0
  pr <- v
  na <- which(is.na(pr))
  pr[na] <- mean(pr, na.rm=TRUE)
  pr <- pmax(zlim[1], pmin(zlim[2], pr))
  pp <- (pr-zlim[1])/diff(zlim)
  e <- floor(pp*(n-1)) + 1
  hc <-col(n)
  out <- hc[e]
  out[na] <- na.col
  out
}

