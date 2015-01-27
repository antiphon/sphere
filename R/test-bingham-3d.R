#' Test of isotropy of antipodally symmetric data
#'
#' @param x unit vectors, [x1,...,xn] 
#'
#' @export
bingham.test3d <- function(x) {
  if(nrow(x)!=3) stop("x should be a 3xn matrix giving the unit directions.")
  #' making sure units:
  r <- sqrt(colSums(x^2))
  x <- t(t(x)/r)
  #' hat
  xx <- x%*%t(x)
  mom<-svd(xx)
  n <- sum(mom$d)
  xu2 <- sum( (mom$d - n/3)^2 ) * 15/(2*n)
  #' p-value
  df <- 5
  p <- 1-pchisq(xu2, df)
  #' that is all
  list(statistic=xu2, p.value = p, n=n)
}


