#' Uniform sample on the surface of an ellipsoid
#'
#' @param R rotation matrix
#' @export
rellipsoid <- function(n, axes=c(1,1,1), noise=0, R=NULL){
  # transform points on the surface of an circle
  # points on a sphere/circle
  d <- length(axes)
  x <- if(d==3) runifsphere(n) else runifcircle(n)
  # transform
  if(is.null(R)) R <- diag(1, d)
  M <- diag(axes)
  y <- t(R%*%M%*%t(x))
  
  
  if(noise) y <- y + rnorm(nrow(y)*d, 0, noise)
  y
}

###########################
#############################################################
#' OLS fitting of an ellipsoid to d-dimensional point cloud
#'
#' @param x point coordinate matrix
#' 
#' @export
ellipsoid_OLS  <- function(x) {
  n <- nrow(x)
  d <- ncol(x)
  if(n < 3) stop("Trying to fit an ellipsoid to less than 3 points.")
  D  <- matrix(2, ncol=d, nrow=d)
  diag(D) <- 1
  D <- D[upper.tri(D,T)]
  ## Y:
  Y<-apply(x, 1, function(x){
    X  <- x%*%t(x)
    vX <- X[upper.tri(X, T)]
    y <- c(D*vX, c(x), 1)
    y
  })
  Y <- t(Y)
  #
  nd <- (d*(d+1)/2)
  nb <- nd + d + 1 
  H <- diag(1, nb)
  Hi <- solve(H)
  # svd
  USV <- svd(Y%*%Hi)
  beta <- Hi%*%USV$v[,nb]
  
  # ellipse parameters
  A <- diag(0, d)
  A[upper.tri(A,T)] <- beta[1:nd]  
  A[lower.tri(A)] <- A[upper.tri(A)]
  b <- beta[(nd+1):(nd+d)]
  dhat <- beta[nb]
  
  chat <- -0.5 * solve(A)%*%b
  Ahat <- A / (t(chat) %*% A %*% chat - dhat)[1]
  
  # rotation and axes:
  ev <- eigen(Ahat)
  axes_len <- 1/sqrt(ev$value) # the semi-axes lengths
  R <- ev$vector # this holds the rotation...
  M <- R%*%diag(axes_len)
  angles <- NULL
  if(d==2) {
    f <- R %*% c(1,0)
    angles <- atan2(f[2],f[1])
  }else if(d==3){
    angles <- rotationMatrix2EulerAngles(R)
  }
  
  # check if we got a valid fit
  valid <- all(!is.na(c(axes_len, angles)))
  # compile
  
  res <- list(center=chat, A=Ahat, semi_axes=axes_len, rot=R, M=M, dim=d, ndata=n, 
              rot_angle=angles, valid=valid)
  class(res) <- "ellipsoid"
  res
}

#############################################################
#' Summarise an ellipsoid
#' 
#' @exportMethod summary
#' @export
summary.ellipsoid <- function(x, ...){
  print(x)
  angtxt <- ifelse(x$dim==2, format(x$rot_angle),paste(c("azimuth:", "inclination:"), format(x$rot_angle), collapse=" "))
  
  
  
  cat("\nEstimates:\n Center:\t \t", paste0("(", paste0(format(x$center), collapse=", "), ")\n"))
  cat(" Semi-axes lengths:\t ", paste0(format(x$semi_axes), collapse=" : "), "\n")
  
  cat(" Rotation angles (rad):\t ", angtxt)
}





############################################################
#' Plot an ellipsoid
#' 
#' @exportMethod plot
#' @export
plot.ellipsoid <- function(x, add=TRUE, persp3d=FALSE, ...){
  if(x$dim==2){
    a <- c(seq(0, 2*pi, length=100))
    y <- cbind(cos(a),sin(a))
    z <- t(x$M%*%t(y))
    lines(z, ...)
  }
  else{
    if(persp3d){
      
    }
    else{
      y <- ellipsoid_shape(axes=x$semi_axes, R=x$rot)
      shade3d(y, ...)
    }
  }
}

####################################################################
#' Ellipsoid shape for 3d plotting
#'
#' Refine an icosahedron, then transform
#' 
#' @param R Rotation matrix.
#' 
#' @export
ellipsoid_shape <- function(N=2, axes=c(1,1,1), R=NULL){
  ico<-icosahedron3d()
  for(i in 1:N) ico <- subdivision3d(ico)
  D <- diag(axes)
  xy <- t(D%*%ico$vb[1:3,])
  exy <- if(is.null(R)) xy else t(R%*%t(xy))
  ico$vb[1:3,] <- t(exy)
  ico$vb <- t( t(ico$vb)/apply(ico$vb, 2, function(v) sqrt(sum(v^2))))
  ico
}

############################################################
#' Print ellipsoid
#' 
#' @exportMethod print
#' @export
print.ellipsoid <- function(x, ...){
  type <- ifelse(x$dim==2, "2D ellipse", "3D ellipsoid")
  if(!is.null(x$ave))
    cat(paste0("Average ", type, ", computed from ", x$nellipses, " ", type, "s.\n"))
  else cat(type, "fitted to", x$n, "points.\n")
}


#############################################################
#' Mean ellipsoid using quaternion average
#' @export
mean_ellipsoids_quat <- function(x, ...){
  # average ellipsoid using average quaternion
  
  if(!is.list(x)) stop("x needs to be a list of ellipsoids.")
  if(!all(sapply(x, class)=="ellipsoid")) stop("x needs to be a list of ellipsoids.")
  
  d <- x[[1]]$dim
  
  semis <- sapply(x, getElement, "semi_axes")
  quats <- sapply(lapply(x, getElement, "rot"), rotationMatrix2quaternion)
  cents <- sapply(x, getElement, "center")
  
  # average
  semi <- rowMeans(semis, na.rm=TRUE)
  rot <- quaternion2rotationMatrix(c( rowMeans(quats, na.rm=TRUE)) )
  M <- rot%*%diag(semi)
  center <- rowMeans(cents, na.rm=TRUE)
  angles <- NULL
  if(d==2) {
    f <- rot %*% c(1,0)
    angles <- atan2(f[2],f[1])
  }else if(d==3){
    angles <- rotationMatrix2EulerAngles(rot)
  }
  
  ave <- list(M=M, rot=rot, semi_axes=semi, center=center, ndata=NA, dim=d, 
              rot_angle=angles, ave=1, nellipses=length(x))
  class(ave) <- "ellipsoid"
  ave
}

#############################################################
#' Mean ellipsoid using empirical mean 
#' @export
mean_ellipsoids <- function(x, nsim=100, ...){
  # average ellipsoid using sampling
  
  if(!is.list(x)) stop("x needs to be a list of ellipsoids.")
  if(!all(sapply(x, class)=="ellipsoid")) stop("x needs to be a list of ellipsoids.")
  
  d <- x[[1]]$dim
  # generate "data"
  ok <- sapply(x, getElement, "valid")
  dats <- lapply(x[ok], function(e) with(e, rellipsoid(n=nsim , axes=semi_axes, R=rot)))
  dats <- do.call(rbind, dats)
  # fit
  ave <- ellipsoid_OLS(dats)
  ave$ndata <- nrow(dats)
  ave$ave <- 1
  ave$nellipses <- sum(ok)
  ave$nsim <- nsim
  ave
}
