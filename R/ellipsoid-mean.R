#############################################################
#' Mean ellipsoid rotation using quaternion average
#' @export
mean_rotation_quat <- function(x, ...){
  # average ellipsoid using average quaternion
  if(!is.list(x)) stop("x needs to be a list of ellipsoids.")
  if(!all(sapply(x, class)=="ellipsoid")) stop("x needs to be a list of ellipsoids.")
  
  # determine suitable ellipsoids: Valid and Proper
  ok <- which(sapply(x, getElement, "valid") &  sapply(x, function(v) det(v$rot))>0    )
  if(!length(ok)) stop("No suitable rotations for quaternions (either too little data or improper rotations)")
  x <- x[ok]
  
  d <- x[[1]]$dim
  
  quats <- sapply(lapply(x, getElement, "rot"), rotationMatrix2quaternion)
  
  # average
  m_quat <- c( rowMeans(quats, na.rm=TRUE)) 
  m_quat <- m_quat/sqrt(sum(m_quat^2))
  rot <- quaternion2rotationMatrix(m_quat)
  rot
}

#############################################################
#' Mean ellipsoid using empirical mean 
#' 
#' The sample size is proportional to estimated s2
#' 
#' @param x a list of ellipsoids
#' @param nsim Total number of points to simulate from ellipsoids for average estimation
#' @param weight Weigth the ellipsoid samples. Default is to use standard deviation of the OLS estimates.
#' @param just_rotation Rescale each ellipse to have same determinant?
#' @param add_noise Simulate noise using the s2 ols estimate?
#' 
#' @export
mean_ellipsoids <- function(x, nsim=1000, weight, just_rotation=FALSE, add_noise=FALSE, 
                            keep_data=FALSE, ...){
  # average ellipsoid using sampling
  
  if(!is.list(x)) stop("x needs to be a list of ellipsoids.")
  if(!all(sapply(x, class)=="ellipsoid")) stop("x needs to be a list of ellipsoids.")
  
  d <- x[[1]]$dim
  # determine suitable ellipsoids
  ok <- which(sapply(x, getElement, "valid"))
  x_ok <- x[ok]
  # determine sample size distribution using standard deviation
  if(missing(weight)){
    s2 <- sapply(x_ok, function(e) e$ols_fit$s2)
    nv <- sapply(x_ok, function(e) e$ndata  )
    if(!all(s2>0)){ 
      warning("Can't extract s2 estimates, reverting to equal weighting.")
      weight <- rep(1, length(ok))
    }
    else weight <- sqrt(nv/s2)
  }
  else{
    weight <- rep(1, length(ok))
  }
  # the sampling weights
  weight <- weight/sum(weight)
  nsimulations <- round(nsim * weight)
  #
  # check if we need generation
  generate <- sapply(x_ok, function(e) is.null(e$data))
  # generate "data" for those who don't have it
  if(just_rotation){
    # scale each ellipsoid so that determinant = 1
    dets<-sapply(x_ok, function(e) prod(e$semi_axes))
    sca <- 1/dets^(1/d)
    #x_ok <- lapply(x_ok, function(e) {e$semi_axes <- e$semi_axes/prod(e$semi_axes);e})
    for(i in 1:length(x_ok)) x_ok[[i]]$semi_axes <- x_ok[[i]]$semi_axes * sca[i]
    # rescale data in case given
    for(g in which(!generate)){
      x_ok[[g]]$data <- x_ok[[g]]$data * sca[g]
    }
  }
  dats <- if(any(!generate)) do.call(rbind, lapply(x_ok[which(!generate)], getElement, "data")) else NULL
  if(add_noise)
      for(i in which(generate)) 
        dats<-rbind(dats, with(x_ok[[i]], rellipsoid(n=nsimulations[i], axes=semi_axes, R=rot, noise = sqrt(ols_fit$s2))))
    else 
      for(i in which(generate))
        dats<-rbind(dats, with(x_ok[[i]], rellipsoid(n=nsimulations[i], axes=semi_axes, R=rot)))
  # fit
  ave <- ellipsoid_OLS(dats, ...)
  ave$ndata <- nrow(dats)
  ave$ave <- 1
  ave$nellipses <- length(ok)
  if(keep_data) ave$data <- dats
  ave$nsims <- nsimulations
  ave
}



#############################################################
#' Mean of ellipses
#' 
#' @param x a list of ellipses (2d)
#' 
#' @export
mean_ellipse <- function(x, weight, ...){
  if(!is.list(x)) stop("x needs to be a list of ellipses.")
  if(!all(sapply(x, class)=="ellipsoid")) stop("x needs to be a list of ellipses.")
  
  d <- x[[1]]$dim
  if(d!=2) stop("x needs to be a list of ellipses (2d).")
  # determine suitable ellipses
  ok <- which(sapply(x, getElement, "valid"))
  x_ok <- x[ok]
  # determine weights
  if(missing(weight)){
    s2 <- sapply(x_ok, function(v)v$ols_fit$s2)
    n <- sapply(x_ok, getElement, "ndata")
    weight <- n/sqrt(s2)
  }
  weight <- weight/sum(weight)
  if(length(weight)!=length(ok))stop("Weights are not right.")
  #### angles
  ang <- sapply(x_ok, getElement, "rot_angle")
  # map to [-pi/2, pi/2] 
  i <- ang < -pi/2
  ang[i] <- ang[i] + pi
  i <- ang > pi/2
  ang[i] <- ang[i] - pi
  ### semi-axes
  semis <- t( sapply(x_ok, getElement, "semi_axes") )
  
  # 
  ave <- list(valid=FALSE, dim=d)
  R<-function(a)cbind(c(cos(a),sin(a)),c(-sin(a),cos(a)))
  # Averages:
  ave$rot_angle <- sum(ang*weight)
  ave$rot <- R(ave$rot_angle)
  ave$semi_axes <- colSums(semis*weight)
  ave$M <- ave$rot %*% diag(ave$semi_axes)
  ave$A <- ave$rot %*% diag(1/ave$semi_axes^2) %*% t(ave$rot)
  # additional info
  ave$ave <- 1
  ave$nellipses <- length(ok)
  class(ave) <- "ellipsoid"
  ave
}
