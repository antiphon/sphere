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
#' @param origin Fix center of the ellipsoid to origin?
#' 
#' @export
ellipsoid_OLS  <- function(x, origin=FALSE) {
  n <- nrow(x)
  d <- ncol(x)
  if(n < 3+d) stop("Trying to fit an ellipsoid to too little amount of points.")
  D  <- matrix(2, ncol=d, nrow=d)
  diag(D) <- 1
  D <- D[upper.tri(D,T)]
  ## Y:
  Y<-apply(x, 1, function(x){
    X  <- x%*%t(x)
    vX <- X[upper.tri(X, T)]
    y <- c(D*vX, if(origin) NULL else c(x), 1)
    y
  })
  Y <- t(Y)
  #
  nd <- (d*(d+1)/2)
  nb <- nd + (d * !origin) + 1 
  H <- diag(1, nb)
  #if(d==2) H[2,2]<-sqrt(2) else H[cbind(c(2,4,5),c(2,4,5))] <- sqrt(2)
  
  Hi <- solve(H)
  # svd
  USV <- svd(Y%*%Hi)
  # the estimate of parameters:
  beta <- c( Hi%*%USV$v[,nb] )
  if(origin){ # add center parameters to comply with other functions
    beta <- c(beta[1:nd], rep(0, d), beta[(nd+1):nb])
  }
  
  # convert to center and matrix -formulation
  res <- ellipsoid_from_beta(beta, d)
  res$ndata <- n
  ####### add the error variance
  # the error variance
  r_data <- sqrt(rowSums(x^2))
  u_data <- x/r_data
  r_pred <- predict(res, u_data)
  resi <- r_data - r_pred
  if(d==2){
    s2 <- sum(resi^2)/(n-1) # this is approximation with 0 curvature TODO better
  }
  else{
    s2 <- sum(resi^2)/(n-1) # this is approximation to 0 curvature TODO BETTER
  }
  # the variance covariance
  warn <- FALSE; m <- 1e-9 # hedging
  ss <- try(S0 <- solve(t(Y)%*%Y + diag(m, ncol(Y))), TRUE)
  while("try-error"%in% is(ss)){
    ss <- try(S0 <- solve(t(Y)%*%Y + diag(m <- m * 2, ncol(Y))), TRUE)
    warn <- TRUE
  }
  if(warn) warning(paste("Hedging needed for solving covariance matrix. (diag ", m, ")"))
  S <- s2 * S0
  #
  if(origin){ # add dummy for center
    S0 <- diag(1e-9, nb+d)
    S0[1:(nb-1), 1:(nb-1)] <- S[-nb,-nb]
    aa <- c(1:(nb-1), nb+d)
    S0[aa, nb+d] <- S0[nb+d, aa] <- S[nb,]
    S0 <- 0.5 * (S0 + t(S0))
    S <- S0
  }
  # make sure symmetric varcov

  if(max(S-t(S))>0.01) warning("varcov not good, big asymmetry.")
  S <- 0.5 * (S + t(S))
  #
  res$ols_fit <- list(varcov=S, beta_est=beta, s2=s2)
  res$origin <- origin
  #
  #' done
  res
}

#############################################################
#' Summarise an ellipsoid
#' 
#' @exportMethod summary
#' @export
summary.ellipsoid <- function(x, ...){
  print(x)
  ang <- round(x$rot_angle, 3)
  angd <- round(x$rot_angle * 180/pi, 1)
  semi_axes <- x$semi_axes
  semi_axes_rel <- round(semi_axes/semi_axes[1], 3)
  angtxt <- ifelse(x$dim==2, format(ang), paste(c("heading", "attitude", "bank"), format(ang), collapse=" "))
  angtxtd <- ifelse(x$dim==2, format(angd), paste(c("heading", "attitude", "bank"), format(angd), collapse=" "))
  cat("\nEstimates:\n Center:\t \t ", paste0("(", paste0(format(x$center), collapse=", "), ")\n"))
  cat(" Semi-axes lengths (absolute):\t ", paste0(format(semi_axes), collapse=" : "), "\n")
  cat(" Semi-axes lengths (relative):\t ", paste0(format(semi_axes_rel), collapse=" : "), "\n")
  
  cat(" Rotation angles (rad):\t ", angtxt,"\n")
  cat(" Rotation angles (deg):\t ", angtxtd,"\n")
  cat(" Error variance: \t ", x$ols_fit$s2, "\n")
  invisible(NULL)
}

#################################################################
#' Confidence interval for ellipsoid parameters using simulation
#' 
#' Assuming normality of errors.
#' 
#' @param tol tolerance for absolute  deviation in the ||beta||=1 constraint.
#' 
#' @exportMethod confint
#' @import mvtnorm
#' @export
confint.ellipsoid <- function(x, fun, nsim=1000, probs=c(0.025, 0.975),  ...){
  S <- x$ols_fit$varcov
  m <- x$ols_fit$beta_est
  f <- qt(0.975, x$ndata-2)
  s <- sqrt(diag(S))
  L <- m - s*f
  U <- m + s*f
  names(m) <- paste0("beta", 1:length(m))
  pval <- round(2*pt(-abs(m/s) , x$ndata - 2) , 5)
  df_orig <- data.frame(mean=m, median=m, sd=s, L,U, p=pval)
  # simulate a bunch for semi_axes and rotations
  
  betas <- sample_ellipse_beta(x, nsim, ...)
  
  
  efs <- apply(betas, 1, function(b) ellipse_form(b, x$dim) )
  rots <- lapply(efs, function(e) ellipse_solve_rota(e$A) )
  axes <- t(sapply(rots, function(b) b$axes ))
  #
  # summary for one vector
  e <- function(v,h0=0) {
    # drop bad values
    v <- v[!is.infinite(v) & !is.na(v)]
    z <- (h0-mean(v))/sd(v)
    p <- round(2*pt(-abs(z), x$ndata-2), 5)
    c(mean=mean(v), median=median(v), sd=sd(v), CI=quantile(v, probs=probs[1]), CI=quantile(v,probs=probs[2]), p=p)
  }
  # basics for axes
  df_axes <- t(apply(axes, 2, e))
  #
  df_fun <- df_ex <- df_diff <- NULL
  if(x$dim==2){
    # eccentricity
    df_ex <- e(sqrt(1-(axes[,1]/axes[,2])^2 )) # eccentricity
    df_ex[6] <- NA # these p-values are not correct
    # custom function
    if(missing(fun))
      df_diff <- e(ellipse_contrast_2d(efs))
  }
  
  if(!missing(fun)) df_fun <- e(fun(efs))
  
  names(df_orig) <- colnames(df_axes)
  rbind(beta=df_orig, axes=df_axes, eccentricity=df_ex, isotropy=df_diff, custom=df_fun)
}


############################################################
#' Plot an ellipsoid
#' 
#' @exportMethod plot
#' @export
plot.ellipsoid <- function(x, add=TRUE, persp3d=FALSE, res=201, scale=1, ...){
  if(x$dim==2){
    a <- c(seq(0, 2*pi, length=res))
    y <- cbind(cos(a),sin(a))
    z <- y * predict(x, y) * scale
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
#' Predict i.e. give the length of a direction to be on the ellipsoid
#'
#' Returns the distance from ellipsoid center to the ellipsoid surface in
#' the given directions.
#'
#'@exportMethod predict
#'@export
predict.ellipsoid <- function(x, u, ...){
  if(missing(u)) stop("direction(s) u needed")
  d <- 1 / diag(u %*% x$A %*% t(u) )
  r <- sqrt(d) 
  r
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

####################################################################
#' Ellipse center and matrix from general parameter form
#' 
#' @param beta OLS estimates
#' @param d dimension
#' @export
ellipse_form <- function(beta, d){
  nd <- (d*(d+1)/2)
  nb <- nd + d + 1 
  # ellipse parameters
  A <- diag(0, d)
  A[upper.tri(A,T)] <- beta[1:nd]  
  A[lower.tri(A)] <- A[upper.tri(A)]
  b <- beta[(nd+1):(nd+d)]
  dhat <- beta[nb]
  m <- 0
  v <- try(Sa <- solve(A + diag(m, ncol(A))))
  while("try-error"%in% is(v))  v <- try(Sa <- solve(A + diag(m<-m+5e-8, ncol(A))))
  chat <- -0.5 * Sa%*%b
  Ahat <- A / (t(chat) %*% A %*% chat - dhat)[1]
  # check definitenes
  e <- eigen(Ahat)
  if(any(e$values < 0)){
    i <- which(e$values > 0)
    Av <- diag(0, ncol(A))
    for(j in i) Av <- Av + e$val[j] * e$vec[,j]%*%t(e$vec[,j])
    Ahat <- Av
  }
  list(c=chat, A=Ahat)
}


#########################################################################
#' Solve rotation and semi-axes from general transform
#' @param A the trasformation matrix in ellipsoid equation
#' @export

ellipse_solve_rota <- function(A){
  ev <- eigen(A)
  # make sure working with a definite:
  if(any(ev$value<0)){
    i <- which(ev$value > 0)
    S <- diag(0, ncol(A))
    for(j in i) S <- S + ev$value[j] * ev$vec[,i]%*%t(ev$vec[,i])
    ev <- eigen(S)
  }
  axes_len <- 1/sqrt(ev$value) # the semi-axes lengths
  R <- ev$vector # this holds the rotation...
  # check proper
#   if(det(R)<0){
#     # mirror
#     R <- Euler2rotationMatrix(rotationMatrix2Euler(R))
#   }
  list(axes=axes_len, R=R)
}


##########################################################################
#' Sample from the OLS estimate of the beta parameters
#'
#' @import mvtnorm 
#' @export
sample_ellipse_beta <- function(x, nsim=100, tol=0, maxiter=5000){
  d <- x$dim
  nb  <- (d*(d+1)/2) + d + 1
  
  b <- rmvnorm(nsim, x$ols_fit$beta_est, x$ols_fit$varcov)
  if(maxiter==0){ tol<-0; b <- b/sqrt(rowSums(b^2))}
  if(tol>0){
    dev <-  abs(sqrt(rowSums(b^2))-1)
    ok <- dev < tol
    b <- b[ok,]
    it <- 0
    failed <- FALSE
    while(length(b)/nb < nsim){
      b <- rbind(b, rmvnorm( 2*(nsim-length(b)/nb), x$ols_fit$beta_est, x$ols_fit$varcov))
      ok <- abs(sqrt(rowSums(b^2))-1) < tol
      b <- b[ok,]
      it<-it+1
      if(it>maxiter) {it<-0; tol <- tol * 10; failed<-TRUE}
    }
    if(failed) warning(paste("beta sampling tolerance was increased to", tol))
    b<-b[1:nsim,]
  }
  b
}


#####################################################################################
#' convert beta vector to ellipsoid object
#' 
#' @param beta beta vector, the coefficients in quadratic form
#' @param d dimension
#' @export
ellipsoid_from_beta <- function(beta, d){
  elform <- ellipse_form(beta, d)
  chat <- elform$c
  Ahat <- elform$A
  # solve rotation and axes:
  rota <- ellipse_solve_rota(Ahat)
  R <- rota$R
  axes_len <- rota$axes
  M <- R%*%diag(axes_len)
  angles <- NULL
  if(d==2) {
    f <- R %*% c(1,0)
    angles <- atan2(f[2],f[1])
  }else if(d==3){
    angles <- NULL#rotationMatrix2EulerAngles(R)
  }
  
  # check if we got a valid fit
  valid <- all(!is.infinite(c(axes_len, angles)) & !is.na(c(axes_len, angles)))
  # compile
  
  res <- list(center=chat, A=Ahat, 
              semi_axes=axes_len, 
              rot=R, 
              M=M, 
              rot_angle=angles, 
              valid=valid, dim=d)
  class(res) <- "ellipsoid"
  res
}