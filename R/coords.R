#' OLD for plotting on the globe, we need to trasform as the texture is fuckd up
#' @export
xyz2globe <- function(xyz) {
  xyzrotate(xyz, ax=pi/2, ay=-pi/2)
}


#' latlon to 3d coordinates
#' @export
ll2xyz <- function(latlon) {
  azi <- lon2azi(latlon[,2])
  inc <- lat2inc(latlon[,1])
  ai2xyz(cbind(azi,inc))
}

#' (azi,incl) to 3d coordinates
#' @export
ai2xyz <- function(aziinc) {
  azi <- aziinc[,1]
  inc <- aziinc[,2]
  wx<-sin(inc)*cos(azi)
  wy<-sin(inc)*sin(azi)
  wz<-cos(inc)
  cbind(x=wx,y=wy,z=wz)
}

#' latitude to inclination
#' @export
lat2inc <- function(lat) {
  pi/2 - lat
}
#' inclination to latitude
#' @export
inc2lat <- function(inc){
  pi/2 - inc
}
#' azimuth to longitude
#' @export
azi2lon <- function(azi){
  lon <- azi
  lon[lon>pi] <- -(2*pi-azi[azi>pi])
  lon
}
#' longitude to azimuth
#' @export
lon2azi <- function(lon){
  azi <- lon
  azi[azi < 0] <- 2*pi + azi[azi < 0]
  azi
}

#' antipode lat and lon
#' @export
antipode <- function(latlon){
  latlon<-cbind(latlon)
  cbind(-latlon[,1], latlon[,2]-sign(latlon[,2])*pi)
}

#' xyz to lat lon
#' @export
xyz2ll <- function(xyz) {
  ai <- xyz2ai(xyz)
  inc <- ai[,2]
  azi <- ai[,1]
  cbind(lat=inc2lat(inc), lon=azi2lon(azi))
}

#' azi-inc to lat lon
#' @export
ai2ll <- function(aziinc){
  aziinc <- rbind(aziinc)
  cbind(lat=inc2lat(aziinc[,2]), lon=azi2lon(aziinc[,1]))
}

#'latlon to aziinc
#' @export
ll2ai <- function(latlon){
  latlon <- rbind(latlon)
  cbind(azi=lon2azi(latlon[,2]), inc=lat2inc(latlon[,1]))
}


#' xyz to azi inc
#' @export
xyz2ai <- function(xyz) {
  xyz <- rbind(xyz)
  r <- apply(xyz, 1, function(a)sqrt(sum(a^2)))
  inc <- acos(xyz[,3]/r)
  azi <- atan2(xyz[,2],xyz[,1])
  cbind(azi=azi, inc=inc)
}


#' rotate xyz coordinates
#' @export 
xyzrotate <- function(xyz, ...) {
  R <- rotationMatrix(...)
  xyz %*% R
}

#' Product of elementary rotation matrices
#' 
#' @param order Order in which to to the product. Aimed for pre-product
#' Default: 3-2-1 i.e. around z, then around y, then around x. 
#' @export 
rotationMatrix <- function(ax=0, ay=0, az=0, order=c(3,2,1)) {
  R<-list()
  R[[1]] <- cbind(c(1,0,0), c(0, cos(ax), sin(ax)), c(0, -sin(ax), cos(ax)) )
  R[[2]] <- cbind(c(cos(ay), 0, -sin(ay)), c(0, 1, 0), c(sin(ay), 0, cos(ay)) )
  R[[3]] <- cbind(c(cos(az), sin(az), 0), c(-sin(az), cos(az), 0), c(0,0,1) )
  R[[order[3]]]%*%R[[order[2]]]%*%R[[order[1]]]
}

#' Rotation using azimuth-inclination.
#' 
#' @export
aziinc2rotationMatrix <- function(rot_ai=c(0,0)){
  xyzrotationMatrix(az=rot_ai[1])%*%(xyzrotationMatrix(ay=rot_ai[2]))
}


#' Rotation matrix to quaternion
#' 
#' @export
rotationMatrix2quaternion <- function(R){
  q4 <- 0.5 * sqrt(1 + R[1,1] + R[2,2] + R[3,3])
  if(is.na(q4)) q4<-0 # improper rotation
  if(q4 < 1e-5){ # try other way
    q1 <- 0.5 * sqrt(1+R[1,1]-R[2,2]-R[3,3])
    if(q1 < 1e-5){ # try other way
      q2 <- 0.5 * sqrt(1-R[1,1]+R[2,2]-R[3,3])
      if(q2 < 1e-5){
        q3 <- 0.5 * sqrt(1-R[1,1]-R[2,2]+R[3,3])
        if(q3 < 1e-5) stop("Can't convert matrix to quaternion.")
        m <- 1/(4*q3) # q3 ok
        q1 <- m * (R[3,1]+R[1,3])
        q2 <- m * (R[3,2]+R[2,3])
        q4 <- m * (R[2,1]-R[1,2])
      }else{ # q2 ok
        m <- 1/(4*q2)
        q1 <- m * (R[2,1]+R[1,2])
        q3 <- m * (R[3,2]+R[2,3])
        q4 <- m * (R[1,3]-R[3,1])
      }
    }else{ # q1 ok 
      m <- 1/(4*q1)
      q2 <- m * (R[1,2]+R[2,1])
      q3 <- m * (R[1,3]+R[3,1])
      q4 <- m * (R[3,2]-R[2,3])
    }
  }else{ # q4 ok
    m <- 1/(4*q4)
    q1 <- m * (R[3,2]-R[2,3])
    q2 <- m * (R[1,3]-R[3,1])
    q3 <- m * (R[2,1]-R[1,2])
  }
  c(q1,q2,q3,q4)
}

#' Quaternion to rotation matrix
#' 
#' @export
quaternion2rotationMatrix <- function(q){
  q3 <- q[-4]
  Q <- matrix(c(0,q[3],-q[2],-q[3],0,q[1],q[2],-q[1],0), 3)
  c(q[4]^2-t(q3)%*%q3) * diag(1,3) + 2*q3%*%t(q3) + 2*q[4]*Q
}

#' Quaternion to Euler axis/angle
#' @export
quaternion2Euler <- function(q){
  e <- q[-4]/sqrt(sum(q[-4]^2))
  theta <- 2 * acos(q[4])
  c(e=e, theta=theta)
}


#' Euler axis/angle to quaternion
#' @param euler c(e1,e2,e3,angle) where e* is the unit axis of rotation
#' 
#' @export
Euler2quaternion <- function(euler){
  theta <- euler[4]
  m <- sin(theta/2)
  c(q=euler[1:3]*m, q4=cos(theta/2))
}

#' Quaternion to Euler angles
#' 
#' from euclideanspace.com,
#' 
#' @details
#' 
#' Convention: Euler is (Heading, Attitude, Bank) so that
#' Heading = rotation  around y, Attitude = rotation around z, 
#' Bank = rotation around x, in that order.
#' 
#' @export
quaternion2EulerAngles <- function(q){
#   if(q[2]*q[3]+q[3]*q[1]==0.5){
#     phi <- 2*atan2(q[2],q[1])
#     
#   }
  phi <- atan2(2*(q[1]*q[3]-q[2]*q[4]), 1-2*(q[3]^2+q[4]^2))
  theta <- asin(2*(q[2]*q[3]+q[1]*q[4]))
  psi <- atan2(2*(q[1]*q[2]-q[3]*q[4]), 1-2*(q[2]^2+q[4]^2))
  c(phi=phi, theta=theta, psi=psi)
}

#' Euler angles to quaternions
#' @export
EulerAngles2quaternion <- function(angle){
  cc <- cos(angle/2)
  ss <- sin(angle/2)
  q1 <- cc[1]*cc[2]*cc[3] - ss[1]*ss[2]*ss[3]
  q2 <- ss[1]*ss[2]*cc[3] + cc[1]*cc[2]*ss[3]
  q3 <- ss[1]*cc[2]*cc[3] + cc[1]*ss[2]*ss[3]
  q4 <- cc[1]*ss[2]*cc[3] - ss[1]*cc[2]*ss[3]
  c(q1,q2,q3,q4)
}


#' Rotation matrix (prop or improp) to Euler axis-angle
#' 
#' @export
rotationMatrix2Euler <- function(R){ 
  e <- det(R)
  Tr <- sum(diag(R))
  a <- 3-e*Tr
  b <- 1+e*Tr
  if(a==0){ 
    theta <- acos(e)
    n <- c(1,0,0)
  }
  else if(b==0){
    theta <- acos(-e)
    eR <- e*R
    e1 <- e2 <- e3 <- 1
    n1 <- sqrt(1+eR[1,1])
    n2 <- sqrt(1+eR[2,2])
    n3 <- sqrt(1+eR[3,3])
    if(n1){
      if(n2) e2 <- eR[1,2]/(n1*n2)
      if(n3) e3 <- eR[1,3]/(n1*n3)
    }else{
      if(n2 & n3) e2 <- e3 <- eR[2,3]/(n2*n3)
    }
    n <- c(n1,n2,n3)*c(e1,e2,e3)
    n <- n/sum(n)
  }
  else{
    theta <- acos(0.5 * (Tr-e))
    n <-  c(R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2]) /sqrt(a*b)
  }
  v <- c(axis=n, angle=theta)
  if(e < 0) attr(v, "proper") <- FALSE
  v
}

#' Euler axis-angle to Rotation matrix
#' 
#' @export
Euler2rotationMatrix <- function(v, proper=TRUE){
  n <- v[1:3]
  a <- v[4]
  co <- cos(a)
  si <- sin(a)
  if(proper){
    r1 <- c(co + n[1]^2*(1-co), n[1]*n[2]*(1-co)-n[3]*si, n[1]*n[3]*(1-co)+n[2]*si)
    r2 <- c(n[1]*n[2]*(1-co)+n[3]*si, co+n[2]^2*(1-co), n[2]*n[3]*(1-co)-n[1]*si)
    r3 <- c(n[1]*n[3]*(1-co)-n[2]*si, n[2]*n[3]*(1-co)+n[1]*si, co+n[3]^2*(1-co))
  } else{
    r1 <- c(co - n[1]^2*(1+co), -n[1]*n[2]*(1+co)-n[3]*si, -n[1]*n[3]*(1+co)+n[2]*si)
    r2 <- c(-n[1]*n[2]*(1+co)+n[3]*si, co-n[2]^2*(1+co), -n[2]*n[3]*(1+co)-n[1]*si)
    r3 <- c(-n[1]*n[3]*(1+co)-n[2]*si, -n[2]*n[3]*(1+co)+n[1]*si, co-n[3]^2*(1+co))  
  }
  
  R <- rbind(r1, r2, r3)
  unname(R)
}

#' Rotation matrix to Euler angles
#' 
#' @export
rotationMatrix2EulerAngles <- function(R){ 
  quaternion2EulerAngles(rotationMatrix2quaternion(R))
}

#' Euler angles to Rotation matrix
#' @export
EulerAngles2rotationMatrix <- function(angle){
  quaternion2rotationMatrix(EulerAngles2quaternion(angle))
}


#' perkele

# xyz<-u <- unit.vecs(3)
# e<-xyz2ll(u)
# f <- ll2xyz(e)
# all.equal(u,f)
# 

