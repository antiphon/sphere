#' test bingham's test
library(devtools)
load_all(".")

pv<-NULL
for(i in 1:5000){
x <- runifsphere(10)
pb <- bingham.test3d(t(x))
pu <- test.unifsphere(xyz2ai(x))
  pv <- rbind(pv, data.frame(pb=pb$p.value, 
                             plat=pu[1], plon=pu[2]))
}

par(mfrow=c(3,1))
res<-NULL
alpha <- 0.05
for(i in 1:3) {
  hist(pv[,i])
  nrej <- sum(pv[,i]<alpha)
  m <- length(pv[,i])
  res <- rbind(res, data.frame(m=m, rejected=nrej, rate=nrej/m, method=i))
}

print(res)
