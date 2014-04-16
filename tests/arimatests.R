s <- 12
phis <- 0.99
phi1 <- 0.0001
phi <- c(phi1,rep(0,s-2),phis,-phi1*phis)
theta <- 0.7
out <- makeARIMA(phi,theta,NULL,SSinit="Ross")
min(eigen(out$Pn)$value)

set.seed(1)
 x <- arima.sim(100,model=list(ar=phi,ma=theta))
KalmanLike(x,mod=out)
out<-arima(x,c(13,0,1)
model<-SSModel(x~-1+SSMarima(ar=phi,ma=theta),H=0)
logLik(model)
arima(x,order=c(1,0,1),seasonal=list(period=12,order=c(1,0,0)),include.mean=FALSE,init=c(phi1,theta,phis),method='ML',
      SSinit="Ros")

Q0ter <- function(phi,theta){
  p <- length(phi)
  q <- length(theta)
  r <- max(p,q+1)
  T <- matrix(0,r,r)
  if (p) T[1:p,1] <- phi
  if (r) T[1:(r-1),2:r] <- diag(r-1)
  V <- matrix(0,r,r)
  ttheta <- c(1,theta)
  V[1:(q+1),1:(q+1)] <- ttheta%x%t(ttheta)
  V <- matrix(V,ncol=1)
  S <- diag(r*r)-T%x%T
  Q0 <- solve(S,V)
  Q0 <- matrix(Q0,ncol=r)
}
microbenchmark(makeARIMA(phi,theta,NULL,SSinit="Ross")$Pn,SSMarima(ar=phi,ma=theta)$Q,Q0ter(phi,theta))
