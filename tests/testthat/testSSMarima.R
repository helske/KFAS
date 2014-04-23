test_that("arimaSSM works properly",{
s <- 12
phis <- 0.99
phi1 <- 0.0001
phi <- c(phi1,rep(0,s-2),phis,-phi1*phis)
theta <-0.7
out <- makeARIMA(phi,theta,NULL,SSinit="Ross")
min(eigen(out$Pn)$value)
out2<-SSMarima(phi,theta)
expect_equivalent(out$Pn,out2$P1)
expect_equivalent(out$T,out2$T)
expect_equivalent(out$Z,c(out2$Z))
})