test_that("arimaSSM works properly",{
  tol<-1e-3
  s <- 12
  phis <- 0.99
  phi1 <- 0.0001
  phi <- c(phi1,rep(0,s-2),phis,-phi1*phis)
  theta <-0.7
  out <- makeARIMA(phi,theta,NULL,SSinit="Ross")
  min(eigen(out$Pn)$value)
  out2<-SSMarima(phi,theta)
  expect_equal(out$Pn,out2$P1,tolerance=tol,check.attributes=FALSE)
  expect_equal(out$T,out2$T,tolerance=tol,check.attributes=FALSE)
  expect_equal(out$Z,c(out2$Z),tolerance=tol,check.attributes=FALSE)
})
