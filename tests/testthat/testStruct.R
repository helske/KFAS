context("Structural time series tests")

tol<-1e-5

test_that("StructTS and KFS give equivalent results",{  
  # This only works for models with one state, as StructTS defines initial covariance matrix weirdly
  y <- Nile
  y[c(2,15,16,50,100)] <-NA
  fit <- StructTS(y, type = "level")
  # Construct model per StructTS:
  model<-SSModel(y ~ SSMtrend(1, fit$coef[1], P1 = fit$model0$P, a1 = fit$model0$a), H = fit$coef[2])
  out<-KFS(model,filtering="state",smoothing="none")
  expect_equal(c(fitted(fit)),coef(out,filtered=TRUE)[-1],tolerance=tol, check.attributes=FALSE)
  
  model<-SSModel(y~SSMtrend(degree=1,Q=NA),H=NA)
  fit2<-fitSSM(model,inits=c(7,10),method="BFGS")
  expect_equal(fit2$model["Q",1],fit$coef[1],tolerance=tol, check.attributes=FALSE)
  expect_equal(fit2$model["H",1],fit$coef[2],tolerance=tol, check.attributes=FALSE)

  pred <- predict(fit, n.ahead = 5)
  expect_equal(pred$pred, predict(fit2$model, n.ahead = 5),tolerance=tol, check.attributes=FALSE)  
})
