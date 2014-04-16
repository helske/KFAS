context("Structural time series tests")

test_that("SSModel works properly",{
  set.seed(123)
  t12<-ts(cbind(t1=rnorm(100),t2=rnorm(100)))
  t12[sample(size=50,1:200)]<-NA
  expect_that( 
    model<-SSModel(t12~SSMcycle(period=10, type='common',Q=2)
                   +SSMcycle(period=10, type='distinct',P1=diag(c(1,1,2,2)),Q=diag(1:2))
                   +SSMtrend(2,type="common",Q=diag(c(1,0.5)))
                   +SSMtrend(2,type="distinct",Q=list(diag(0.1,2),diag(0.01,2)),P1=diag(c(0.1,0.01,0.1,0.01)))
                   +SSMseasonal(period=4,type="common")
                   +SSMseasonal(period=4,type="distinct",Q=diag(c(2,3)),P1=diag(c(2,2,2,3,3,3)))
                   +SSMseasonal(period=5,type="common",sea.type="trig")
                   +SSMseasonal(period=5,type="distinct",sea.type="trig",Q=diag(c(0.1,0.2)),
                                P1=diag(rep(c(0.1,0.2),each=4)))
                   +SSMarima(ar=0.9,ma=0.2)+SSMregression(~-1+t2,index=1,Q=1)
    ),not(gives_warning()))
  expect_that(print(model),not(gives_warning()))
  expect_that(logLik(model),not(gives_warning()))
  expect_equivalent(logLik(model),-442.270604625564)
  expect_that(out<-KFS(model,filtering=c("state","mean"),smoothing=c("state","mean","disturbance")),
              not(gives_warning()))
  expect_equal(out$d,11)
  expect_equal(out$j,1)
})

test_that("StructTS and KFS give equivalent results",{
require(graphics)
trees <- window(treering, start = 0)
fit <- StructTS(trees, type = "level")
# Construct model per StructTS:
model<-SSModel(trees~SSMtrend(degree=1,Q=fit$coef[1],P1=fit$model0$P,a1=fit$model0$a),H=fit$coef[2])
out<-KFS(model)
expect_equivalent(c(fitted(fit)),coef(out,filtered=TRUE)[-1])
expect_equivalent(fit$loglik,logLik(model))

fit2 <- StructTS(trees, type = "level",fixed=c(NA,1e-15))
model<-SSModel(trees~SSMtrend(degree=1,Q=fit2$coef[1],P1=fit2$model0$P,a1=fit2$model0$a),H=fit2$coef[2])
expect_less_than(abs(fit2$loglik-logLik(model)),1e-4)
fit2 <- StructTS(trees, type = "level",fixed=c(1e-15,NA))
model<-SSModel(trees~SSMtrend(degree=1,Q=fit2$coef[1],P1=fit2$model0$P,a1=fit2$model0$a),H=fit2$coef[2])
expect_less_than(abs(fit2$loglik-logLik(model)),1e-11)

model<-SSModel(trees~SSMtrend(degree=1,Q=NA),H=NA)
fit2<-fitSSM(model,inits=c(0,0),method="BFGS")
expect_equivalent(c(fit2$model$H),fit$coef[2])
expect_less_than(abs(c(fit2$model$Q)-fit$coef[1]),1e-9)
