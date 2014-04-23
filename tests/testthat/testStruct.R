context("Structural time series tests")



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
})