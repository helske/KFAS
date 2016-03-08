test_that("SSModel works properly",{
  tol<-1e-3
  set.seed(123)
  d<-data.frame(x=rnorm(100))
  t12<-ts(cbind(t1=rnorm(100)+d$x,t2=rnorm(100)))
  t12[sample(size=50,1:200)]<-NA
  expect_warning(
    model<-SSModel(t12~SSMcycle(period=10, type='common',Q=2)
                   +SSMcycle(period=10, type='distinct',P1=diag(c(1,1,2,2)),Q=diag(1:2))
                   +SSMtrend(2,type="common",Q=diag(c(1,0.5)))
                   +SSMtrend(2,type="distinct",Q=list(diag(0.1,2),diag(0.01,2)),
                             P1=diag(c(0.1,0.01,0.1,0.01)))
                   +SSMseasonal(period=4,type="common")
                   +SSMseasonal(period=4,type="distinct",Q=diag(c(2,3)),P1=diag(c(2,2,2,3,3,3)))
                   +SSMseasonal(period=5,type="common",sea.type="trig")
                   +SSMseasonal(period=5,type="distinct",sea.type="trig",Q=diag(c(0.1,0.2)),
                                P1=diag(rep(c(0.1,0.2),each=4)))
                   +SSMarima(ar=0.9,ma=0.2)+SSMregression(~-1+x,index=1,Q=1,data=d)
    ), NA)
  expect_warning(print(model), NA)
  expect_warning(logLik(model), NA)
  expect_equal(logLik(model),-442.006705531500,tolerance=tol,check.attributes=FALSE)
  expect_warning(out<-KFS(model,filtering=c("state","mean"),smoothing=c("state","mean","disturbance")), NA)
  expect_equal(out$d,11)
  expect_equal(out$j,1)
})
