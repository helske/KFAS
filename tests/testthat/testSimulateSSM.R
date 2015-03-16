test_that("FitSSM, transformSSM and simulateSSM works",{
  tol<-1e-4
  
  model<-SSModel(log(cbind(front,rear))~ -1 + log(PetrolPrice) + log(kms)
                 + SSMregression(~law,data=Seatbelts,index=1)
                 + SSMcustom(Z=diag(2),T=diag(2),R=matrix(1,2,1),
                             Q=matrix(1),P1=diag(10,2))
                 + SSMseasonal(period=12,sea.type='trigonometric'),
                 data=Seatbelts,H=matrix(NA,2,2)) 
  model$y[c(6,205,50,50+192)] <- NA
  updatefn <- function(pars, model,...){
    diag(model$H[,,1])<-exp(0.5*pars[1:2])
    model$H[1,2,1]<-model$H[2,1,1]<-tanh(pars[3])*prod(sqrt(exp(0.5*pars[1:2])))
    model$R[28:29]<-exp(pars[4:5])   
    model
  }
  expect_that(fit<-fitSSM(model=model,updatefn=updatefn,inits=c(-10,-10,0.7,-4,-4),method="BFGS"),not(gives_warning()))
  expect_equal(fit$optim.out$p,c(-10.431,-9.460,0.732,-4.182,-4.218),tolerance=1e-4,check.attributes=FALSE)
  expect_equal(fit$optim.out$value,-337.3011,tolerance=1e-6,check.attributes=FALSE)
  
  expect_that(mod <- transformSSM(fit$model, "augment"),not(gives_warning()))
  
  set.seed(123)  
  
  expect_that(out <- KFS(mod,smoothing="state"),not(gives_warning()))  
  expect_that(sim <- simulateSSM(mod,"states",nsim=25,antithetics=TRUE),not(gives_warning()))
  
  expect_equal(rowMeans(sim[192,,]), out$alpha[192,],tolerance=tol)
  expect_equal(rowMeans(sim[1,,]), out$alpha[1,],tolerance=tol)
  expect_equal(cov(t(sim[1,,])), out$V[,,1],tolerance=1e-2,check.attributes=FALSE)
  expect_equal(var(sim[1,29,]), 1.102994,tolerance=tol,check.attributes=FALSE)
  expect_equal(cov(t(sim[192,,])), out$V[,,192],tolerance=1e-2,check.attributes=FALSE)
  expect_equal(var(sim[192,29,]), 1.219971, tolerance=tol,check.attributes=FALSE)
  
  expect_that(out <- KFS(fit$model,smoothing=c("state","disturbance")),not(gives_warning()))
  expect_that(sim <- simulateSSM(fit$model,"disturbances",nsim=25,antithetics=T),not(gives_warning()))
  
  expect_equal(rowMeans(sim[192,,]), c(out$eps[192,],0), check.attributes=FALSE)
  expect_equal(rowMeans(sim[1,,]), c(out$eps[1,],out$eta[1,]), check.attributes=FALSE)
  expect_equal(apply(sim[1,1:2,],1,var), out$V_eps[,1],tolerance=0.1,check.attributes=FALSE)
  expect_equal(var(sim[1,3,]), 1.122758,tolerance=tol,check.attributes=FALSE)
  expect_equal(apply(sim[192,1:2,],1,var), out$V_eps[,192],tolerance=0.1,check.attributes=FALSE)
  expect_equal(var(sim[192,3,]), 0.7674087,tolerance=tol,check.attributes=FALSE)
})
