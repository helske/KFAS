context("KFAS and glm comparison")

# Test for Gaussian
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
glm.gaussian <- glm(weight ~ group)
model.gaussian <- SSModel(weight~group)
tmp<-KFS(model.gaussian,filtering="state",smoothing="none")
model.gaussian <- SSModel(weight~group,H=mean(c(tmp$v[1:11][tmp$Finf==0]^2/tmp$F[1:11][tmp$Finf==0],tmp$v[12:20]^2/tmp$F[12:20])))
kfas.gaussian <-KFS(model.gaussian,smoothing=c('state','signal','mean'))

# Test for Poisson GLM
# Same example as in ?glm
d <- data.frame(treatment= gl(3,3), outcome = gl(3,1,9), counts = c(18,17,15,20,10,20,25,13,12))
glm.poisson <- glm(counts ~ outcome + treatment, data=d, family = poisson(),control=list(epsilon=1e-15))
model.poisson<-SSModel(counts ~ outcome + treatment, data=d, distribution = 'poisson')
kfas.poisson<-KFS(model.poisson,smoothing=c('state','signal','mean'))

## Test for Binomial GLM
## example from Venables and Ripley (2002, pp. 190-2.)
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive = 20-numdead)
glm.binomial <- glm(SF ~ sex*ldose, family = binomial,control=list(epsilon=1e-15))
model.binomial<-SSModel(numdead ~ sex*ldose, u=20,distribution = 'binomial')
kfas.binomial<-KFS(model.binomial,smoothing=c('state','signal','mean'))

## Test for Gamma GLM
# A Gamma example from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
  u = c(5,10,15,20,30,40,60,80,100),
  lot1 = c(118,58,42,35,27,25,21,19,18),
  lot2 = c(69,35,26,21,18,16,13,12,12))
glm.gamma<-glm(lot1 ~ log(u), data = clotting, family = Gamma("log"),control=list(epsilon=1e-15))
#dispersion=1
model.gamma1<-SSModel(lot1 ~ log(u), data = clotting, distribution = 'gamma')
kfas.gamma1<-KFS(model.gamma1,smoothing=c('state','signal','mean'))
#dispersion from gamma.glm
model.gamma2<-SSModel(lot1 ~ log(u), u = 1/summary(glm.gamma)$dispersion, data = clotting, distribution = 'gamma')
kfas.gamma2<-KFS(model.gamma2,smoothing=c('state','signal','mean'))  

## Test for NB GLM
## From MASS library, ?glm.nb
glm.NB <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine,control=glm.control(epsilon=1e-12))
#u from theta
model.NB<-SSModel(Days ~ Sex/(Age + Eth*Lrn), u=glm.NB$theta, data = quine, distribution = 'negative binomial')
kfas.NB<-KFS(model.NB,smoothing=c('state','signal','mean'))

test_that("Gaussian GLM fitting works properly",{
  for(i in 1:length(model.gaussian$y))
    expect_equal(coef(kfas.gaussian,start=i,end=i),coef(glm.gaussian),info=paste("Error in time step",i))
  for(i in 1:length(model.gaussian$y))
    expect_equivalent(kfas.gaussian$V[,,i],summary(glm.gaussian)$cov.s,info=paste("Error in time step",i))
  expect_equivalent(c(fitted(kfas.gaussian)),fitted(glm.gaussian))
  expect_equivalent(c(kfas.gaussian$V_mu),predict(glm.gaussian,type='response',se.fit=TRUE)$se.fit^2)
  expect_equal(deviance(kfas.gaussian),deviance(glm.gaussian))
})

test_that("Poisson GLM fitting works properly",{
  for(i in 1:length(model.poisson$y))
    expect_equivalent(coef(kfas.poisson,start=i,end=i),coef(glm.poisson),info=paste("Error in time step",i))
  for(i in 1:length(model.poisson$y))
    expect_equivalent(kfas.poisson$V[,,i],summary(glm.poisson)$cov.s,info=paste("Error in time step",i))
  ## linear
  expect_equivalent(c(kfas.poisson$thetahat),glm.poisson$linear.predictor)
  ## predictions on response scale
  expect_equivalent(c(fitted(kfas.poisson)),fitted(glm.poisson))
  # prediction variances on link scale
  expect_equivalent(c(kfas.poisson$V_theta),predict(glm.poisson,type='link',se.fit=TRUE)$se.fit^2)
  # prediction variances on response scale
  expect_equivalent(c(kfas.poisson$V_mu),predict(glm.poisson,type='response',se.fit=TRUE)$se.fit^2)
  expect_equal(deviance(kfas.poisson),deviance(glm.poisson))
})


test_that("binomial GLM fitting works properly",{
  for(i in 1:length(model.binomial$y))
    expect_equal(coef(kfas.binomial,start=i,end=i),coef(glm.binomial),info=paste("Error in time step",i))
  for(i in 1:length(model.binomial$y))
    expect_equivalent(kfas.binomial$V[,,i],summary(glm.binomial)$cov.s,info=paste("Error in time step",i))
  ## linear
  expect_equivalent(c(kfas.binomial$thetahat),glm.binomial$linear.predictor)
  ## predictions on response scale
  expect_equivalent(c(fitted(kfas.binomial)),fitted(glm.binomial))
  # prediction variances on link scale
  expect_equivalent(c(kfas.binomial$V_theta),predict(glm.binomial,type='link',se.fit=TRUE)$se.fit^2)
  # prediction variances on response scale
  expect_equivalent(c(kfas.binomial$V_mu),predict(glm.binomial,type='response',se.fit=TRUE)$se.fit^2)
  expect_equal(deviance(kfas.binomial),deviance(glm.binomial))
})


test_that("gamma GLM fitting works properly",{
  for(i in 1:length(model.gamma1$y))
    expect_equal(coef(kfas.gamma1,start=i,end=i),coef(glm.gamma),info=paste("Error in time step",i, ", dispersion=1"))
  for(i in 1:length(model.gamma1$y))
    expect_equivalent(kfas.gamma1$V[,,i],summary(glm.gamma,dispersion=1)$cov.s,
                      info=paste("Error in time step",i, ", dispersion=1"))  
  expect_equivalent(c(kfas.gamma1$thetahat),glm.gamma$linear.predictor,
                    info=paste("dispersion=1"))
  ## predictions on response scale
  expect_equivalent(c(fitted(kfas.gamma1)),fitted(glm.gamma),info=paste("dispersion=1"))
  # prediction variances on link scale
  expect_equivalent(c(kfas.gamma1$V_theta),
                    predict(glm.gamma,type='link',se.fit=TRUE,dispersion=1)$se.fit^2,
                    info=paste("dispersion=1"))
  # prediction variances on response scale
  expect_equivalent(c(kfas.gamma1$V_mu),predict(glm.gamma,type='response',se.fit=TRUE,dispersion=1)$se.fit^2,
                    info=paste("dispersion=1"))
  for(i in 1:length(model.gamma2$y))
    expect_equal(coef(kfas.gamma2,start=i,end=i),coef(glm.gamma),info=paste("Error in time step",i))
  for(i in 1:length(model.gamma2$y))
    expect_equivalent(kfas.gamma2$V[,,i],summary(glm.gamma)$cov.s,info=paste("Error in time step",i))
  ## linear
  expect_equivalent(c(kfas.gamma2$theta),glm.gamma$linear.predictor)
  ## predictions on response scale
  expect_equivalent(c(fitted(kfas.gamma2)),fitted(glm.gamma))
  # prediction variances on link scale
  expect_equivalent(c(kfas.gamma2$V_theta),predict(glm.gamma,type='link',se.fit=TRUE)$se.fit^2)
  # prediction variances on response scale
  expect_equivalent(c(kfas.gamma2$V_mu),predict(glm.gamma,type='response',se.fit=TRUE)$se.fit^2)
  expect_equal(deviance(kfas.gamma2),deviance(glm.gamma))
})

test_that("negative binomial GLM fitting works properly",{
  
  for(i in 1:length(model.NB$y))
    expect_equal(coef(kfas.NB,start=i,end=i),coef(glm.NB),info=paste("Error in time step",i))
  for(i in 1:length(model.NB$y))
    expect_equivalent(kfas.NB$V[,,i],summary(glm.NB)$cov.s,
                      info=paste("Error in time step",i))
  
  expect_equivalent(c(kfas.NB$thetahat),glm.NB$linear.predictor)
  ## predictions on response scale
  expect_equivalent(c(fitted(kfas.NB)),fitted(glm.NB))
  # prediction variances on link scale
  expect_equivalent(c(kfas.NB$V_theta),
                    predict(glm.NB,type='link',se.fit=TRUE,dispersion=1)$se.fit^2)
  # prediction variances on response scale
  expect_equivalent(c(kfas.NB$V_mu),predict(glm.NB,type='response',se.fit=TRUE,)$se.fit^2)
  
  likfn<-function(pars,model,estimate=TRUE){ 
    model$u[]<-exp(pars[1])
    model$P1inf[]<-0
    model$a1[]<-pars[2:15]
    if(estimate)
      return(-logLik(model,check=TRUE,nsim=0))
    model
  }
  
  fit<-optim(f=likfn,p=c(log(glm.NB$theta),glm.NB$coef),model=model.NB,method="BFGS")
  expect_equal(c(exp(fit$p[1]),fit$p[-1]),c(glm.NB$theta,glm.NB$coef))
  expect_equal(deviance(kfas.NB),deviance(glm.NB))
})


test_that("Residuals for Gaussian GLM works properly",{
  expect_equivalent(as.numeric(residuals(kfas.gaussian,type="deviance")),
                    residuals(glm.gaussian,type="deviance"))
  expect_equivalent(as.numeric(residuals(kfas.gaussian,type="pearson")),
                    residuals(glm.gaussian,type="pearson"))
  expect_equivalent(as.numeric(residuals(kfas.gaussian,type="response")),
                    residuals(glm.gaussian,type="response"))
  
  expect_equivalent(as.numeric(rstandard(kfas.gaussian,type="pearson")),
                    rstandard(glm.gaussian,type="pearson"))
  expect_equivalent(as.numeric(rstandard(kfas.gaussian,type="deviance")),
                    rstandard(glm.gaussian,type="deviance"))
})


test_that("Residuals for Poisson GLM works properly",{
  expect_equivalent(as.numeric(residuals(kfas.poisson,type="deviance")),
                    residuals(glm.poisson,type="deviance"))
  expect_equivalent(as.numeric(residuals(kfas.poisson,type="pearson")),
                    residuals(glm.poisson,type="pearson"))
  expect_equivalent(as.numeric(residuals(kfas.poisson,type="response")),
                    residuals(glm.poisson,type="response"))
  
  expect_equivalent(as.numeric(rstandard(kfas.poisson,type="pearson")),
                    rstandard(glm.poisson,type="pearson"))
  expect_equivalent(as.numeric(rstandard(kfas.poisson,type="deviance")),
                    rstandard(glm.poisson,type="deviance"))
})


test_that("Residuals for Binomial GLM works properly",{
  expect_equivalent(as.numeric(residuals(kfas.binomial,type="deviance")),
                    residuals(glm.binomial,type="deviance"))
  expect_equivalent(as.numeric(residuals(kfas.binomial,type="pearson")),
                    residuals(glm.binomial,type="pearson"))
  expect_equivalent(as.numeric(residuals(kfas.binomial,type="response")),
                    residuals(glm.binomial,type="response"))
  
  expect_equivalent(as.numeric(rstandard(kfas.binomial,type="pearson")),
                    rstandard(glm.binomial,type="pearson"))
  expect_equivalent(as.numeric(rstandard(kfas.binomial,type="deviance")),
                    rstandard(glm.binomial,type="deviance"))
})


test_that("Residuals for Gamma GLM works properly",{
  expect_equivalent(as.numeric(residuals(kfas.gamma2,type="deviance")),
                    residuals(glm.gamma,type="deviance"))
  expect_equivalent(as.numeric(residuals(kfas.gamma2,type="pearson")),
                    residuals(glm.gamma,type="pearson"))
  expect_equivalent(as.numeric(residuals(kfas.gamma2,type="response")),
                    residuals(glm.gamma,type="response"))
  
  expect_equivalent(as.numeric(rstandard(kfas.gamma2,type="pearson")),
                    rstandard(glm.gamma,type="pearson"))
  expect_equivalent(as.numeric(rstandard(kfas.gamma2,type="deviance")),
                    rstandard(glm.gamma,type="deviance"))
})


test_that("Residuals for negative binomial GLM works properly",{
  expect_equivalent(as.numeric(residuals(kfas.NB,type="deviance")),
                    residuals(glm.NB,type="deviance"))
  expect_equivalent(as.numeric(residuals(kfas.NB,type="pearson")),
                    residuals(glm.NB,type="pearson"))
  expect_equivalent(as.numeric(residuals(kfas.NB,type="response")),
                    residuals(glm.NB,type="response"))
  
  expect_equivalent(as.numeric(rstandard(kfas.NB,type="pearson")),
                    rstandard(glm.NB,type="pearson"))
  expect_equivalent(as.numeric(rstandard(kfas.NB,type="deviance")),
                    rstandard(glm.NB,type="deviance"))
})


test_that("Predictions for GLM works properly",{
  pred.glm.gaussian.link<-predict(glm.gaussian,type="link",se.fit=TRUE)
  pred.kfas.gaussian.link<-predict(model.gaussian,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.gaussian.response<-predict(glm.gaussian,type="response",se.fit=TRUE)
  pred.kfas.gaussian.response<-predict(model.gaussian,type="response",se.fit=TRUE,interval="confidence")  
  expect_equivalent(c(pred.kfas.gaussian.link[,"fit"]),pred.glm.gaussian.link$fit)
  expect_equivalent(c(pred.kfas.gaussian.link[,"se.fit"]),pred.glm.gaussian.link$se.fit)
  expect_equivalent(c(pred.kfas.gaussian.response[,"fit"]),pred.glm.gaussian.response$fit)
  expect_equivalent(c(pred.kfas.gaussian.response[,"se.fit"]),pred.glm.gaussian.response$se.fit)
  
  
  pred.glm.poisson.link<-predict(glm.poisson,type="link",se.fit=TRUE)
  pred.kfas.poisson.link<-predict(model.poisson,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.poisson.response<-predict(glm.poisson,type="response",se.fit=TRUE)
  pred.kfas.poisson.response<-predict(model.poisson,type="response",se.fit=TRUE,interval="confidence")  
  expect_equivalent(pred.glm.poisson.link$fit,obj=c(pred.kfas.poisson.link[,"fit"]))
  expect_equivalent(pred.glm.poisson.link$se.fit,obj=c(pred.kfas.poisson.link[,"se.fit"]))
  expect_equivalent(pred.glm.poisson.response$fit,obj=c(pred.kfas.poisson.response[,"fit"]))
  expect_equivalent(pred.glm.poisson.response$se.fit,obj=c(pred.kfas.poisson.response[,"se.fit"]))
  
  
  pred.glm.binomial.link<-predict(glm.binomial,type="link",se.fit=TRUE)
  pred.kfas.binomial.link<-predict(model.binomial,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.binomial.response<-predict(glm.binomial,type="response",se.fit=TRUE)
  pred.kfas.binomial.response<-predict(model.binomial,type="response",se.fit=TRUE,interval="confidence")  
  expect_equivalent(pred.glm.binomial.link$fit,obj=c(pred.kfas.binomial.link[,"fit"]))
  expect_equivalent(pred.glm.binomial.link$se.fit,obj=c(pred.kfas.binomial.link[,"se.fit"]))
  expect_equivalent(pred.glm.binomial.response$fit,obj=c(pred.kfas.binomial.response[,"fit"]))
  expect_equivalent(pred.glm.binomial.response$se.fit,obj=c(pred.kfas.binomial.response[,"se.fit"]))
  
  
  pred.glm.gamma.link<-predict(glm.gamma,type="link",se.fit=TRUE)
  pred.kfas.gamma.link<-predict(model.gamma2,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.gamma.response<-predict(glm.gamma,type="response",se.fit=TRUE)
  pred.kfas.gamma.response<-predict(model.gamma2,type="response",se.fit=TRUE,interval="confidence")  
  expect_equivalent(pred.glm.gamma.link$fit,obj=c(pred.kfas.gamma.link[,"fit"]))
  expect_equivalent(pred.glm.gamma.link$se.fit,obj=c(pred.kfas.gamma.link[,"se.fit"]))
  expect_equivalent(pred.glm.gamma.response$fit,obj=c(pred.kfas.gamma.response[,"fit"]))
  expect_equivalent(pred.glm.gamma.response$se.fit,obj=c(pred.kfas.gamma.response[,"se.fit"]))
  
  
  pred.glm.NB.link<-predict(glm.NB,type="link",se.fit=TRUE)
  pred.kfas.NB.link<-predict(model.NB,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.NB.response<-predict(glm.NB,type="response",se.fit=TRUE)
  pred.kfas.NB.response<-predict(model.NB,type="response",se.fit=TRUE,interval="confidence")  
  expect_equivalent(pred.glm.NB.link$fit,obj=c(pred.kfas.NB.link[,"fit"]))
  expect_equivalent(pred.glm.NB.link$se.fit,obj=c(pred.kfas.NB.link[,"se.fit"]))
  expect_equivalent(pred.glm.NB.response$fit,obj=c(pred.kfas.NB.response[,"fit"]))
  expect_equivalent(pred.glm.NB.response$se.fit,obj=c(pred.kfas.NB.response[,"se.fit"]))
  
  #confidence intervals
  glm.gaussian.CI<-pred.glm.gaussian.link$fit +qnorm(0.025)*pred.glm.gaussian.link$se.fit%o%c(1,-1)
  expect_equivalent(unclass(pred.kfas.gaussian.link[,c("lwr","upr")]),glm.gaussian.CI)
  expect_equivalent(unclass(pred.kfas.gaussian.response[,c("lwr","upr")]),gaussian()$linkinv(glm.gaussian.CI))
  
  glm.poisson.CI<-pred.glm.poisson.link$fit +qnorm(0.025)*pred.glm.poisson.link$se.fit%o%c(1,-1)
  expect_equivalent(unclass(pred.kfas.poisson.link[,c("lwr","upr")]),glm.poisson.CI)
  expect_equivalent(unclass(pred.kfas.poisson.response[,c("lwr","upr")]),poisson()$linkinv(glm.poisson.CI))
  
  glm.binomial.CI<-pred.glm.binomial.link$fit +qnorm(0.025)*pred.glm.binomial.link$se.fit%o%c(1,-1)
  expect_equivalent(unclass(pred.kfas.binomial.link[,c("lwr","upr")]),glm.binomial.CI)
  expect_equivalent(unclass(pred.kfas.binomial.response[,c("lwr","upr")]),binomial()$linkinv(glm.binomial.CI))
  
  glm.gamma.CI<-pred.glm.gamma.link$fit +qnorm(0.025)*pred.glm.gamma.link$se.fit%o%c(1,-1)
  expect_equivalent(unclass(pred.kfas.gamma.link[,c("lwr","upr")]),glm.gamma.CI)
  expect_equivalent(unclass(pred.kfas.gamma.response[,c("lwr","upr")]),Gamma(link="log")$linkinv(glm.gamma.CI))
  
  glm.NB.CI<-pred.glm.NB.link$fit +qnorm(0.025)*pred.glm.NB.link$se.fit%o%c(1,-1)
  expect_equivalent(unclass(pred.kfas.NB.link[,c("lwr","upr")]),glm.NB.CI)
  expect_equivalent(unclass(pred.kfas.NB.response[,c("lwr","upr")]),glm.NB$family$linkinv(glm.NB.CI))
})
##

# #Estimate dispersion:
# 
# x<-rnorm(10000)
# beta<-2
# dispersion<-1
# y<-rgamma(n=length(x),shape=1/dispersion,scale=exp(beta*x)*dispersion)
# 
# glmfit<-glm(y~x,family=Gamma("log"),control=list(epsilon=1e-15,maxit=10000))
# 
# model<-SSModel(y~x,distribution="gamma")
# likfn_u<-function(pars,model,estimate=TRUE){ 
#   model$u[]<-exp(pars)
#   if(estimate)
#     return(-logLik(model,check=TRUE))
#   model
# }
# fit_u<-nlm(f=likfn_u,p=log(dispersion),model=model)
# sort(abs(c(glm=summary(glmfit)$disp,MASS=1/gamma.shape(glmfit,it.lim=1000,eps=1e-15)$alpha,
#            KFAS=exp(fit_u$e))-dispersion))

#imp<-importanceSSM(model.gamma2,nsim=10000,antit=T,type="signal")
##imp$s<-exp(imp$s)
#weighted.hist(imp$s[1,1,],imp$w)

# imp<-importanceSSM(model.binomial,nsim=10000,antit=TRUE,type="signal")
# imp$s<-exp(imp$s)/(1 + exp(imp$s))
# weighted.hist(imp$s[10,1,],imp$w)
