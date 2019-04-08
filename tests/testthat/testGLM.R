context("KFAS and glm comparison")

# Poor starting values
y <- rep(0:1,c(15,10))
model<-SSModel(y~1,dist="binomial")
expect_equivalent(coef(KFS(model,theta=7.1)),coef(KFS(model)))
expect_equivalent(coef(KFS(model,theta=-4.6)),coef(KFS(model)))

tol<-1e-4
require(MASS)
# Test for Gaussian
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
glm.gaussian <- glm(weight ~ group)
model.gaussian <- SSModel(weight~group)
tmp<-KFS(model.gaussian,filtering="state",smoothing="none")
model.gaussian <- SSModel(weight~group,H=mean(c(tmp$v[1:tmp$d][tmp$Finf==0]^2/tmp$F[1:tmp$d][tmp$Finf==0],
                                                tmp$v[-(1:tmp$d)]^2/tmp$F[-(1:tmp$d)])))
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
kfas.binomial<-KFS(model.binomial,smoothing=c('state','signal','mean'),maxiter=1000,convtol=1e-15)
kfas.binomial2<-KFS(model.binomial,smoothing=c('state','signal','mean'),maxiter=1000)

## Test for Gamma GLM
# A Gamma example from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
  u = c(5,10,15,20,30,40,60,80,100),
  lot1 = c(118,58,42,35,27,25,21,19,18),
  lot2 = c(69,35,26,21,18,16,13,12,12))
glm.gamma<-glm(lot1 ~ log(u), data = clotting, family = Gamma("log"),control=list(epsilon=1e-8))
#dispersion=1
model.gamma1<-SSModel(lot1 ~ log(u), data = clotting, distribution = 'gamma')
kfas.gamma1<-KFS(model.gamma1,smoothing=c('state','signal','mean'))
#dispersion from gamma.glm
model.gamma2<-SSModel(lot1 ~ log(u), u = 1/summary(glm.gamma)$dispersion, data = clotting, distribution = 'gamma')
kfas.gamma2<-KFS(model.gamma2,smoothing=c('state','signal','mean'),maxiter=1000,convtol=1e-8)  

## Test for NB GLM
## From MASS library, ?glm.nb
glm.NB <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine,control=glm.control(epsilon=1e-12))

# estimate theta
theta0<-1
model.NB<-SSModel(Days ~ Sex/(Age + Eth*Lrn), u=theta0, data = quine, distribution = 'negative binomial')
kfas.NB<-KFS(model.NB,smoothing='mean')
theta<-theta.ml(y=quine$Days,mu=c(fitted(kfas.NB)))
maxit<-10
i<-0
while(abs(theta-theta0)/theta0>1e-7 && i<maxit){
  model.NB$u[]<-theta0<-theta
  kfas.NB<-KFS(model.NB,smoothing='mean')
  theta<-theta.ml(y=quine$Days,mu=c(fitted(kfas.NB)))
  i<-i+1
}
model.NB$u[]<-theta
kfas.NB<-KFS(model.NB,smoothing=c('state','signal','mean'),convtol=1e-15)

test_that("Gaussian GLM fitting works properly",{
  for(i in 1:length(model.gaussian$y))
    expect_equal(coef(kfas.gaussian,start=i,end=i),coef(glm.gaussian),
                 tolerance=tol, check.attributes=FALSE, info=paste("Error in time step",i))
  for(i in 1:length(model.gaussian$y))
    expect_equal(kfas.gaussian$V[,,i],summary(glm.gaussian)$cov.s,tolerance=tol, check.attributes=FALSE,
                 info=paste("Error in time step",i))
  expect_equal(c(fitted(kfas.gaussian)),fitted(glm.gaussian),tolerance=tol, check.attributes=FALSE)
  expect_equal(c(kfas.gaussian$V_mu),predict(glm.gaussian,type='response',se.fit=TRUE)$se.fit^2
               ,tolerance=tol, check.attributes=FALSE)
})

test_that("Poisson GLM fitting works properly",{
  for(i in 1:length(model.poisson$y))
    expect_equal(coef(kfas.poisson,start=i,end=i),coef(glm.poisson),info=paste("Error in time step",i)
                 ,tolerance=tol, check.attributes=FALSE)
  for(i in 1:length(model.poisson$y))
    expect_equal(kfas.poisson$V[,,i],summary(glm.poisson)$cov.s,info=paste("Error in time step",i)
                 ,tolerance=tol, check.attributes=FALSE)
  ## linear
  expect_equal(c(kfas.poisson$thetahat),glm.poisson$linear.predictor,tolerance=tol, check.attributes=FALSE)
  ## predictions on response scale
  expect_equal(c(fitted(kfas.poisson)),fitted(glm.poisson),tolerance=tol, check.attributes=FALSE)
  # prediction variances on link scale
  expect_equal(c(kfas.poisson$V_theta),predict(glm.poisson,type='link',se.fit=TRUE)$se.fit^2
               ,tolerance=tol, check.attributes=FALSE)
  # prediction variances on response scale
  expect_equal(c(kfas.poisson$V_mu),predict(glm.poisson,type='response',se.fit=TRUE)$se.fit^2
               ,tolerance=tol, check.attributes=FALSE)
})


test_that("binomial GLM fitting works properly",{
  for(i in 1:length(model.binomial$y))
    expect_equal(coef(kfas.binomial,start=i,end=i),coef(glm.binomial),
                 info=paste("Error in time step",i),tolerance=tol, check.attributes=FALSE)
  for(i in 1:length(model.binomial$y))
    expect_equal(kfas.binomial$V[,,i],summary(glm.binomial)$cov.s,
                 info=paste("Error in time step",i),tolerance=tol, check.attributes=FALSE)
  ## linear
  expect_equal(c(kfas.binomial$thetahat),glm.binomial$linear.predictor,tolerance=tol, check.attributes=FALSE)
  ## predictions on response scale
  expect_equal(c(fitted(kfas.binomial)),fitted(glm.binomial),tolerance=tol, check.attributes=FALSE)
  # prediction variances on link scale
  expect_equal(c(kfas.binomial$V_theta),
               predict(glm.binomial,type='link',se.fit=TRUE)$se.fit^2,tolerance=tol, check.attributes=FALSE)
  # prediction variances on response scale
  expect_equal(c(kfas.binomial$V_mu),
               predict(glm.binomial,type='response',se.fit=TRUE)$se.fit^2,tolerance=tol, check.attributes=FALSE)
})


test_that("gamma GLM fitting works properly",{
  for(i in 1:length(model.gamma1$y))
    expect_equal(coef(kfas.gamma1,start=i,end=i),coef(glm.gamma),
                 info=paste("Error in time step",i, ", dispersion=1"),tolerance=tol, check.attributes=FALSE)
  for(i in 1:length(model.gamma1$y))
    expect_equal(kfas.gamma1$V[,,i],summary(glm.gamma,dispersion=1)$cov.s,
                 info=paste("Error in time step",i, ", dispersion=1"),tolerance=tol, check.attributes=FALSE)  
  expect_equal(c(kfas.gamma1$thetahat),glm.gamma$linear.predictor,
               info=paste("dispersion=1"),tolerance=tol, check.attributes=FALSE)
  ## predictions on response scale
  expect_equal(c(fitted(kfas.gamma1)),fitted(glm.gamma),info=paste("dispersion=1"),tolerance=tol, check.attributes=FALSE)
  # prediction variances on link scale
  expect_equal(c(kfas.gamma1$V_theta),
               predict(glm.gamma,type='link',se.fit=TRUE,dispersion=1)$se.fit^2,
               info=paste("dispersion=1"),tolerance=tol, check.attributes=FALSE)
  # prediction variances on response scale
  expect_equal(c(kfas.gamma1$V_mu),predict(glm.gamma,type='response',se.fit=TRUE,dispersion=1)$se.fit^2,
               info=paste("dispersion=1"),tolerance=tol, check.attributes=FALSE)
  for(i in 1:length(model.gamma2$y))
    expect_equal(coef(kfas.gamma2,start=i,end=i),coef(glm.gamma),
                 info=paste("Error in time step",i),tolerance=tol, check.attributes=FALSE)
  for(i in 1:length(model.gamma2$y))
    expect_equal(kfas.gamma2$V[,,i],summary(glm.gamma)$cov.s,
                 info=paste("Error in time step",i),tolerance=tol, check.attributes=FALSE)
  ## linear
  expect_equal(c(kfas.gamma2$theta),glm.gamma$linear.predictor,tolerance=tol, check.attributes=FALSE)
  ## predictions on response scale
  expect_equal(c(fitted(kfas.gamma2)),fitted(glm.gamma),tolerance=tol, check.attributes=FALSE)
  # prediction variances on link scale
  expect_equal(c(kfas.gamma2$V_theta),
               predict(glm.gamma,type='link',se.fit=TRUE)$se.fit^2,tolerance=tol, check.attributes=FALSE)
  # prediction variances on response scale
  expect_equal(c(kfas.gamma2$V_mu),
               predict(glm.gamma,type='response',se.fit=TRUE)$se.fit^2,tolerance=tol, check.attributes=FALSE)
})

test_that("negative binomial GLM fitting works properly",{
  
  for(i in 1:length(model.NB$y))
    expect_equal(coef(kfas.NB,start=i,end=i),coef(glm.NB),info=paste("Error in time step",i),tolerance=tol, check.attributes=FALSE)
  for(i in 1:length(model.NB$y))
    expect_equal(kfas.NB$V[,,i],summary(glm.NB)$cov.s,
                 info=paste("Error in time step",i),tolerance=tol, check.attributes=FALSE)
  
  expect_equal(c(kfas.NB$thetahat),glm.NB$linear.predictor,tolerance=tol, check.attributes=FALSE)
  ## predictions on response scale
  expect_equal(c(fitted(kfas.NB)),fitted(glm.NB),tolerance=tol, check.attributes=FALSE)
  # prediction variances on link scale
  expect_equal(c(kfas.NB$V_theta),
               predict(glm.NB,type='link',se.fit=TRUE,dispersion=1)$se.fit^2,tolerance=tol, check.attributes=FALSE)
  # prediction variances on response scale
  expect_equal(c(kfas.NB$V_mu),predict(glm.NB,type='response',se.fit=TRUE,)$se.fit^2,tolerance=tol, check.attributes=FALSE)
  
  expect_equal(model.NB$u[1],glm.NB$theta,tolerance=tol, check.attributes=FALSE)
  
  likfn<-function(pars,model,estimate=TRUE){ 
    model$u[]<-exp(pars[1])
    model$P1inf[]<-0
    model$a1[]<-pars[2:15]
    if(estimate)
      return(-logLik(model,check=TRUE,nsim=0))
    model
  }
  
  fit<-optim(f=likfn,p=c(log(glm.NB$theta),glm.NB$coef),model=model.NB)
  expect_equal(c(exp(fit$p[1]),fit$p[-1]),c(glm.NB$theta,glm.NB$coef),tolerance=tol, check.attributes=FALSE)
})


test_that("Residuals for Gaussian GLM works properly",{
  expect_equal(as.numeric(residuals(kfas.gaussian,type="pearson")),
               residuals(glm.gaussian,type="pearson"),tolerance=tol, check.attributes=FALSE)
  expect_equal(as.numeric(residuals(kfas.gaussian,type="response")),
               residuals(glm.gaussian,type="response"),tolerance=tol, check.attributes=FALSE)
  
  expect_equal(as.numeric(rstandard(kfas.gaussian,type="pearson")),
               rstandard(glm.gaussian,type="pearson"),tolerance=tol, check.attributes=FALSE)
})


test_that("Residuals for Poisson GLM works properly",{
  expect_equal(as.numeric(residuals(kfas.poisson,type="pearson")),
               residuals(glm.poisson,type="pearson"),tolerance=tol, check.attributes=FALSE)
  expect_equal(as.numeric(residuals(kfas.poisson,type="response")),
               residuals(glm.poisson,type="response"),tolerance=tol, check.attributes=FALSE)
  
  expect_equal(as.numeric(rstandard(kfas.poisson,type="pearson")),
               rstandard(glm.poisson,type="pearson"),tolerance=tol, check.attributes=FALSE)
})


test_that("Residuals for Binomial GLM works properly",{
  expect_equal(as.numeric(residuals(kfas.binomial,type="pearson")),
               residuals(glm.binomial,type="pearson"),tolerance=tol, check.attributes=FALSE)
  expect_equal(as.numeric(residuals(kfas.binomial,type="response")),
               residuals(glm.binomial,type="response"),tolerance=tol, check.attributes=FALSE)
  
  expect_equal(as.numeric(rstandard(kfas.binomial,type="pearson")),
               rstandard(glm.binomial,type="pearson"),tolerance=tol, check.attributes=FALSE)
})


test_that("Residuals for Gamma GLM works properly",{
  expect_equal(as.numeric(residuals(kfas.gamma2,type="pearson")),
               residuals(glm.gamma,type="pearson"),tolerance=tol, check.attributes=FALSE)
  expect_equal(as.numeric(residuals(kfas.gamma2,type="response")),
               residuals(glm.gamma,type="response"),tolerance=tol, check.attributes=FALSE)
  
  expect_equal(as.numeric(rstandard(kfas.gamma2,type="pearson")),
               rstandard(glm.gamma,type="pearson"),tolerance=tol, check.attributes=FALSE)
})


test_that("Residuals for negative binomial GLM works properly",{
  expect_equal(as.numeric(residuals(kfas.NB,type="pearson")),
               residuals(glm.NB,type="pearson"),tolerance=tol, check.attributes=FALSE)
  expect_equal(as.numeric(residuals(kfas.NB,type="response")),
               residuals(glm.NB,type="response"),tolerance=tol, check.attributes=FALSE)
  
  expect_equal(as.numeric(rstandard(kfas.NB,type="pearson")),
               rstandard(glm.NB,type="pearson"),tolerance=tol, check.attributes=FALSE)
})


test_that("Predictions for GLM works properly",{
  pred.glm.gaussian.link<-predict(glm.gaussian,type="link",se.fit=TRUE)
  pred.kfas.gaussian.link<-predict(model.gaussian,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.gaussian.response<-predict(glm.gaussian,type="response",se.fit=TRUE)
  pred.kfas.gaussian.response<-predict(model.gaussian,type="response",se.fit=TRUE,interval="confidence")  
  expect_equal(c(pred.kfas.gaussian.link[,"fit"]),pred.glm.gaussian.link$fit,tolerance=tol, check.attributes=FALSE)
  expect_equal(c(pred.kfas.gaussian.link[,"se.fit"]),pred.glm.gaussian.link$se.fit,tolerance=tol, check.attributes=FALSE)
  expect_equal(c(pred.kfas.gaussian.response[,"fit"]),pred.glm.gaussian.response$fit,tolerance=tol, check.attributes=FALSE)
  expect_equal(c(pred.kfas.gaussian.response[,"se.fit"]),pred.glm.gaussian.response$se.fit,tolerance=tol, check.attributes=FALSE)
  
  
  pred.glm.poisson.link<-predict(glm.poisson,type="link",se.fit=TRUE)
  pred.kfas.poisson.link<-predict(model.poisson,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.poisson.response<-predict(glm.poisson,type="response",se.fit=TRUE)
  pred.kfas.poisson.response<-predict(model.poisson,type="response",se.fit=TRUE,interval="confidence")  
  expect_equal(pred.glm.poisson.link$fit,obj=c(pred.kfas.poisson.link[,"fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.poisson.link$se.fit,obj=c(pred.kfas.poisson.link[,"se.fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.poisson.response$fit,obj=c(pred.kfas.poisson.response[,"fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.poisson.response$se.fit,obj=c(pred.kfas.poisson.response[,"se.fit"]),tolerance=tol, check.attributes=FALSE)
  
  
  pred.glm.binomial.link<-predict(glm.binomial,type="link",se.fit=TRUE)
  pred.kfas.binomial.link<-predict(model.binomial,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.binomial.response<-predict(glm.binomial,type="response",se.fit=TRUE)
  pred.kfas.binomial.response<-predict(model.binomial,type="response",se.fit=TRUE,interval="confidence")  
  expect_equal(pred.glm.binomial.link$fit,obj=c(pred.kfas.binomial.link[,"fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.binomial.link$se.fit,obj=c(pred.kfas.binomial.link[,"se.fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.binomial.response$fit,obj=c(pred.kfas.binomial.response[,"fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.binomial.response$se.fit,obj=c(pred.kfas.binomial.response[,"se.fit"]),tolerance=tol, check.attributes=FALSE)
  
  
  pred.glm.gamma.link<-predict(glm.gamma,type="link",se.fit=TRUE)
  pred.kfas.gamma.link<-predict(model.gamma2,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.gamma.response<-predict(glm.gamma,type="response",se.fit=TRUE)
  pred.kfas.gamma.response<-predict(model.gamma2,type="response",se.fit=TRUE,interval="confidence")  
  expect_equal(pred.glm.gamma.link$fit,obj=c(pred.kfas.gamma.link[,"fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.gamma.link$se.fit,obj=c(pred.kfas.gamma.link[,"se.fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.gamma.response$fit,obj=c(pred.kfas.gamma.response[,"fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.gamma.response$se.fit,obj=c(pred.kfas.gamma.response[,"se.fit"]),tolerance=tol, check.attributes=FALSE)
  
  
  pred.glm.NB.link<-predict(glm.NB,type="link",se.fit=TRUE)
  pred.kfas.NB.link<-predict(model.NB,type="link",se.fit=TRUE,interval="confidence")
  pred.glm.NB.response<-predict(glm.NB,type="response",se.fit=TRUE)
  pred.kfas.NB.response<-predict(model.NB,type="response",se.fit=TRUE,interval="confidence")  
  expect_equal(pred.glm.NB.link$fit,obj=c(pred.kfas.NB.link[,"fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.NB.link$se.fit,obj=c(pred.kfas.NB.link[,"se.fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.NB.response$fit,obj=c(pred.kfas.NB.response[,"fit"]),tolerance=tol, check.attributes=FALSE)
  expect_equal(pred.glm.NB.response$se.fit,obj=c(pred.kfas.NB.response[,"se.fit"]),tolerance=tol, check.attributes=FALSE)
  
  #confidence intervals
  glm.gaussian.CI<-pred.glm.gaussian.link$fit +qnorm(0.025)*pred.glm.gaussian.link$se.fit%o%c(1,-1)
  expect_equal(unclass(pred.kfas.gaussian.link[,c("lwr","upr")]),glm.gaussian.CI,tolerance=tol, check.attributes=FALSE)
  expect_equal(unclass(pred.kfas.gaussian.response[,c("lwr","upr")]),gaussian()$linkinv(glm.gaussian.CI),tolerance=tol, check.attributes=FALSE)
  
  glm.poisson.CI<-pred.glm.poisson.link$fit +qnorm(0.025)*pred.glm.poisson.link$se.fit%o%c(1,-1)
  expect_equal(unclass(pred.kfas.poisson.link[,c("lwr","upr")]),glm.poisson.CI,tolerance=tol, check.attributes=FALSE)
  expect_equal(unclass(pred.kfas.poisson.response[,c("lwr","upr")]),poisson()$linkinv(glm.poisson.CI),tolerance=tol, check.attributes=FALSE)
  
  glm.binomial.CI<-pred.glm.binomial.link$fit +qnorm(0.025)*pred.glm.binomial.link$se.fit%o%c(1,-1)
  expect_equal(unclass(pred.kfas.binomial.link[,c("lwr","upr")]),glm.binomial.CI,tolerance=tol, check.attributes=FALSE)
  expect_equal(unclass(pred.kfas.binomial.response[,c("lwr","upr")]),binomial()$linkinv(glm.binomial.CI),tolerance=tol, check.attributes=FALSE)
  
  glm.gamma.CI<-pred.glm.gamma.link$fit +qnorm(0.025)*pred.glm.gamma.link$se.fit%o%c(1,-1)
  expect_equal(unclass(pred.kfas.gamma.link[,c("lwr","upr")]),glm.gamma.CI,tolerance=tol, check.attributes=FALSE)
  expect_equal(unclass(pred.kfas.gamma.response[,c("lwr","upr")]),Gamma(link="log")$linkinv(glm.gamma.CI),tolerance=tol, check.attributes=FALSE)
  
  glm.NB.CI<-pred.glm.NB.link$fit +qnorm(0.025)*pred.glm.NB.link$se.fit%o%c(1,-1)
  expect_equal(unclass(pred.kfas.NB.link[,c("lwr","upr")]),glm.NB.CI,tolerance=tol, check.attributes=FALSE)
  expect_equal(unclass(pred.kfas.NB.response[,c("lwr","upr")]),glm.NB$family$linkinv(glm.NB.CI),tolerance=tol, check.attributes=FALSE)
})
