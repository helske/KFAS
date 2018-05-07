context("Testing ARIMA models")

test_that("artransform works properly",{
x <- seq(-5, 5, by = 0.01)
expect_equal(Vectorize(artransform)(x), tanh(x))
expect_equal(artransform(c(100, -100)), c(2, -1))
expect_equal(artransform(c(0, 0)), c(0, 0))
expect_equal(artransform(c(1, 1, 1)), c(-0.3984572, 0.6233126, 0.7615942),
             tolerance = 1e-6)
})

test_that("arimaSSM works properly",{
  tol <- 1e-3
  s <- 12
  phis <- 0.99
  phi1 <- 0.0001
  phi <- c(phi1,rep(0,s-2),phis,-phi1*phis)
  theta <- 0.7
  out <- makeARIMA(phi,theta,NULL,SSinit="Ross")
  out2<-SSMarima(phi,theta)
  expect_equal(out$Pn,out2$P1,tolerance=tol,check.attributes=FALSE)
  expect_equal(out$T,out2$T,tolerance=tol,check.attributes=FALSE)
  expect_equal(out$Z,c(out2$Z),tolerance=tol,check.attributes=FALSE)
})
test_that("AR modelling works properly",{
  tol <- 1e-5
  # remove intercept from model as it is estimated differently in arima
  model <- SSModel(lh~-1+SSMarima(ar=c(0,0,0),Q=1),H=0)

  likfn <- function(pars, model, estimate=TRUE){
    pars[1:3] <- artransform(pars[1:3])
    tmp <- try(SSMarima(pars[1:3], Q = exp(pars[4])),silent=TRUE)
    if(!inherits(tmp,"try-error")){
      model["T","arima"] <- tmp$T
      model["Q","arima"] <- tmp$Q
      model["P1","arima"] <- tmp$P1
      if(estimate){
        -logLik(model)
      } else model
    } else {
      if(estimate){
        1e100
      } else model
    }
  }

  fit_kfas <- optim(par = c(rep(0, 3), log(var(lh,na.rm = TRUE))), model = model,
                    fn = likfn, method = "BFGS")
  model <- likfn(fit_kfas$par,model, FALSE)
  fit_arima <- arima(lh, order = c(3, 0, 0), include = FALSE)

  expect_equal(fit_arima$coef,model$T[,1,1],check.attributes=FALSE,tol=tol)
  expect_equal(fit_arima$sigma,model$Q[1],check.attributes=FALSE,tol=tol)
})

test_that("ARIMA modelling works properly",{
  tol <- 1e-4
  data(Nile)
  model <- SSModel(Nile ~ SSMarima(ar=c(0.1,0.2),ma=c(0.1),d=1), H=0)

  #due to stationarity checks it's easier to use own objective function and optim directly
  likfn <- function(pars, model, estimate=TRUE){
    # use artransform to transform parameters so the model is stationary and invertible
    tmp <- try(SSMarima(artransform(pars[1:2]),
                        artransform(pars[3]),d=1,Q = exp(pars[4])),silent=TRUE)
    if(!inherits(tmp,"try-error")){
      model["T","arima"] <- tmp$T
      model["R","arima"] <- tmp$R
      model["P1","arima"] <- tmp$P1
      model["Q","arima"] <- tmp$Q
      if(estimate){
        -logLik(model)
      } else model
    } else {
      if(estimate){
        1e100
      } else model
    }
  }

  inits <- c(0.1,0.5,0.5,log(15000))
  fit <- optim(inits, likfn, model=model, method='BFGS')
  model <- likfn(fit$par,model,FALSE)

  fit_arima <- arima(x=Nile, order = c(2,1,1))
  expect_equal(fit_arima$coef,c(model$T[2:3,2,1],model$R[3]),check.attributes=FALSE,tol=tol)
  expect_equal(fit_arima$sigma,model$Q[1],check.attributes=FALSE,tol=tol)
})
