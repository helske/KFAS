#' KFAS: Functions for Exponential Family State Space Models
#' 
#' Package KFAS contains functions for Kalman filtering, smoothing and 
#' simulation of linear state space models with exact diffuse initialization.
#' 
#' The linear Gaussian state space model is given by
#' 
#' \deqn{y_t = Z_t \alpha_t + \epsilon_t,}{y[t] = Z[t]\alpha[t] + \epsilon[t],
#' (observation equation)}
#' 
#' \deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t,}{\alpha[t+1] = T[t]\alpha[t]
#' + R[t]\eta[t], (transition equation)}
#' 
#' where \eqn{\epsilon_t \sim N(0,H_t)}{\epsilon[t] ~ N(0,H[t])}, \eqn{\eta_t
#' \sim N(0,Q_t)}{\eta[t] ~ N(0,Q[t])} and \eqn{\alpha_1 \sim
#' N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])} independently of each other.
#' 
#' All system and covariance matrices \code{Z}, \code{H}, \code{T}, \code{R} and
#' \code{Q} can be time-varying, and partially or totally missing observations
#' \eqn{y_t}{y[t]} are allowed.
#' 
#' Covariance matrices H and Q has to be positive semidefinite (although this is
#' not checked).
#' 
#' Dimensions of system matrices are
#'
#' \tabular{rl}{
#'  \code{Z} \tab \eqn{p \times m \times 1}{p*m*1} or \eqn{p \times m \times n}{p*m*n} in time varying case \cr
#'  \code{H} \tab \eqn{p \times p \times 1}{p*p*1} or \eqn{p \times p \times n}{p*p*n} in time varying case (Omitted in non-gaussian models) \cr
#'  \code{T} \tab \eqn{m \times m \times 1}{m*m*1} or \eqn{m \times m \times n}{m*m*n} in time varying case \cr
#'  \code{R} \tab \eqn{m \times k \times 1}{m*k*1} or \eqn{m \times k \times n}{m*k*n} in time varying case \cr
#'  \code{Q} \tab \eqn{k \times k \times 1}{k*k*1} or \eqn{k \times k \times n}{k*k*n} in time varying case \cr
#'  \code{u} \tab \eqn{n \times p}{p*n} (Omitted in gaussian models) \cr
#'  }
#' 
#' In case of any of the series in model is defined as non-gaussian, the 
#' observation equation is of form \deqn{\prod_i^p 
#' p_i(y_{i,t}|\theta_t)}{\prod[i]^p p(y[i,t]|\theta[t]),} with 
#' \eqn{\theta_{i,t}=Z_{i,t}\alpha_t}{\theta[i,t]=Z[i,t]\alpha[t]} being one of 
#' the following:
#' 
#' If observations \eqn{y_{i,1},\ldots,y_{i,n}}{y[i,1],\ldots,y[i,n]} are 
#' distributed as \eqn{N(\mu_t,u_t)}{N(\mu[t],u[t])}, then
#' \eqn{\theta_t=\mu_t}{\theta[t]=\mu[t]}. Note that now variances are defined
#' using \code{u}, not \code{H}. If the correlation between Gaussian observation
#' equations is needed, one can use \eqn{u_t=0}{u[t]=0} and add correlating
#' disturbances into state equation (although care is then needed when making
#' inferences about signal which contains the error terms also).
#' 
#' If observations are distributed as
#' \eqn{Poisson(u_t\lambda_t)}{Poisson(u[t]\lambda[t])}, where \eqn{u_t}{u[t]}
#' is an offset term, then \eqn{\theta_t = 
#' log(u_t\lambda_t)}{\theta[t]=log(u[t]\lambda[t])}.
#' 
#' If observations are distributed as
#' \eqn{binomial(u_t,\pi_t)}{binomial(u[t],\pi[t])}, then \eqn{\theta_t =
#' log[\pi_t/(1-\pi_t)]}{\theta[t] = log(\pi[t]/(1-\pi[t]))}, where
#' \eqn{\pi_t}{\pi[t]} is the probability of success at time \eqn{t}.
#' 
#' If observations are distributed as
#' \eqn{gamma(u_t,\mu_t)}{gamma(u[t],\mu[t])}, then \eqn{\theta_t =
#' log(\mu_t)}{[\theta[t] = log(\mu[t])]}, where \eqn{\mu_t}{\mu[t]} is the mean
#' parameter and \eqn{u_t}{u[t]} is the shape parameter.
#' 
#' If observations are distributed as \eqn{negative
#' binomial(u_t,\mu_t)}{negative binomial(u[t],\mu[t])} with expected value
#' \eqn{\mu_t}{\mu[t]} and variance \eqn{\mu_t+ \mu_t^2/u_t}\eqn{\mu[t]+
#' \mu[t]^2/u[t]} (see \code{\link{dbinom}}), then \eqn{\theta_t =
#' log[\mu_t]}{\theta[t] = log(\mu[t])}.
#' 
#' For exponential family models \eqn{u_t=1}{u[t]=1} as a default. For completely Gaussian models, parameter is omitted.
#'
#'
#' For the unknown elements of initial state vector \eqn{a_1}{a[1]}, KFAS uses
#' exact diffuse initialization by Koopman and Durbin (2000, 2001, 2003), where
#' the unknown initial states are set to have a zero mean and infinite variance,
#' so \deqn{P_1 = P_{\ast,1} + \kappa P_{\infty,1},}{P[1] = P[*,1] +
#' \kappaP[inf,1],} with \eqn{\kappa} going to infinity and
#' \eqn{P_{\infty,1}}{P[inf,1]} being diagonal matrix with ones on diagonal
#' elements corresponding to unknown initial states.
#' 
#' This method is basically a equivalent of setting uninformative priors for the
#' initial states in a Bayesian setting. Note that although the states are set 
#' as independent a priori, this effect vanishes quickly.
#' 
#' Diffuse phase is continued until rank of \eqn{P_{\infty,t}}{P[inf,t]} becomes
#' zero. Rank of \eqn{P_{\infty}}{P[inf]} decreases by 1, if 
#' \eqn{F_\infty>tol>0}{F[inf]>tol>0}. Usually the number of diffuse time points
#' equals the number unknown elements of initial state vector, but missing 
#' observations or time-varying system matrices can affect this. See Koopman and
#' Durbin (2000, 2001, 2003) for details for exact diffuse and non-diffuse 
#' filtering.  If the number of diffuse states is large compared to the data, it
#' is possible that the model is degenerate in a sense that not enough
#' information is available for leaving the diffuse phase.
#' 
#' To lessen the notation and storage space, KFAS uses letters P, F and K for
#' non-diffuse part of the corresponding matrices, omitting the asterisk in
#' diffuse phase.
#' 
#' All functions of KFAS use the univariate approach (also known as sequential
#' processing, see Anderson and Moore (1979)) which is from Koopman and Durbin
#' (2000, 2001). In univariate approach the observations are introduced one
#' element at the time. Therefore the prediction error variance matrices F and
#' Finf does not need to be non-singular, as there is no matrix inversions in
#' univariate approach algorithm.  This provides more stable and possibly more
#' faster filtering and smoothing than normal multivariate Kalman filter
#' algorithm. If covariance matrix H is not diagonal, it is possible to
#' transform the model by either using LDL decomposition on H, or augmenting the
#' state vector with \eqn{\epsilon} disturbances. See \code{\link{transformSSM}}
#' for more details.
#' 
#' 
#' @references Koopman, S.J. and Durbin J. (2000).  Fast filtering and
#' smoothing for non-stationary time series models, Journal of American
#' Statistical Assosiation, 92, 1630-38.
#'
#' Koopman, S.J. and Durbin J. (2001).  Time Series Analysis by State Space
#' Methods. Oxford: Oxford University Press.
#'
#' Koopman, S.J. and Durbin J. (2003).  Filtering and smoothing of state vector
#' for diffuse state space models, Journal of Time Series Analysis, Vol. 24,
#' No. 1.
#' 
#' #' Shumway, Robert H. and Stoffer, David S. (2006).  Time Series Analysis and
#' Its Applications: With R examples.  \cr
#' @docType package
#' @name KFAS
#' @aliases KFAS
#' @useDynLib KFAS, .registration=TRUE
#' @examples
#'
#' 
# 
#' # Example of local level model for Nile series
#'
#' modelNile<-SSModel(Nile~SSMtrend(1,Q=list(matrix(NA))),H=matrix(NA))
#' modelNile
#' modelNile<-fitSSM(inits=c(log(var(Nile)),log(var(Nile))),model=modelNile,
#'                   method='BFGS',control=list(REPORT=1,trace=1))$model
#' # Filtering and state smoothing
#' out<-KFS(modelNile,filtering='state',smoothing='state') 
#' out
#'
#' # Confidence and prediction intervals for the expected value and the observations.
#' # Note that predict uses original model object, not the output from KFS.
#' conf<-predict(modelNile,interval='confidence')
#' pred<-predict(modelNile,interval='prediction')
#' 
#' ts.plot(cbind(Nile,pred,conf[,-1]),col=c(1:2,3,3,4,4), 
#'         ylab='Predicted Annual flow', main='River Nile')
#'
#'
#' # Missing observations, using same parameter estimates
#'
#' y<-Nile
#' y[c(21:40,61:80)]<-NA
#' modelNile<-SSModel(y~SSMtrend(1,Q=list(modelNile$Q)),H=modelNile$H)
#'
#' out<-KFS(modelNile,filtering='mean',smoothing='mean')
#'
#' # Filtered and smoothed states
#' plot.ts(cbind(y,fitted(out,filtered=TRUE),fitted(out)), plot.type='single', 
#'         col=1:3, ylab='Predicted Annual flow', main='River Nile')
#' 
#' 
#' # Example of multivariate local level model with only one state
#' # Two series of average global temperature deviations for years 1880-1987
#' # See Shumway and Stoffer (2006), p. 327 for details
#'
#' data(GlobalTemp)
#' 
#' model<-SSModel(GlobalTemp~SSMtrend(1,Q=NA,type='common'),H=matrix(NA,2,2))        
#'
#' # Estimating the variance parameters
#' inits<-chol(cov(GlobalTemp))[c(1,4,3)]
#' inits[1:2]<-log(inits[1:2])
#' fit<-fitSSM(inits=c(0.5*log(.1),inits),model=model,method='BFGS')
#' 
#' out<-KFS(fit$model)
#'
#' ts.plot(cbind(model$y,coef(out)),col=1:3)
#' legend('bottomright',legend=c(colnames(GlobalTemp), 'Smoothed signal'), col=1:3, lty=1)
#'
#'
#'
#' # Seatbelts data
#' \dontrun{
#' model<-SSModel(log(drivers)~SSMtrend(1,Q=list(NA))+
#'                SSMseasonal(period=12,sea.type='trigonometric',Q=NA)+
#'                log(PetrolPrice)+law,data=Seatbelts,H=NA)
#' 
#' # As trigonometric seasonal contains several disturbances which are all 
#' # identically distributed, default behaviour of fitSSM is not enough, 
#' # as we have constrained Q. We can either provide our own 
#' # model updating function with fitSSM, or just use optim directly:
#'
#' # option 1:
#' ownupdatefn<-function(pars,model,...){
#'   model$H[]<-exp(pars[1])
#'   diag(model$Q[,,1])<-exp(c(pars[2],rep(pars[3],11)))
#'   model #for option 2, replace this with -logLik(model) and call optim directly
#' }
#' 
#' fit<-fitSSM(inits=log(c(var(log(Seatbelts[,'drivers'])),0.001,0.0001)),
#'             model=model,updatefn=ownupdatefn,method='BFGS')
#' 
#' out<-KFS(fit$model,smoothing=c('state','mean'))
#' out
#' ts.plot(cbind(out$model$y,fitted(out)),lty=1:2,col=1:2,
#' main='Observations and smoothed signal with and without seasonal component')
#' lines(signal(out,states=c('regression','trend'))$signal,col=4,lty=1)
#' legend('bottomleft',
#' legend=c('Observations', 'Smoothed signal','Smoothed level'), 
#' col=c(1,2,4), lty=c(1,2,1))
#' 
#' 
#' # Multivariate model with constant seasonal pattern,
#' # using the the seat belt law dummy only for the front seat passangers, 
#' # and restricting the rank of the level component by using custom component 
#' 
#' model<-SSModel(log(cbind(front,rear))~ -1 + log(PetrolPrice) + log(kms)
#'                + SSMregression(~law,data=Seatbelts,index=1)
#'                + SSMcustom(Z=diag(2),T=diag(2),R=matrix(1,2,1),
#'                            Q=matrix(1),P1inf=diag(2))
#'                + SSMseasonal(period=12,sea.type='trigonometric'),
#'                  data=Seatbelts,H=matrix(NA,2,2)) 
#' 
#' # An alternative way for defining the rank deficient trend component:
#' \dontrun{#' 
#' model<-SSModel(log(cbind(front,rear))~ -1 + log(PetrolPrice) + log(kms)
#'                + SSMregression(~law,data=Seatbelts,index=1)
#'                + SSMtrend(degree = 1, Q=list(matrix(0,2,2)))
#'                + SSMseasonal(period=12,sea.type='trigonometric'),
#'                  data=Seatbelts,H=matrix(NA,2,2)) 
#' # Modify model manually:                  
#' model$Q<-array(1,c(1,1,1))
#' model$R<-model$R[,-2,,drop=FALSE]
#' attr(model,"k")<-as.integer(1)
#' attr(model,"eta_types")<-attr(model,"eta_types")[1]
#' }
#' 
#' likfn<-function(pars,model,estimate=TRUE){
#'   diag(model$H[,,1])<-exp(0.5*pars[1:2])
#'   model$H[1,2,1]<-model$H[2,1,1]<-tanh(pars[3])*prod(sqrt(exp(0.5*pars[1:2]))) 
#'   model$R[28:29]<-exp(pars[4:5])
#'   if(estimate) return(-logLik(model))
#'   model
#' }  
#' 
#' model<-SSModel(log(cbind(front,rear))~ -1 + log(PetrolPrice) + log(kms)
#'                + SSMregression(~law,data=Seatbelts,index=1)
#'                + SSMtrend(degree = 1, Q=list(matrix(0,2,2)))
#'                + SSMseasonal(period=12,sea.type='trigonometric'),
#'                  data=Seatbelts,H=matrix(NA,2,2)) 
#' 
#' 
#' 
#' likfn<-function(pars,model,estimate=TRUE){
#'   diag(model$H[,,1])<-exp(0.5*pars[1:2])
#'   model$H[1,2,1]<-model$H[2,1,1]<-tanh(pars[3])*prod(sqrt(exp(0.5*pars[1:2]))) 
#'   model$R[28:29]<-exp(pars[4:5])
#'   if(estimate) return(-logLik(model))
#'   model
#' }        
# 
#' fit<-optim(f=likfn,p=c(-7,-7,1,-1,-3),method='BFGS',model=model)
#' model<-likfn(fit$p,model,estimate=FALSE)
#' model$R[28:29,,1]%*%t(model$R[28:29,,1])
#' model$H
#' 
#' out<-KFS(model)
#' out
#' ts.plot(cbind(signal(out,states=c('custom','regression'))$signal,model$y),col=1:4)
#' 
#' # For confidence or prediction intervals, use predict on the original model
#' pred <- predict(model,states=c('custom','regression'),interval='prediction')
#' ts.plot(pred$front,pred$rear,model$y,col=c(1,2,2,3,4,4,5,6),lty=c(1,2,2,1,2,2,1,1))
#' }
#' 
#'  \dontrun{
#' # Poisson model
#' model<-SSModel(VanKilled~law+SSMtrend(1,Q=list(matrix(NA)))+
#'                SSMseasonal(period=12,sea.type='dummy',Q=NA),
#'                data=Seatbelts, distribution='poisson')
#'
#' # Estimate variance parameters
#' fit<-fitSSM(inits=c(-4,-7), model=model,method='BFGS')
#'
#' model<-fit$model
#' 
#' # use approximating model, gives posterior mode of the signal and the linear predictor
#' out_nosim<-KFS(model,smoothing=c('signal','mean'),nsim=0)
#' # State smoothing via importance sampling
#' out_sim<-KFS(model,smoothing=c('signal','mean'),nsim=1000)
#' 
#' out_nosim
#' out_sim
#' }
#' 
#' # Example of generalized linear modelling with KFS
#'
#' # Same example as in ?glm 
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' print(d.AD <- data.frame(treatment, outcome, counts))
#' glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
#' 
#'
#' model<-SSModel(counts ~ outcome + treatment, data=d.AD, 
#'                distribution = 'poisson')
#' 
#' out<-KFS(model)
#' coef(out,start=1,end=1)
#' coef(glm.D93)
#' 
#' summary(glm.D93)$cov.s
#' out$V[,,1]
#' 
#' outnosim<-KFS(model,smoothing=c('state','signal','mean'))
#' set.seed(1)
#' outsim<-KFS(model,smoothing=c('state','signal','mean'),nsim=1000)
#' 
#' 
#' ## linear 
#' # GLM
#' glm.D93$linear.predictor
#' # approximate model, this is the posterior mode of p(theta|y)
#' c(outnosim$thetahat)            
#' # importance sampling on theta,  gives E(theta|y)            
#' c(outsim$thetahat)                         
#' 
#' 
#' 
#' ## predictions on response scale
#' # GLM
#' fitted(glm.D93) 
#' # approximate model with backtransform, equals GLM
#' c(fitted(outnosim))                
#' # importance sampling on exp(theta)       
#' fitted(outsim)                                 
#' 
#' # prediction variances on link scale
#' # GLM
#' as.numeric(predict(glm.D93,type='link',se.fit=TRUE)$se.fit^2)
#' # approx, equals to GLM results
#' c(outnosim$V_theta)                                    
#' # importance sampling on theta        
#' c(outsim$V_theta)                                            
#' 
#' 
#' # prediction variances on response scale
#' # GLM
#' as.numeric(predict(glm.D93,type='response',se.fit=TRUE)$se.fit^2) 
#' # approx, equals to GLM results
#' c(outnosim$V_mu)                             
#' # importance sampling on theta                    
#' c(outsim$V_mu)                                                   
#' 
#' # approxSSM uses modified step-halving for more robust convergence than glm:
#' y <- rep (0:1,c(15,10))
#' suppressWarnings(glm(formula = y ~ 1, family = binomial(link = "logit"), start = 2))
#' model<-SSModel(y~1,dist="binomial")
#' KFS(model,theta=2)
#' 
#' \dontrun{
#' data(sexratio)
#' model<-SSModel(Male~SSMtrend(1,Q=list(NA)),u=sexratio[,'Total'],data=sexratio,
#'                distribution='binomial')
#' fit<-fitSSM(model,inits=-15,method='BFGS',control=list(trace=1,REPORT=1))
#' fit$model$Q #1.107652e-06
#' 
#' # Computing confidence intervals in response scale 
#' # Uses importance sampling on response scale (4000 samples with antithetics)
#' 
#' pred<-predict(fit$model,type='response',interval='conf',nsim=1000) 
#'  
#' ts.plot(cbind(model$y/model$u,pred),col=c(1,2,3,3),lty=c(1,1,2,2))
#' 
#' # Now with sex ratio instead of the probabilities:
#' imp<-importanceSSM(fit$model,nsim=1000,antithetics=TRUE)
#' sexratio.smooth<-numeric(length(model$y))
#' sexratio.ci<-matrix(0,length(model$y),2)
#' w<-imp$w/sum(imp$w)
#' for(i in 1:length(model$y)){
#'  sexr<-exp(imp$sample[i,1,])
#'  sexratio.smooth[i]<-sum(sexr*w)
#'  oo<-order(sexr)
#'  sexratio.ci[i,]<-c(sexr[oo][which.min(abs(cumsum(w[oo]) - 0.05))],
#'                       sexr[oo][which.min(abs(cumsum(w[oo]) - 0.95))])
#' }
#' 
#' # Same by direct transformation:
#' out<-KFS(fit$model,smoothing='signal',nsim=1000)
#' sexratio.smooth2 <- exp(out$thetahat)
#' sexratio.ci2<-exp(c(out$thetahat) + qnorm(0.025) * sqrt(drop(out$V_theta))%o%c(1, -1))
#' 
#' ts.plot(cbind(sexratio.smooth,sexratio.ci,sexratio.smooth2,sexratio.ci2),
#'         col=c(1,1,1,2,2,2),lty=c(1,2,2,1,2,2))
#'}
#' # Example of Cubic spline smoothing
#' \dontrun{
#' require(MASS)
#' data(mcycle)
#' 
#' model<-SSModel(accel~-1+SSMcustom(Z=matrix(c(1,0),1,2),
#'                                  T=array(diag(2),c(2,2,nrow(mcycle))),
#'                                  Q=array(0,c(2,2,nrow(mcycle))),
#'                                  P1inf=diag(2),P1=diag(0,2)),data=mcycle)
#' 
#' model$T[1,2,]<-c(diff(mcycle$times),1)
#' model$Q[1,1,]<-c(diff(mcycle$times),1)^3/3
#' model$Q[1,2,]<-model$Q[2,1,]<-c(diff(mcycle$times),1)^2/2
#' model$Q[2,2,]<-c(diff(mcycle$times),1)
#' 
#' 
#' updatefn<-function(pars,model,...){ 
#'   model$H[]<-exp(pars[1])
#'   model$Q[]<-model$Q[]*exp(pars[2])
#'   model
#' }
#' 
#' fit<-fitSSM(model,inits=c(4,4),updatefn=updatefn,method='BFGS')
#'
#' pred<-predict(fit$model,interval='conf',level=0.95)
#' plot(x=mcycle$times,y=mcycle$accel,pch=19)
#' lines(x=mcycle$times,y=pred[,1])
#' lines(x=mcycle$times,y=pred[,2],lty=2)
#' lines(x=mcycle$times,y=pred[,3],lty=2)
#' }
#' 
#' 
NULL
#' Oxford-Cambridge boat race results 1829-2000
#'
#' Results of the annual boat race between universities of Oxford (0) and Cambridge (1).
#'
#' @name boat
#' @docType data
#' @format A time series object containing 172 observations.
#' @references  Koopman, S.J. and Durbin J. (2001).  Time Series Analysis by State Space Methods. Oxford: Oxford University Press.
#' @source http://www.ssfpack.com/DKbook.html
#' @keywords datasets
NULL
#' Two series of average global temperature deviations for years 1880-1987
#'
#' This data set contains two series of average global temperature deviations
#' for years 1880-1987. These series are same as used in Shumway and Stoffer
#' (2006), where they are known as HL and Folland series. For more details, see
#' Shumway and Stoffer (2006, p. 327).
#'
#'
#' @name GlobalTemp
#' @docType data
#' @format A time series object containing 108 times 2 observations.
#' @references Shumway, Robert H. and Stoffer, David S. (2006). Time Series
#' Analysis and Its Applications: With R examples.
#' @source http://lib.stat.cmu.edu/general/stoffer/tsa2/
#' @keywords datasets
NULL
#' Number of males and females born in Finland from 1751 to 2011
#'
#' A time series object containing the number of males and females born in Finland from 1751 to 2011.
#' 
#' @name sexratio
#' @docType data
#' @format A time series object containing the number of males and females born in Finland from 1751 to 2011.
#' @source Statistics Finland
#' @keywords datasets
NULL 
