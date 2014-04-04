\name{fitSSM}
\alias{fitSSM}
\title{Maximum Likelihood Estimation of a State Space Model}
\usage{
fitSSM(model, inits, updatefn, checkfn, ...)
}
\arguments{
  \item{inits}{Initial values for \code{optim}}

  \item{model}{Model object of class \code{SSModel}.}

  \item{updatefn}{User defined function which updates the
  model given the parameters. Must be of form
  \code{updatefn(pars, model,...)}, i.e. must contain
  ellipsis \code{...}. If not supplied, a default function
  is used, which estimates the values marked as NA in time
  invariant covariance matrices Q and H.}

  \item{checkfn}{Optional function for model checking. If
  supplied, after updating the model, if
  \code{checkfn(model)} returns TRUE, -log-likelihood is
  computed, otherwise \code{.Machine$double.xmax} is
  returned. See examples. If not supplied, check.model=TRUE
  is used for checking possible NA or Inf values, see
  ?logLik.SSModel.}

  \item{...}{Further arguments for functions \code{optim},
  \code{updatefn} and \code{logLik.SSModel}, such as
  \code{method='BFGS'}.}
}
\value{
A list with elements \item{optim.out}{Output from function
\code{optim}. } \item{model}{Model with estimated
parameters. }
}
\description{
Function \code{fitSSM} finds the maximum likelihood
estimates for unknown parameters of an arbitary state space
model, given the user-defined model updating function.
}
\examples{
# Example function for updating covariance matrices H and Q
# (also used as a default function in fitSSM)

updatefn <- function(pars,model,...){
Q<-as.matrix(model$Q[,,1])
naQd  <- which(is.na(diag(Q)))
naQnd <- which(upper.tri(Q[naQd,naQd]) & is.na(Q[naQd,naQd]))
Q[naQd,naQd][lower.tri(Q[naQd,naQd])] <- 0
diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
Q[naQd,naQd][naQnd] <- pars[length(naQd)+1:length(naQnd)]
model$Q[naQd,naQd,1] <- crossprod(Q[naQd,naQd])
if(!identical(model$H,'Omitted')){
   H<-as.matrix(model$H[,,1])
   naHd  <- which(is.na(diag(H)))
   naHnd <- which(upper.tri(H[naHd,naHd]) & is.na(H[naHd,naHd]))
   H[naHd,naHd][lower.tri(H[naHd,naHd])] <- 0
   diag(H)[naHd] <- exp(0.5 * pars[length(naQd)+length(naQnd)+1:length(naHd)])
   H[naHd,naHd][naHnd] <- pars[length(naQd)+length(naQnd)+length(naHd)+1:length(naHnd)]
   model$H[naHd,naHd,1] <- crossprod(H[naHd,naHd])
 }

 model
}

# Example function for checking the validity of covariance matrices.

checkfn <- function(model){
  #test positive semidefiniteness of H and Q
  inherits(try(ldl(model$H[,,1]),TRUE),'try-error') ||
  inherits(try(ldl(model$Q[,,1]),TRUE),'try-error')
}
}

