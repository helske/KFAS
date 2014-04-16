\name{logLik.SSModel}
\alias{logLik}
\alias{logLik.SSModel}
\title{Log-likelihood of the State Space Model.}
\usage{
\method{logLik}{SSModel}(object, nsim = 0, antithetics = TRUE, theta,
  check.model = FALSE, transform = c("ldl", "augment"), maxiter = 50,
  seed, convtol = 1e-15, ...)
}
\arguments{
  \item{object}{State space model of class \code{SSModel}.}

  \item{nsim}{Number of independent samples used in
  estimating the log-likelihood of the non-Gaussian state
  space model. Default is 0, which gives good starting
  value for optimization. Only used for non-Gaussian
  model.}

  \item{antithetics}{Logical. If TRUE, two antithetic
  variables are used in simulations, one for location and
  another for scale. Default is TRUE. Only used for
  non-Gaussian model.}

  \item{theta}{Initial values for conditional mode theta.
  Only used for non-Gaussian model.}

  \item{check.model}{Logical. If TRUE, function
  \code{is.SSModel} is called before computing the
  likelihood. Default is FALSE.}

  \item{transform}{How to transform the model in case of
  non-diagonal covariance matrix \eqn{H}. Defaults to
  \code{'ldl'}. See function \code{\link{transformSSM}} for
  details.}

  \item{maxiter}{The maximum number of iterations used in
  linearisation. Default is 50. Only used for non-Gaussian
  model.}

  \item{seed}{The value is used as a seed via set.seed
  function. Only used for non-Gaussian model.}

  \item{convtol}{Tolerance parameter for convergence checks
  for Gaussian approximation.  Iterations are continued
  until
  \eqn{tol>abs(dev_{old}-dev_{new})/(abs(dev_{new})+0.1))}.}

  \item{...}{Ignored.}
}
\value{
\item{}{log-likelihood of the state space model.}
}
\description{
Function \code{logLik.SSmodel} computes the log-likelihood
value of a state space model.
}

