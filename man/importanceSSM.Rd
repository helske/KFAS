\name{importanceSSM}
\alias{importanceSSM}
\title{Importance Sampling of Exponential Family State Space Model}
\usage{
importanceSSM(model, type = c("states", "signals"), filtered = FALSE,
  nsim = 1000, save.model = FALSE, theta, antithetics = FALSE,
  maxiter = 50)
}
\arguments{
  \item{model}{Exponential family state space model of
  class \code{SSModel}.}

  \item{type}{What to simulate, \code{'states'} or
  \code{'signals'}. Default is \code{'states'}}

  \item{filtered}{Simulate from
  \eqn{p(\alpha_t|y_{t-1},...,y_1)} instead of
  \eqn{p(\alpha|y)}.}

  \item{nsim}{Number of independent samples. Default is
  1000.}

  \item{save.model}{Return the original model with the
  samples. Default is FALSE.}

  \item{theta}{Initial values for conditional mode theta.}

  \item{antithetics}{Logical. If TRUE, two antithetic
  variables are used in simulations, one for location and
  another for scale. Default is FALSE.}

  \item{maxiter}{Maximum number of iterations used in
  linearisation. Default is 50.}
}
\value{
A list containing elements \code{samples}, \code{weights}
and \code{model} (if \code{save.model==TRUE}).
}
\description{
Importance Sampling of Exponential Family State Space
Model.
}
\details{
Function \code{importanceSSM} simulates states or signals
of the exponential family state space model conditioned
with the observations, returning the simulated samples of
the states/signals with the corresponding importance
weights.

Function can use two antithetic variables, one for location
and other for scale, so output contains four blocks of
simulated values which correlate which each other (ith
block correlates negatively with (i+1)th block, and
positively with (i+2)th block etc.).
}

