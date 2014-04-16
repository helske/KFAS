\name{predict.SSModel}
\alias{predict}
\alias{predict.SSModel}
\title{State Space Model Predictions}
\usage{
\method{predict}{SSModel}(object, newdata, interval = c("none", "confidence",
  "prediction"), level = 0.95, type = c("response", "link"),
  states = NULL, se.fit = FALSE, nsim = 0, prob = TRUE, maxiter = 50,
  ...)
}
\arguments{
  \item{object}{Object of class \code{SSModel}.}

  \item{newdata}{A compatible \code{SSModel} object to be
  added in the end of the old object for which the
  predictions are required.}

  \item{interval}{Type of interval calculation.}

  \item{level}{Confidence level for intervals.}

  \item{type}{Scale of the prediction, \code{'response'} or
  \code{'link'}.}

  \item{states}{Which states are used in computing the
  predictions. Either a numeric vector containing the
  indices of the corresponding states, or a character
  vector defining the types of the corresponding states.
  Possible choices are \dQuote{all}, \dQuote{arima},
  \dQuote{custom}, \dQuote{cycle}, \dQuote{seasonal},
  \dQuote{trend}, or \dQuote{regression}. These can be
  combined. Default is \dQuote{all}.}

  \item{nsim}{Number of independent samples used in
  importance sampling. Used only for non-Gaussian models.}

  \item{se.fit}{If TRUE, standard errors are computed.
  Default is FALSE.}

  \item{prob}{if TRUE (default), the predictions in
  binomial/poisson case are probabilities/rates instead of
  counts.}

  \item{\dots}{Ignored.}
}
\value{
A matrix or list of matrices containing the predictions,
and optionally standard errors.
}
\description{
Function \code{predict.SSModel} predicts the future
observations of a state space model of class
\code{\link{SSModel}}
}
\details{
For non-Gaussian models, the results depend whether
importance sampling is used (\code{nsim>0}). without
simulations, the confidence intervals in response scale are
computed in linear predictor scale, and then transformed to
response scale. The prediction intervals are not supported.
With importance sampling, the confidence intervals are
computed as the empirical quantiles from the weighted
sample, whereas the prediction intervals contain additional
step of simulating the response variables from the sampling
distribution \eqn{p(y|\theta^i)}.

If no simulations are used, the standard errors in response
scale are computed using delta method.
}

