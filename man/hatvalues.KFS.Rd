\name{hatvalues.KFS}
\alias{hatvalues.KFS}
\title{Extract Hat Values from KFS Output}
\usage{
\method{hatvalues}{KFS}(model, ...)
}
\arguments{
  \item{model}{An object of class \code{KFS}.}

  \item{\dots}{Ignored.}
}
\value{
Multivariate time series containing hat values.
}
\description{
Extract hat values from KFS output, when \code{KFS} was run
with signal (non-Gaussian case) or mean smoothing (Gaussian
case).
}
\details{
Hat values are the diagonal elements of \code{V_t/H_t}
where V_t is the covariance matrix of signal/mean at time t
and H_t is the covariance matrix of disturbance vector
\eqn{\epsilon} of (approximating) Gaussian model at time t.
}

