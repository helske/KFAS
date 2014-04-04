\name{residuals.KFS}
\alias{residuals.KFS}
\title{Extract Residuals of KFS output}
\usage{
\method{residuals}{KFS}(object, type = c("recursive", "deviance", "pearson",
  "response", "state"), ...)
}
\arguments{
  \item{object}{KFS object}

  \item{type}{Character string defining the type of
  residuals.}

  \item{...}{Ignored.}
}
\description{
Extract Residuals of KFS output
}
\details{
For object of class KFS, several types of residuals can be
computed:

\itemize{ \item 'recursive': One-step ahead prediction
residuals \deqn{v_{t,i}),} with residuals being undefined
in diffuse phase. Only supported for fully Gaussian models.

\item 'response': Data minus fitted values,
\eqn{y-E(y)}{y-E(y)}.

\item 'pearson':
\deqn{(y_{t,i}-\theta_{t,i})/\sqrt{V(\mu)_{t,i}}, \quad
i=1,\ldots,p,t=1,\ldots,n,}{(y[t,i]-\theta[t,i])V(\mu)[t,i]^(-0.5),
i=1,\ldots,p, t=1,\ldots,n,} where
\eqn{V(\mu_{t,i})}{V(\mu[t,i])} is the variance function of
the model.

\item 'state': Residuals based on the smoothed disturbance
terms \eqn{\eta} are defined as \deqn{\hat \eta_t, \quad
t=1,\ldots,n,}{L^{-1}[t] \eta[t], t=1,\ldots,n}.

\item 'deviance': Deviance residuals. }
}

