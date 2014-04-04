\name{signal}
\alias{signal}
\title{Extracting the Partial Signal Of a State Space Model}
\usage{
signal(object, states = "all", filtered = FALSE)
}
\arguments{
  \item{object}{Object of class \code{KFS}.}

  \item{states}{Which states are combined? Either a numeric
  vector containing the indices of the corresponding
  states, or a character vector defining the types of the
  corresponding states. Possible choices are \dQuote{all},
  \dQuote{arima}, \dQuote{custom}, \dQuote{cycle},
  \dQuote{seasonal}, \dQuote{trend}, or
  \dQuote{regression}. These can be combined. Default is
  \dQuote{all}.}

  \item{filtered}{If TRUE, filtered signal is used.
  Otherwise smoothed signal is used.}
}
\value{
\item{signal}{Time series object of filtered signal
\eqn{Z_ta_t}{Z[t]a[t]} or smoothed signal
\eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]} using only the defined
states.  } \item{variance}{Cov(\eqn{Z_ta_t}{Z[t]a[t]}) or
Cov(\eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]}) using only the
defined states. For the covariance matrices of the filtered
signal, only the non-diffuse part P is used.  }
}
\description{
Function \code{signal} returns the signal of a state space
model using only subset of states.
}

