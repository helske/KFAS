\name{[<-.SSModel}
\alias{[.SSModel}
\alias{[<-.SSModel}
\alias{subset.SSModel}
\alias{subset<-}
\alias{subset<-.SSModel}
\title{Extract or Replace Parts of a State Space Model}
\usage{
\method{[}{SSModel}(x, element, states, etas, series, times) <- value

\method{[}{SSModel}(x, element, states, etas, series, times)

\method{subset}{SSModel}(x, element, states, etas, series, times) <- value

subset(x) <- value

\method{subset}{SSModel}(x, element, states, etas, series, times, ...)
}
\arguments{
  \item{x}{Object of class \code{SSModel}.}

  \item{element}{Which element is chosen. Possible choices
  are 'y','Z','H','T','R','Q','a1','P1','P1inf', and 'u'.}

  \item{states}{Which states are chosen. Either a numeric
  vector containing the indices of the corresponding
  states, or a character vector defining the types of the
  corresponding states. Possible choices are \dQuote{all},
  \dQuote{arima}, \dQuote{custom}, \dQuote{cycle},
  \dQuote{seasonal}, \dQuote{trend}, or
  \dQuote{regression}. These can be combined. Default is
  \dQuote{all}.}

  \item{etas}{Which disturbances eta are chosen. Used for
  elements \dQuote{R} and \dQuote{Q}.Either a numeric
  vector containing the indices of the corresponding etas,
  or a character vector defining the types of the
  corresponding etas. Possible choices are \dQuote{all},
  \dQuote{arima}, \dQuote{custom}, \dQuote{cycle},
  \dQuote{seasonal}, \dQuote{trend}, or
  \dQuote{regression}. These can be combined.}

  \item{series}{Numeric. Which series are chosen. Used for
  elements \dQuote{y}, \dQuote{Z}, and \dQuote{u}.}

  \item{times}{Numeric. Which time points are chosen.}

  \item{value}{A value to be assigned to x.}

  \item{...}{ignored.}
}
\value{
A selected subset of the chosen element or a value.
}
\description{
S3 methods for extracting or replacing parts of objects of
class \code{SSModel}.
}

