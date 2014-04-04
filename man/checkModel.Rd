\name{is.SSModel}
\alias{is.SSModel}
\title{Test whether object is a valid \code{SSModel} object}
\usage{
is.SSModel(object, na.check = FALSE, return.logical = TRUE)
}
\arguments{
  \item{object}{An object to be tested.}

  \item{na.check}{Test the system matrices for NA and
  infinite values. Default is \code{FALSE}.}

  \item{return.logical}{If \code{FALSE}, error is given if
  the the model is not a valid \code{SSModel} object.
  Otherwise logical value is returned. Defaults to
  \code{FALSE}.}
}
\value{
Logical value or nothing, depending on the value of
\code{return.logical}.
}
\description{
Function \code{is.SSModel} tests whether the object is a
valid \code{SSModel} object.
}
\details{
Note that the validity of the values in \code{y} and
\code{Z} are not tested. These can contain NA values (but
not infinite values), with condition that when
\code{Z[i,,t]} contains NA value, the corresponding
\code{y[t,i]} must also have NA value. In this case
\code{Z[i,,t]} is not referenced in filtering and
smoothing, and algorithms works properly. Note also that
this does result NA values in \code{thetahat}, so it could
be beneficial to use for example zeroes in place of NA
values in Z, making first sure that the above condition is
met.
}

