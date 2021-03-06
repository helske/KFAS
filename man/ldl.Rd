% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ldl.R
\name{ldl}
\alias{ldl}
\title{LDL Decomposition of a Matrix}
\usage{
ldl(x, tol)
}
\arguments{
\item{x}{Symmetrix matrix.}

\item{tol}{Tolerance parameter for LDL decomposition, determines which
diagonal values are counted as zero. Same value is used in isSymmetric function. 
Default is \code{max(100, max(abs(diag(as.matrix(x))))) * .Machine$double.eps}.}
}
\value{
Transformed matrix with D in diagonal, L in strictly lower diagonal
and zeros on upper diagonal.
}
\description{
Function \code{ldl} computes the LDL decomposition of a positive semidefinite matrix.
}
\examples{
# Positive semidefinite matrix, example matrix taken from ?chol
x <- matrix(c(1:5, (1:5)^2), 5, 2)
x <- cbind(x, x[, 1] + 3*x[, 2])
m <- crossprod(x)
l <- ldl(m, tol = 1e-8) # arm64 Mac setup in CRAN fails with default tolerance
d <- diag(diag(l))
diag(l) <- 1
all.equal(l \%*\% d \%*\% t(l), m, tol = 1e-15)
}
