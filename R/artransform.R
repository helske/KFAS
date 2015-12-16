#' Mapping real valued parameters to stationary region
#'
#' Function \code{artransform} transforms \eqn{p} real valued parameters to
#' stationary region of \eqn{p}th order autoregressive process using
#' parametrization suggested by Jones (1980). Fortran code is a converted from
#' \code{stats} package's C-function \code{partrans}.
#'
#' @export
#' @param param Real valued parameters for the transformation.
#' @return transformed The parameters satisfying the stationary constrains.
#' @references Jones, R. H (1980). Maximum likelihood fitting
#' of ARMA models to time series with missing observations, Technometrics
#' Vol 22. p. 389--395.
#' @examples
#' artransform(1:3)
artransform <- function(param) {
  .Fortran("fartransform", phi = tanh(param), length(param))$phi
}
