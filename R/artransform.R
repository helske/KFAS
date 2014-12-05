#' Mapping real valued parameters to stationary region
#'
#' Function \code{artransform} transforms \eqn{p} real valued parameters to 
#' stationary region of \eqn{p}th order autoregressive process using 
#' parametrization suggested by Jones (1980). Fortran code is a converted from
#' stats package's C-function partrans.
#'
#' @export
#' @param param Real valued parameters for the transformation.
#' @return transformed The parameters satisfying the stationary constrains.
artransform <- 
  function(param) {
    param <- tanh(param)
    p <- length(param)
    .Fortran("fartransform", as.double(param), phi = param, 
             as.integer(p))$phi
} 
