#' Deviance of a State Space Model
#' 
#' Returns the deviance of a object of class \code{KFS} defined as the sum of squared deviance residuals.
#' 
#' Note that this is NOT a $-2*(logL-logL*)$ where $L$ is the likelihood and $L*$ is the saturated likelihood.
#' Instead this is based on the conditional likelihood $p(y|theta)$ i.e. it disregards the effect of hidden states.
#' Therefore the value returned by this function might not make sense for non-GLM setting.
#' it is just 
#' Thus this deviance is defined as -2*
#' @export
#' @param object An object of class \code{KFS}.
#' @param \dots Ignored.
#' @return The value of the deviance extracted from object.
deviance.KFS <- 
  function(object, ...) {
    sum(residuals(object, type = "deviance")^2, na.rm = TRUE)
} 
