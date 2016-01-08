#' Defunct Functions of Package KFAS
#'
#' The function listed here are removed from KFAS.
#'
#' Subset-based methods were removed as they were somewhat confusingly
#' named as the \code{subset} generic in \code{base}, and most likely not
#' that useful (compared to \code{\link{[<-.SSModel}}).
#'
#' Deviance.KFS was removed as it was mostly useless. The value was not a \eqn{-2*(logL-logL*)}
#' where \eqn{L} is the likelihood and \eqn{L*} is the saturated likelihood.
#' Instead it was based on the conditional likelihood \eqn{p(y|theta)} i.e. it disregards
#' the effect of hidden states. Therefore the value
#' returned by this function did not make much sense in non-GLM setting.
#'
#' @export
#' @aliases KFAS-defunct
#' @name KFAS-defunct
#' @rdname KFAS-defunct
#' @keywords internal
deviance.KFS <- function(object, ...) {
  .Defunct()
  sum(residuals(object, type = "deviance") ^ 2, na.rm = TRUE)
}

#' @method subset<- SSModel
#' @export
#' @rdname KFAS-defunct
`subset<-.SSModel` <- function(x, element, states, etas, series, times,
                               ..., value) {

  .Defunct(new = "[<-.SSModel")
  x[element, states, etas, series, times, ...] <- value
  x
}
#' @rdname KFAS-defunct
#' @export
`subset<-` <- function(x, ..., value) {
  .Defunct(new = "[<-")
  UseMethod("subset<-")
}

#' @method subset SSModel
#' @export
#' @rdname KFAS-defunct
subset.SSModel <- function(x, element, states, etas, series, times, ...) {
  .Defunct(new = "[.SSModel")
  x[element, states, etas, series, times, ...]

}
