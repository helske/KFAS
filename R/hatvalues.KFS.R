#' Extract Hat Values from KFS Output
#'
#' Extract hat values from KFS output, when \code{KFS} was run with signal
#' (non-Gaussian case) or mean smoothing (Gaussian case).
#'
#' @details
#'
#' Hat values in \code{KFAS} are defined as the diagonal elements of \code{V_t/H_t} where V_t
#' is the covariance matrix of signal/mean at time t and H_t is the covariance
#' matrix of disturbance vector \eqn{\epsilon} of (approximating) Gaussian model
#' at time t. This definition gives identical results with the standard
#' definition in case of GLMs. Note that it is possible to construct a state
#' space model where this definition is not meaningful (for example the
#' covariance matrix H_t can contain zeros on diagonal).
#'
#' @export
#' @importFrom stats hatvalues
#' @param model An object of class \code{KFS}.
#' @param \dots Additional arguments to \code{approxSSM}.
#' @return Multivariate time series containing hat values.
#' @examples
#' model <- SSModel(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
#' out <- KFS(model, filtering = "state", smoothing = "none")
#' # estimate sigma2
#' model["H"] <- mean(c(out$v[1:out$d][out$Finf==0]^2/out$F[1:out$d][out$Finf==0],
#'                      out$v[-(1:out$d)]^2/out$F[-(1:out$d)]))
#' c(hatvalues(KFS(model)))
#'
hatvalues.KFS <- function(model, ...) {
  if (any(model$model$distribution != "gaussian")) {
    app <- approxSSM(model$model, ...)
    if (is.null(model$V_theta))
      stop("KFS was run without signal smoothing, cannot compute hat values.")
    hatv <- matrix(apply(model$V_theta/app$H, 3, diag), attr(model$model, "n"),
      attr(model$model, "p"), byrow = TRUE)
  } else {
    if (is.null(model$V_mu))
      stop("KFS was run without mean smoothing, cannot compute hat values.")
    hatv <- matrix(apply(model$V_mu, 3, diag),
      attr(model$model, "n"), attr(model$model,"p"), byrow = TRUE)/
      matrix(apply(model$model$H, 3, diag), attr(model$model, "n"),
        attr(model$model, "p"), byrow = TRUE)
  }
  attributes(hatv) <- attributes(model$model$y)
  hatv[is.na(model$model$y)]<-NA
  hatv
}
