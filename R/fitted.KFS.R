#' Smoothed Estimates or One-step-ahead Predictions of Fitted Values
#'
#' Computes fitted values from output of \code{KFS}
#' (or using the \code{SSModel} object), i.e. one-step-ahead
#' predictions  \eqn{f(\theta_t | y_{t-1}, \ldots, y_1)}{
#' f(\theta[t] | y[t-1], ... , y[1]),} (\code{m}) or smoothed estimates
#' \eqn{f(\theta_t | y_n, \ldots, y_1)}{f(\theta[t] | y[n], ... , y[1]),} (\code{muhat}),
#' where \eqn{f} is the inverse of the link function
#' (identity in Gaussian case), except in case of Poisson distribution where
#' \eqn{f} is multiplied with the exposure \eqn{u_t}{u[t]}.
#'
#' @export
#' @importFrom stats fitted
#' @name fitted.SSModel
#' @inheritParams coef.SSModel
#' @return Multivariate time series containing fitted values.
#' @seealso \code{\link{signal}} for partial signals and their covariances.
#' @examples
#' data("sexratio")
#' model <- SSModel(Male ~ SSMtrend(1,Q = list(NA)),u = sexratio[, "Total"],
#'   data = sexratio, distribution = "binomial")
#' model <- fitSSM(model,inits = -15, method = "BFGS")$model
#' out <- KFS(model)
#' identical(drop(out$muhat), fitted(out))
#'
#' fitted(model)
fitted.KFS <- function(object, start = NULL, end = NULL, filtered = FALSE, ...) {
  if (!filtered) {
    if (!is.null(object$muhat)) {
      tmp <- object$muhat
    } else stop("Input does not contain smoothed estimates for means, rerun KFS with mean smoothing.")
  } else {
    if (!is.null(object[["m", exact = TRUE]])) {
      tmp <- object$m
    } else stop("Input does not contain filtered estimates for means, rerun KFS with mean filtering.")
  }
  tmp <- window(tmp, start = start, end = end)
  if (!is.null(start) && identical(start, end)) {
    tmp[1, ]
  } else {
    drop(tmp)
  }
}
#' @export
#' @rdname fitted.SSModel
fitted.SSModel <- function(object, start = NULL, end = NULL, filtered = FALSE, nsim = 0, ...) {
  if (filtered) {
    out <- KFS(object, filtering = "mean", smoothing = "none", nsim = nsim, ...)
  } else {
    out <- KFS(object, filtering = "none", smoothing = "mean", nsim = nsim, ...)
  }
  fitted(out, start, end, filtered)
}
