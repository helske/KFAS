#' Diagnostic Plots of State Space Models
#'
#' Diagnostic plots based on standardized residuals for objects of class \code{SSModel}.
#'
#' @export
#' @param x Object of class \code{SSModel}.
#' @param nsim The number of independent samples used in importance sampling.
#' Only used for non-Gaussian model. Default is 0, which computes the
#' approximating Gaussian model by \code{\link{approxSSM}} and performs the
#' usual Gaussian filtering/smoothing so that the smoothed state estimates
#' equals to the conditional mode of \eqn{p(\alpha_t|y)}{p(\alpha[t]|y)}.
#' In case of \code{nsim = 0}, the mean estimates and their variances are computed using
#' the Delta method (ignoring the covariance terms).
#' @param plot_args Additional arguments to \code{\link{plot.ts}}
#' @param acf_args Additional arguments to \code{\link{acf}}.
#' @param ... Ignored.
#' @examples
#' modelNile <- SSModel(Nile ~ SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
#' modelNile <- fitSSM(inits = c(log(var(Nile)),log(var(Nile))), model = modelNile,
#'  method = "BFGS")$model
#'
#' if (interactive()) {
#'   plot(modelNile)
#' }

plot.SSModel <- function(x, nsim = 0, plot_args = NULL, acf_args = NULL, ...) {

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  acf2 <- function(x, acf_args) {
    tmp <- acf(unclass(x), plot = FALSE, na.action = na.pass, acf_args)
    tmp$n.used <- sum(!is.na(x))
    plot(tmp)
  }

  par(ask = TRUE)
  if (all(x$distribution == "gaussian")) {
    out <- KFS(x, smoothing = c("mean", "disturbance"), filtering = c("state","mean"))
    res <- rstandard(out)
    if (attr(x, "p") == 1) {
      do.call(plot,
        list(x = cbind(
          observations = x$y,
          recursive = res,
          irregular = rstandard(out, "pearson")),
          main = "recursive (one-step ahead) and irregular (smoothed) residuals", plot_args))
      do.call(acf2,
        list(x = res, acf_args = acf_args))
      qqnorm(as.numeric(res))
      qqline(as.numeric(res))
    } else {

      do.call(plot,
        list(x = x$y,
          main = "Observations", plot_args))
      do.call(plot,
        list(x = res,
          main = "Recursive (one-step ahead) residuals", plot_args))
      do.call(plot,
        list(x = rstandard(out, "pearson"),
          main = "Irregular (smoothed) residuals", plot_args))
      do.call(acf2,
        list(x = res, acf_args = acf_args))
      for (i in 1:attr(x, "p")) {
        qqnorm(res[, i], main = paste("Normal Q-Q plot for",colnames(res)[i]))
        qqline(res[, i])
      }

    }

  } else {
    out <- KFS(x, smoothing = "mean", filtering = "mean", nsim = nsim)
    res <- rstandard(out)

    if (attr(x, "p") == 1) {
      do.call(plot,
        list(x = cbind(
          observations = x$y,
          recursive = res,
          irregular = rstandard(out, "pearson")),
          main = "recursive (one-step ahead) and irregular (smoothed) residuals", plot_args))
      do.call(acf2,
        list(x = res, acf_args = acf_args))
    } else {

      do.call(plot,
        list(x = x$y,
          main = "Observations", plot_args))
      do.call(plot,
        list(x = res,
          main = "Recursive (one-step ahead) residuals", plot_args))
      do.call(plot,
        list(x = rstandard(out, "pearson"),
          main = "Irregular (smoothed) residuals", plot_args))
      do.call(acf2,
        list(x = res, acf_args = acf_args))
    }
  }


}
