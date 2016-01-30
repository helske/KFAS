#' Extract Residuals of KFS output
#'
#' @details For object of class KFS, several types of residuals can be computed:
#'
#'   \itemize{ \item  \code{"recursive"}: One-step-ahead prediction residuals
#'   \eqn{v_{t,i}=y_{t,i}-Z_{t,i}a_{t,i}}{v[t,i]=y[t,i]-Z[t,i]a[t,i]}. For non-Gaussian case recursive
#'   residuals are computed as \eqn{y_{t}-f(Z_{t}a_{t})}{y[t]-Z[t]a[t]}, i.e.
#'   non-sequentially. Computing recursive residuals for large non-Gaussian
#'   models can be time consuming as filtering is needed.
#'
#'   \item \code{"pearson"}:  \deqn{(y_{t,i}-\theta_{t,i})/\sqrt{V(\mu_{t,i})},
#'   \quad i=1,\ldots,p,t=1,\ldots,n,}{(y[t,i]-\theta[t,i])V(\mu[t,i])^(-0.5),
#'   i=1,\ldots,p, t=1,\ldots,n,} where \eqn{V(\mu_{t,i})}{V(\mu[t,i])} is the
#'   variance function of the series \eqn{i}
#'
#'   \item \code{"response"}: Data minus fitted values, \eqn{y-E(y)}{y-E(y)}.
#'
#'   \item \code{"state"}:  Residuals based on the smoothed disturbance terms
#'   \eqn{\eta} are defined as \deqn{\hat \eta_t, \quad t=1,\ldots,n.}{\hat \eta[t], t=1,\ldots,n.}
#'   Only defined for fully Gaussian models.
#'
#'   \item \code{"deviance"}: Deviance residuals. Deprecated. This option was
#'   meant to be used only for the GLM comparisons, as their generalization to
#'   other models is lacking, but these will be completely removed in future in
#'   order to avoid misleading results in non-GLM settings. }
#' @export
#' @importFrom stats residuals
#' @param object KFS object
#' @param type Character string defining the type of residuals.
#' @param ... Ignored.
residuals.KFS <-  function(object,
  type = c("recursive", "pearson", "response", "state", "deviance"), ...) {

  type <- match.arg(type)

  if(type == "deviance")
    .Deprecated(msg="Argument type=\"deviance\" is deprecated.")


  if ((type == "state") && any(object$model$distribution !=  "gaussian"))
    stop("State residuals are only supported for fully Gaussian models.")

  series <- switch(type,
    recursive = {
      if(any(object$model$distribution !=  "gaussian") && is.null(object[["m", exact = TRUE]]))
        stop("KFS object does not contain filtered means. ")
      if (all(object$model$distribution ==  "gaussian") && is.null(object[["v", exact = TRUE]]))
        stop("KFS object does not contain prediction errors v. ")
      if(all(object$model$distribution ==  "gaussian")){
        series <- object$v
      }else {
        if (sum(bins <- object$model$distribution == "binomial") > 0)
          series[, bins] <- series[, bins]/object$model$u[, bins]
        series <- object$model$y-object[["m", exact = TRUE]]
      }
      series
    },
    response = {
      series <- object$model$y
      if (sum(bins <- object$model$distribution == "binomial") > 0){
        series[, bins] <- series[, bins]/object$model$u[, bins]
      }
      series - fitted(object)
    },
    state = {
      if (is.null(object$etahat)) {
        stop("KFS object needs to contain smoothed estimates of state disturbances eta.")
      } else {
        object$etahat
      }
    },
    pearson = {
      varianceFunction <- function(object) {
        vars <- object$model$y
        for (i in 1:length(object$model$distribution)){
          vars[, i] <- switch(object$model$distribution[i],
            gaussian = 1,
            poisson = object$muhat[, i],
            binomial = object$muhat[, i] * (1 - object$muhat[, i])/object$model$u[, i],
            gamma = object$muhat[, i]^2,
            `negative binomial` = object$muhat[, i] + object$muhat[, i]^2/object$model$u[, i])
        }
        vars
      }
      series <- object$model$y
      if (sum(bins <- object$model$distribution == "binomial") > 0)
        series[, bins] <- series[, bins]/object$model$u[, bins]
      (series - fitted(object))/sqrt(varianceFunction(object))
    },
    deviance = {
      series <- object$model$y
      if (sum(bins <- object$model$distribution == "binomial") > 0)
        series[, bins] <- series[,bins]/object$model$u[, bins]
      for (i in 1:attr(object$model, "p"))
        series[, i] <-
        ifelse(series[, i] > object$muhat[, i], 1, -1) *
        sqrt(switch(object$model$distribution[i],
          gaussian = (series[, i] - object$muhat[, i])^2,
          poisson = 2 * (series[, i] * log(ifelse(series[, i] == 0, 1, series[, i]/object$muhat[, i])) - series[, i] + object$muhat[, i]),
          binomial = 2 * object$model$u[, i] *
            (series[, i] * log(ifelse(series[, i] == 0, 1, series[, i]/object$muhat[, i])) +
                (1 - series[, i]) * log(ifelse(series[, i] == 1 | object$muhat[, i] == 1, 1, (1 - series[, i])/(1 - object$muhat[, i])))),
          gamma = -2 * (log(ifelse(series[, i] == 0, 1, series[, i]/object$muhat[, i])) -
              (series[, i] - object$muhat[, i])/object$muhat[, i]), `negative binomial` = 2 *
            (series[, i] * log(pmax(1, series[, i])/object$muhat[, i]) -
                (series[, i] + object$model$u[, i]) * log((series[, i] + object$model$u[, i])/(object$muhat[, i] + object$model$u[, i])))))
      series
    }
  )
  ts(drop(series), start = start(object$model$y), frequency = frequency(object$model$y),
    names = colnames(object$model$y))
}
