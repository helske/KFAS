#' Extract Standardized Residuals from KFS output
#' @export
#' @importFrom stats rstandard
#' @details For object of class KFS with fully Gaussian observations, several
#'   types of standardized residuals can be computed. Also the standardization
#'   for multivariate residuals can be done either by Cholesky decomposition
#'   \eqn{L^{-1}_t residual_t,}{L^(-1)[t]residual[t]} or component-wise
#'   \eqn{residual_t/sd(residual_t),}{residual[t]/sd(residual[t])}.
#'
#'   \itemize{
#'   \item "recursive": For Gaussian models the vector standardized one-step-ahead prediction
#'   residuals are defined as
#'   \deqn{v_{t,i}/\sqrt{F_{i,t}},}{v[t,i]/\sqrt{F[i,t]},} with residuals
#'   being undefined in diffuse phase. Note that even in multivariate case these
#'   standardized residuals coincide with the ones obtained from the Kalman
#'   filter without the sequential processing (which is not true for the
#'   non-standardized innovations).
#'   For non-Gaussian models the vector standardized recursive residuals are obtained as
#'   \deqn{L^{-1}_t (y_{t}-\mu_{t}),}{L^(-1)[t](y[t]-\mu[t]),} where
#'   \eqn{L_t}{L[t]} is the lower triangular matrix from Cholesky decomposition
#'   of \eqn{Var(y_t|y_{t-1},\ldots,y_1)}{Var(y[t]|y[t-1],\ldots,y[1])}. Computing these for large
#'   non-Gaussian models can be time consuming as filtering is needed.
#'
#'   For Gaussian models the component-wise standardized one-step-ahead prediction
#'   residuals are defined as
#'   \deqn{v_{t}/\sqrt{diag(F_{t})},}{v[t])/\sqrt{diag(F[t])},} where \eqn{v_{t}}{v[t]} and
#'   \eqn{F_{t}}{F[t]} are based on the standard multivariate processing.
#'   For non-Gaussian models these are obtained as
#'   \deqn{(y_{t}-\mu_{t})/\sqrt{diag(F_t)},}{(y[t]-\mu[t])/\sqrt{diag(F[t])},} where
#'   \eqn{F_t=Var(y_t|y_{t-1},\ldots,y_1)}{F[t]=Var(y[t]|y[t-1],\ldots,y[1])}.
#'
#'   \item "state":  Residuals based on the smoothed state disturbance terms
#'   \eqn{\eta} are defined as \deqn{L^{-1}_t \hat \eta_t, \quad
#'   t=1,\ldots,n,}{L^{-1}[t] \hat \eta[t], t=1,\ldots,n,} where \eqn{L_t}{L[t]} is
#'   either the lower triangular matrix from Cholesky decomposition of
#'   \eqn{Var(\hat\eta_t) = Q_t - V_{\eta,t}}{Var(\hat\eta[t] = Q[t] - V[\eta,t]}, or the diagonal of the same
#'   matrix.
#'
#'   \item "pearson":  Standardized Pearson residuals
#'   \deqn{L^{-1}_t(y_{t}-\theta_{i}), \quad t=1,\ldots,n,}{L^(-1)[t](y[t]-\theta[t]), t=1,\ldots,n,} where
#'   \eqn{L_t}{L[t]} is the lower triangular matrix from Cholesky decomposition
#'   of \eqn{Var(\hat\mu_t) = H_t - V_{\mu,t}}{Var(\hat\mu[t]) = H[t] - V[\mu,t]}, or the diagonal of the same
#'   matrix. For Gaussian models, these coincide with the standardized smoothed
#'   \eqn{\epsilon} disturbance residuals
#'   (as \eqn{V_{\mu,t} = V_{\eps,t}}{V[\mu,t]=V[\eps,t]}),
#'   and for generalized linear models
#'   these coincide with the standardized Pearson residuals (hence the name).
#'  }
#'
#' Note that the variance estimates from \code{KFS} are of form Var(x | y),
#' e.g., \code{V_eps} from \code{KFS} is \eqn{Var(\eps_t | Y)}{Var(\eps_t | Y)}
#' and matches with equation 4.69 in Section 4.5.3 of Durbin and Koopman (2012).
#' (in case of univariate Gaussian model). But for the standardization we need
#' Var(E(x | y)) (e.g., Var(\code{epshat}) which we get with the law
#' of total variance as \eqn{H_t - V_eps}{H[t] - V[eps]} for example.
#'
#' @param model KFS object
#' @param type Type of residuals. See details.
#' @param standardization_type Type of standardization. Either \code{"marginal"}
#'   (default) for marginal standardization or  \code{"cholesky"} for Cholesky-type standardization.
#' @param zerotol Tolerance parameter for positivity checking in standardization. Default is zero.
#' The values of D <= zerotol * max(D, 0) are deemed to zero.
#' @param expected Logical value defining the approximation of H_t in case of Gamma
#' and negative binomial distribution. Default is \code{FALSE} which matches the
#' algorithm of Durbin & Koopman (1997), whereas \code{TRUE} uses the expected value
#' of observations in the equations, leading to results which match with \code{glm} (where applicable).
#' The latter case was the default behaviour of KFAS before version 1.3.8.
#' Essentially this is the difference between observed and expected information.
#' @param ... Ignored.
#' @examples
#' # Replication of residual plot of Section 8.2 of Durbin and Koopman (2012)
#' model <- SSModel(log(drivers) ~ SSMtrend(1, Q = list(NA))+
#'     SSMseasonal(period = 12, sea.type = "trigonometric", Q = NA),
#'   data = Seatbelts, H = NA)
#'
#' updatefn <- function(pars, model){
#'   model$H[] <- exp(pars[1])
#'   diag(model$Q[, , 1]) <- exp(c(pars[2], rep(pars[3], 11)))
#'   model
#' }
#'
#' fit <- fitSSM(model = model,
#'   inits = log(c(var(log(Seatbelts[, "drivers"])), 0.001, 0.0001)),
#'   updatefn = updatefn)
#'
#' # tiny differences due to different optimization algorithms
#' setNames(c(diag(fit$model$Q[,,1])[1:2], fit$model$H[1]),
#'   c("level", "seasonal", "irregular"))
#'
#' out <- KFS(fit$model, smoothing = c("state", "mean", "disturbance"))
#'
#' plot(cbind(
#'   recursive = rstandard(out),
#'   irregular = rstandard(out, "pearson"),
#'   state = rstandard(out, "state")[,1]),
#'   main = "recursive and state residuals", type = "h")
#'
#' # Figure 2.8 in DK2012
#' model_Nile <- SSModel(Nile ~
#'     SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
#' model_Nile <- fitSSM(model_Nile, c(log(var(Nile)), log(var(Nile))),
#'   method = "BFGS")$model
#'
#' out_Nile <- KFS(model_Nile,  smoothing = c("state", "mean", "disturbance"))
#'
#' par(mfrow = c(2, 2))
#' res_p <- rstandard(out_Nile, "pearson")
#' ts.plot(res_p)
#' abline(a = 0, b= 0, lty = 2)
#' hist(res_p, freq = FALSE)
#' lines(density(res_p))
#' res_s <- rstandard(out_Nile, "state")
#' ts.plot(res_s)
#' abline(a = 0, b= 0, lty = 2)
#' hist(res_s, freq = FALSE)
#' lines(density(res_s))
#'
rstandard.KFS <- function(model,
  type = c("recursive", "pearson", "state"),
  standardization_type = c("marginal","cholesky"), zerotol = 0,
  expected = FALSE,...) {

  if (identical(model$model, "not stored")) stop("No model stored as part of KFS, cannot compute residuals.")
  type <- match.arg(type)
  stype <- match.arg(standardization_type)

  if (type == "state" && any(model$model$distribution != "gaussian"))
    stop("State residuals are only supported for fully gaussian models.")

  res_names <- if (type == "state") attr(model$model, "eta_types") else colnames(model$model$y)
  x <- do.call(paste0(type,"_standardized"), list(model, stype, zerotol, expected))
  return(ts(drop(x), start = start(model$model$y), frequency = frequency(model$model$y),
    names = res_names))
}
