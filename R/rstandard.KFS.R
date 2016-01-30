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
#'   \eqn{Q_t - V_{\eta,t}}{Q[t] - V[\eta,t]}, or the diagonal of the same matrix.
#'
#'   \item "pearson":  Standardized Pearson residuals
#'   \deqn{L^{-1}_t(y_{t}-\theta_{i}), \quad t=1,\ldots,n,}{L^(-1)[t](y[t]-\theta[t]), t=1,\ldots,n,} where
#'   \eqn{L_t}{L[t]} is the lower triangular matrix from Cholesky decomposition
#'   of \eqn{Var(y_t|y_{n},\ldots,y_1)}{Var(y[t]|y[n],\ldots,y[1])}, or the diagonal of the same
#'   matrix. For Gaussian models, these coincide with the standardized smoothed
#'   \eqn{\epsilon} disturbance residuals, and for generalized linear models
#'   these coincide with the standardized Pearson residuals (hence the name).
#'
#'
#' \item 'deviance': Deviance residuals. Deprecated. This option was meant to be
#' used only for the GLM comparisons, as their generalization to other models is
#' lacking, but these will be completely removed in future in order to avoid
#' misleading results in non-GLM settings. }
#'
#' @param model KFS object
#' @param type Type of residuals. See details.
#' @param standardization_type Type of standardization. Either \code{"marginal"}
#'   (default) for marginal standardization or  \code{"cholesky"} for Cholesky-type standardization.
#' @param ... Ignored.
#' @examples
#' modelNile <- SSModel(Nile ~ SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
#' modelNile <- fitSSM(inits = c(log(var(Nile)),log(var(Nile))), model = modelNile,
#'   method = "BFGS")$model
#' # Filtering and state smoothing
#' out <- KFS(modelNile, smoothing = c("state", "mean", "disturbance"))
#'
#' plot(cbind(
#'     recursive = rstandard(out),
#'     irregular = rstandard(out, "pearson"),
#'     state = rstandard(out, "state")),
#'   main = "recursive and auxiliary residuals")
rstandard.KFS <- function(model,
  type = c("recursive", "pearson", "state"),
  standardization_type = c("marginal","cholesky"), ...) {

  type <- match.arg(type)
  stype <- match.arg(standardization_type)

  if (type == "deviance")
    .Deprecated(msg="Argument type=\"deviance\" is deprecated.")

  if (type == "state" && any(model$model$distribution != "gaussian"))
    stop("State residuals are only supported for fully gaussian models.")

  x <- do.call(paste0(type,"_standardized"), list(model,stype))
  return(ts(drop(x), start = start(model$model$y), frequency = frequency(model$model$y),
    names = colnames(model$model$y)))
}
