#' @rdname SSModel
#' @param f An user-defined function which should return the output from 
#' \code{\link{SSMcustom}}. See examples.
#' @examples
#' # Example on bespoke function for time-varying trend
#' trend <- function(sigma, n) {
#'   Z <- array(seq_len(n), c(1, 1, n))
#'   T <- R <- matrix(1, 1, 1)
#'   Q <- matrix(sigma^2, 1, 1)
#'   a1 <- 0
#'   P1 <- 10
#'   SSMcustom(Z, T, R, Q, a1, P1, n = n, state_names = "timevarying trend")
#' }
#' 
#' model <- SSModel(Nile ~ SSMbespoke(trend(NA, length(Nile))), H = NA)
#' updatefn <- function(pars, model){
#'   model$Q[1, 1, 1] <- exp(0.5 * pars[1])
#'   model$H[1, 1, 1] <- exp(0.5 * pars[2])
#'   model
#' }
#' 
#' fit <- fitSSM(model, c(1, 20), updatefn, method = "BFGS")
#' conf_intv <- predict(fit$model, interval = "confidence", level = 0.95)
#' 
#' ts.plot(
#'   cbind(Nile, conf_intv), 
#'   col = c(1, 2, 2, 2),
#'   ylab = "Predicted Annual flow", 
#'   main = "River Nile"
#' ) 
#' @export
SSMbespoke <- function(f, index, n) {
  f
}
