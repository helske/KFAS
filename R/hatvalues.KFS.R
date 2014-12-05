#' Extract Hat Values from KFS Output
#'
#' Extract hat values from KFS output, when \code{KFS} was run with signal (non-Gaussian case) 
#' or mean smoothing (Gaussian case). 
#' 
#' @details Hat values are defined as the diagonal elements of 
#' \code{V_t/H_t} where V_t is the covariance matrix of signal/mean at time t and H_t is the 
#' covariance matrix of disturbance vector \eqn{\epsilon} of (approximating) Gaussian model 
#' at time t.
#' @export
#' @import stats
#' @param model An object of class \code{KFS}.
#' @param \dots Additional arguments to \code{approxSSM}.
#' @return Multivariate time series containing hat values.
hatvalues.KFS <- 
  function(model, ...) {
    if (any(model$model$distribution != "gaussian")) {
        app <- approxSSM(model$model, ...)
        if (is.null(model$V_theta)) 
            stop("KFS was run without signal smoothing, cannot compute hat values.")
        hatv <- matrix(apply(model$V_theta/app$H, 3, diag), attr(model$model, "n"), 
            attr(model$model, "p"), byrow = TRUE)
    } else {
        if (is.null(model$V_mu)) 
            stop("KFS was run without mean smoothing, cannot compute hat values.")
        hatv <- matrix(apply(model$V_mu, 3, diag), attr(model$model, "n"), attr(model$model, 
            "p"), byrow = TRUE)/matrix(apply(model$model$H, 3, diag), attr(model$model, 
            "n"), attr(model$model, "p"), byrow = TRUE)
    }
    attributes(hatv) <- attributes(model$model$y)
    hatv
} 
