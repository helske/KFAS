#' Extract Estimated States of State Space Model
#' 
#' Extracts the estimates states from output of \code{KFS}. For non-Gaussian 
#' models without simulation, these are estimates of conditional modes of 
#' states. For Gaussian models and non-Gaussian models with importance sampling,
#' these are estimates of conditional means of states.
#' 
#' @export
#' @param object An object of class \code{KFS}.
#' @param start The start time of the period of interest. Defaults to first time
#'  point of the object.
#' @param end The end time of the period of interest. Defaults to the last time 
#' point of the object.
#' @param filtered Logical, return filtered instead of smoothed estimates of 
#' state vector. Default is \code{FALSE}.
#' @param \dots Ignored.
#' @return Multivariate time series containing estimates states.
coef.KFS <-
  function(object, start = NULL, end = NULL, filtered = FALSE, ...) {
    if (!filtered) {
        if (!is.null(object$alphahat)) {
            tmp <- object$alphahat
        } else stop("Input does not contain smoothed estimates for states, rerun
                    KFS with state smoothing.")
    } else {
        if (!is.null(object[["a", exact = TRUE]])) {
            tmp <- object$a
        } else stop("Input does not contain filtered estimates for states, rerun
                    KFS with state filtering.")
    }
    tmp <- window(tmp, start = start, end = end)
    if (start == end && !is.null(start)) 
        tsp(tmp) <- class(tmp) <- NULL
    drop(tmp)
} 
