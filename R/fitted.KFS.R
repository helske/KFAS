#' Extract Fitted Values of State Space Model
#' 
#' Extracts fitted values from output of \code{KFS}.
#' @export
#' @param object An object of class \code{KFS}.
#' @param start The start time of the period of interest. Defaults to first time point of the object.
#' @param end The end time of the period of interest. Defaults to the last time point of the object.
#' @param filtered Logical, return filtered instead of smoothed estimates of mean vector.
#' Default is \code{FALSE}.
#' @param \dots Ignored.
#' @return Multivariate time series containing fitted values.
fitted.KFS <- 
  function(object, start = NULL, end = NULL, filtered = FALSE, ...) {
    if (!filtered) {
        if (!is.null(object$muhat)) {
            tmp <- object$muhat
        } else stop("Input does not contain smoothed estimates for means, rerun KFS with mean smoothing.")
    } else {
        if (!is.null(object$m)) {
            tmp <- object$m
        } else stop("Input does not contain filtered estimates for means, rerun KFS with mean filtering.")
    }
    tmp <- window(tmp, start = start, end = end)
    if (start == end && !is.null(start)) 
        tsp(tmp) <- class(tmp) <- NULL
    drop(tmp)
} 
