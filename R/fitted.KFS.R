#' Extract Fitted Values of State Space Model
#' 
#' Extracts fitted values from output of \code{KFS}.
#' @export
#' @inheritParams coef.KFS
#' @return Multivariate time series containing fitted values.
fitted.KFS <- 
  function(object, start = NULL, end = NULL, filtered = FALSE, ...) {
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
    if (!is.null(start) && start == end) 
        tsp(tmp) <- class(tmp) <- NULL
    drop(tmp)
} 
