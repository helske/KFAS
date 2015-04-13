#' Extract Estimated States of State Space Model
#' 
#' Extracts the estimated states from output of \code{KFS}. For non-Gaussian 
#' models without simulation, these are the estimates of conditional modes of 
#' states. For Gaussian models and non-Gaussian models with importance sampling,
#' these are the estimates of conditional means of states.
#' 
#' @export
#' @param object An object of class \code{KFS}.
#' @param start The start time of the period of interest. Defaults to first time
#'  point of the object.
#' @param end The end time of the period of interest. Defaults to the last time 
#' point of the object.
#' @param filtered Logical, return filtered instead of smoothed estimates of 
#' state vector. Default is \code{FALSE}.
#' @param states Which states to extract? Either a numeric vector containing
#'  the indices of the corresponding states, or a character vector defining the
#'  types of the corresponding states. Possible choices are
#'   \code{"all"},  \code{"level"}, \code{"slope"}, 
#'   \code{"trend"},  \code{"regression"}, \code{"arima"}, \code{"custom"}, 
#'   \code{"cycle"} or \code{"seasonal"}, where \code{"trend"} extracts states relating to trend.
#'    These can be combined. Default is \code{"all"}.
#' @param \dots Ignored.
#' @return Multivariate time series containing estimates states.
coef.KFS <- function(object, start = NULL, end = NULL, filtered = FALSE, states = "all", ...) {
  if (!filtered) {
    if (!is.null(object$alphahat)) {
      tmp <- object$alphahat
    } else stop("Input does not contain smoothed estimates for states, rerun KFS with state smoothing.")
  } else {
    if (!is.null(object[["a", exact = TRUE]])) {
      tmp <- object$a
    } else stop("Input does not contain filtered estimates for states, rerun KFS with state filtering.")
  }
  states <- match.arg(arg = states, choices = c("all", "arima", "custom", "level","slope",
                                                "cycle", "seasonal", "trend", "regression"),
                      several.ok = TRUE)
  
  if ("all" %in% states) {
    states <- as.integer(1:attr(object$model, "m"))
  } else {
    if("trend" %in% states)
      states <- c(states, "level", "slope")
    states <- which(attr(object$model, "state_types") %in% states)
  }
  tmp <- window(tmp[,states], start = start, end = end)
  if (!is.null(start) && start == end) 
    tsp(tmp) <- class(tmp) <- NULL
  drop(tmp)
} 
