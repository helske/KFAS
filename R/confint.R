#' Confidence Intervals of Smoothed States
#'
#' Extract confidence intervals of the smoothed estimates of states from the 
#' output of \code{KFS}.
#'
#' @export
#' @importFrom stats confint window
#' @name confint.KFS
#' @param object An object of class \code{KFS}.
#' @param parm Which states to extract? Either a numeric vector containing
#'  the indices of the corresponding states, or a character vector defining the
#'  types of the corresponding states. Possible choices are
#'   \code{"all"},  \code{"level"}, \code{"slope"},
#'   \code{"trend"},  \code{"regression"}, \code{"arima"}, \code{"custom"},
#'   \code{"cycle"} or \code{"seasonal"}, where \code{"trend"} extracts all states
#'   relating to trend. These can be combined. Default is \code{"all"}.
#' @param level The confidence level required. Defaults to 0.95.
#' @param \dots Ignored.
#' @return A list of confidence intervals for each state
#' @examples
#'
#' model <- SSModel(log(drivers) ~ SSMtrend(1, Q = list(1)) +
#'  SSMseasonal(period = 12, sea.type = "trigonometric") +
#'  log(PetrolPrice) + law, data = Seatbelts, H = 1)
#' out <- KFS(model)
#'
#' confint(out, parm = "regression", last = TRUE)
#' confint(out, start = c(1983, 12), end = c(1984, 2))
#'
confint.KFS <- function(object, parm = "all", level = 0.95, ...) {
  
  if ("alphahat" %in% names(object)) {
    V <- object$V
  } else stop("Input does not contain smoothed estimates for states, rerun KFS with state smoothing.")
  if (is.numeric(states)) {
    states <- as.integer(states)
    if (min(states) < 1 | max(states) > attr(object$model, "m"))
      stop("Vector states should contain the indices or names (state types) of the states.")
  } else {
    states <- match.arg(arg = states,
      choices = c("all", "arima", "custom", "level","slope", "cycle",
        "seasonal", "trend", "regression"),
      several.ok = TRUE)
    
    if ("all" %in% states) {
      states <- 1:object$dims$m
    } else {
      if (identical(object$model, "not stored")) 
        stop("No model stored as part of KFS, cannot infer state types.")
      if ("trend" %in% states) {
        states <- c(states, "level", "slope")
      }
      states <- which(attr(object$model, "state_types") %in% states)
    }
  }
  tmp <- sqrt(V[states, states, ])
  if (!is.null(start) && identical(start,end)) {
    tmp[1, ]
  } else {
    drop(tmp)
  }
  
}
