#'Extracting the Partial Signal Of a State Space Model
#'
#'Function \code{signal} returns the signal of a state space model using only
#'subset of states.
#'
#'@export
#'@param object Object of class \code{KFS}.
#'@param states Which states are combined? Either a numeric vector containing
#'  the indices of the corresponding states, or a character vector defining the
#'  types of the corresponding states. Possible choices are
#'   \code{"all"},  \code{"level"}, \code{"slope"}, 
#'   \code{"trend"},  \code{"regression"}, \code{"arima"}, \code{"custom"}, 
#'   \code{"cycle"} or \code{"seasonal"}, where \code{"trend"} extracts states relating to trend.
#'    These can be combined. Default is \code{"all"}.
#'@param filtered If \code{TRUE}, filtered signal is used. Otherwise smoothed signal is
#'  used.
#'@return 
#'\item{signal}{Time series object of filtered signal \eqn{Z_ta_t}{Z[t]a[t]} or 
#'smoothed signal \eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]} using only the defined states. } 
#'\item{variance}{Cov(\eqn{Z_ta_t}{Z[t]a[t]}) or Cov(\eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]}) using only the defined states. 
#'For the covariance matrices of the filtered signal, only the non-diffuse part of P is used.  }
signal <- 
  function(object, states = "all", filtered = FALSE) {
    if (!inherits(object, "KFS")) 
        stop("Object must be an output from function KFS.")
    if (is.numeric(states)) {
        states <- as.integer(states)
        if (min(states) < 1 | max(states) > attr(object$model, "m")) 
            stop("Vector states should contain the indices or names of the states which are combined.")
    } else {
        states <- match.arg(arg = states, choices = c("all", "arima", "custom", "level","slope",
                                                      "cycle", "seasonal", "trend", "regression"),
                            several.ok = TRUE)
        if("trend" %in% states)
          states <- c(states, "level", "slope")
        if ("all" %in% states) {
            states <- as.integer(1:attr(object$model, "m"))
        } else states <- which(attr(object$model, "state_types") %in% states)
    }
    if (!isTRUE(length(states) > 0)) 
        stop("Selected states not in the model.")
    if (identical(states, as.integer(1:attr(object$model, "m")))) {
        if (all(object$model$distribution == "gaussian")) {
            if (filtered && !is.null(object$m)) {
                return(list(signal = object$m, variance = object$P_mu))
            } else {
                if (!is.null(object$muhat)) 
                  return(list(signal = object$muhat, variance = object$V_muhat))
            }
        } else {
            if (filtered && !is.null(object$t)) {
                return(list(signal = object$t, variance = object$P_theta))
            } else {
                if (!is.null(object$thetahat)) 
                  return(list(signal = object$thetahat, variance = object$V_theta))
            }
        }
    }
    if (filtered) {
        if (is.null(object[["a", exact = TRUE]])) 
            stop("Object does not contain filtered estimates of states.")
        a <- object$a
        P <- object$P
    } else {
        if (is.null(object$alphahat)) 
            stop("Object does not contain smoothed estimates of states.")
        a <- object$alphahat
        P <- object$V
    }
    signal <- .Fortran(fsignaltheta, NAOK = TRUE, as.integer(dim(object$model$Z)[3] > 
        1), object$model$Z, t(a)[1:attr(object$model, "m"), 1:attr(object$model, 
        "n")], P[1:attr(object$model, "m"), 1:attr(object$model, "m"), 1:attr(object$model, 
        "n")], as.integer(attr(object$model, "p")), as.integer(attr(object$model, 
        "n")), as.integer(attr(object$model, "m")), theta = array(0, c(attr(object$model, 
        "n"), attr(object$model, "p"))), V_theta = array(0, c(attr(object$model, 
        "p"), attr(object$model, "p"), attr(object$model, "n"))), d = as.integer(0), 
        states, as.integer(length(states)))
    attributes(signal$theta) <- attributes(object$model$y)
    list(signal = signal$theta, variance = signal$V_theta)
} 
