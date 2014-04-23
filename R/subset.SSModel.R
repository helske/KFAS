#' @method subset<- SSModel
#' @S3method subset<- SSModel
#' @rdname Extract.SSModel
`subset<-.SSModel` <- function(x, element, states, etas, series, times, ..., value) {
    
    # is.SSModel(x,return.logical=FALSE)
    element <- match.arg(arg = element, choices = c("y", "Z", "H", "T", "R", "Q", "a1", "P1", "P1inf", "u"))
    
    if (!(element %in% c("y", "u", "Q"))) {
        if (missing(states)) {
            states <- 1:attr(x, "m")
        } else {
            if (is.numeric(states)) {
                states <- as.integer(states)
                if (min(states) < 1 | max(states) > attr(x, "m")) 
                  stop("Vector states should contain the indices or names of the states which are modified.")
            } else {
                states <- match.arg(arg = states, choices = c("all", "arima", "custom", "cycle", "seasonal", "trend", "regression"), 
                  several.ok = TRUE)
                if ("all" %in% states) {
                  states <- 1:attr(x, "m")
                } else states <- which(attr(x, "state_types") %in% states)
            }
        }
    }
    if (element %in% c("R", "Q")) {
        if (missing(etas)) {
            etas <- 1:attr(x, "k")
        } else {
            if (is.numeric(etas)) {
                etas <- as.integer(etas)
                if (min(etas) < 1 | max(etas) > attr(x, "k")) 
                  stop("Vector etas should contain the indices or names of the etas which are modified.")
            } else {
                etas <- match.arg(arg = etas, choices = c("all", "arima", "custom", "cycle", "seasonal", "trend", "regression"), 
                  several.ok = TRUE)
                if ("all" %in% etas) {
                  etas <- 1:attr(x, "k")
                } else etas <- which(attr(x, "eta_types") %in% etas)
            }
        }
    }
    if (element %in% c("y", "u", "Z")) {
        if (missing(series)) {
            series <- 1:attr(x, "p")
        } else if (!all(series %in% (1:attr(x, "p")))) 
            stop("Argument series must have values between 1 to p, where p is the number of time series in model. ")
    }
    
    if (missing(times)) {
        switch(element, y = , u = x[[element]][, series] <- value, Z = x[[element]][series, states, ] <- value, H = x[[element]][series, 
            series, ] <- value, T = x[[element]][states, states, ] <- value, R = x[[element]][states, etas, ] <- value, Q = x[[element]][etas, 
            etas, ] <- value, a1 = x[[element]][states, 1] <- value, P1 = x[[element]][states, states] <- value, P1inf = x[[element]][states, 
            states] <- value, )
        
    } else {
        switch(element, y = , u = x[[element]][times, series] <- value, Z = x[[element]][series, states, times] <- value, H = x[[element]][series, 
            series, times] <- value, T = x[[element]][states, states, times] <- value, R = x[[element]][states, etas, times] <- value, 
            Q = x[[element]][etas, etas, times] <- value, a1 = x[[element]][states, 1] <- value, P1 = x[[element]][states, states] <- value, 
            P1inf = x[[element]][states, states] <- value, )
    }
    
    x
}
#' @rdname Extract.SSModel
#' @export
`subset<-` <- function(x, ..., value) UseMethod("subset<-")

#' @method subset SSModel
#' @S3method subset SSModel
#' @rdname Extract.SSModel
#' @param ... ignored.
subset.SSModel <- function(x, element, states, etas, series, times, ...) {
    
    # is.SSModel(x,return.logical=FALSE)
    element <- match.arg(arg = element, choices = c("y", "Z", "H", "T", "R", "Q", "a1", "P1", "P1inf", "u"))
    
    if (!(element %in% c("y", "u", "Q"))) {
        if (missing(states)) {
            states <- 1:attr(x, "m")
        } else {
            if (is.numeric(states)) {
                states <- as.integer(states)
                if (min(states) < 1 | max(states) > attr(x, "m")) 
                  stop("Vector states should contain the indices or types of the states which are modified.")
            } else {
                states <- match.arg(arg = states, choices = c("all", "arima", "custom", "cycle", "seasonal", "trend", "regression"), 
                  several.ok = TRUE)
                if ("all" %in% states) {
                  states <- 1:attr(x, "m")
                } else states <- which(attr(x, "state_types") %in% states)
            }
        }
    }
    if (element %in% c("R", "Q")) {
        if (missing(etas)) {
            etas <- 1:attr(x, "k")
        } else {
            if (is.numeric(etas)) {
                etas <- as.integer(etas)
                if (min(etas) < 1 | max(etas) > attr(x, "k")) 
                  stop("Vector etas should contain the indices or types of the etas which are modified.")
            } else {
                etas <- match.arg(arg = etas, choices = c("all", "arima", "custom", "cycle", "seasonal", "trend", "regression"), 
                  several.ok = TRUE)
                if ("all" %in% etas) {
                  etas <- 1:attr(x, "k")
                } else etas <- which(attr(x, "eta_types") %in% etas)
            }
        }
    }
    if (element %in% c("y", "u", "Z")) {
        if (missing(series)) {
            series <- 1:attr(x, "p")
        } else if (!all(series %in% (1:attr(x, "p")))) 
            stop("Argument series must have values between 1 to p, where p is the number of time series in model. ")
    }
    
    if (missing(times)) {
        switch(element, y = , u = x[[element]][, series,drop=FALSE], Z = x[[element]][series, states, ,drop=FALSE], H = x[[element]][series, series,,drop=FALSE 
            ], T = x[[element]][states, states, ,drop=FALSE], R = x[[element]][states, etas,,drop=FALSE ], Q = x[[element]][etas, etas, ,drop=FALSE], a1 = x[[element]][states, 
            1,drop=FALSE], P1 = x[[element]][states, states,drop=FALSE], P1inf = x[[element]][states, states,drop=FALSE], )
        
    } else {
        switch(element, y = , u = x[[element]][times, series,drop=FALSE], Z = x[[element]][series, states, times,drop=FALSE], H = x[[element]][series, 
            series, times,drop=FALSE], T = x[[element]][states, states, times,drop=FALSE], R = x[[element]][states, etas, times,drop=FALSE], Q = x[[element]][etas, 
            etas, times,drop=FALSE], a1 = x[[element]][states, 1,drop=FALSE], P1 = x[[element]][states, states,drop=FALSE], P1inf = x[[element]][states, states,drop=FALSE], 
            )
    }
    
} 
