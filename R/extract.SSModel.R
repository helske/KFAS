#' Extract or Replace Parts of a State Space Model
#'
#' S3 methods for extracting or replacing parts of objects of class \code{SSModel}. These methods 
#' ensure that dimensions of system matrices are not altered. \code{[} and \code{subset} and 
#' corresponding replacement methods are identical methods with different method names.
#'
#' @method [<- SSModel
#' @S3method [<- SSModel
#' @rdname Extract.SSModel
#' @param x Object of class \code{SSModel}.
#' @param element Which element is chosen. Possible choices are 'y','Z','H','T','R','Q','a1','P1','P1inf', and 'u'.
#' @param states Which states are chosen. Either a numeric vector containing the indices of the corresponding states,
#' or a character vector defining the types of the corresponding states. 
#' Possible choices are \dQuote{all}, \dQuote{arima}, \dQuote{custom}, \dQuote{cycle}, \dQuote{seasonal}, 
#' \dQuote{trend}, or \dQuote{regression}. These can be combined. Default is \dQuote{all}.
#' @param etas Which disturbances eta are chosen. Used for elements \dQuote{R} and \dQuote{Q}.Either a numeric vector containing the indices of the corresponding etas,
#' or a character vector defining the types of the corresponding etas. 
#' Possible choices are \dQuote{all}, \dQuote{arima}, \dQuote{custom}, \dQuote{cycle}, \dQuote{seasonal}, 
#' \dQuote{trend}, or \dQuote{regression}. These can be combined. 
#' @param series Numeric. Which series are chosen. Used for elements \dQuote{y}, \dQuote{Z}, and \dQuote{u}.
#' @param times Numeric. Which time points are chosen.
#' @param value A value to be assigned to x.
#' @return A selected subset of the chosen element or a value.
#' @examples
#' set.seed(1)
#' model<-SSModel(rnorm(10)~1)
#' model["H"]
#' model["H"]<-10
#' # H is still an array:
#' model["H"]
#' logLik(model)
#' model$H<-1
#' # model["H"] throws an error as H is now scalar:
#' model$H
#' logLik(model,check.model=TRUE) #with check.model=FALSE (default) R crashes!

`[<-.SSModel` <- function(x, element, states, etas, series, times, value) {
    
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
#' @method [ SSModel
#' @S3method [ SSModel
#' @rdname Extract.SSModel
`[.SSModel` <- function(x, element, states, etas, series, times) {
    
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
        switch(element, y = , u = x[[element]][, series,drop=FALSE], Z = x[[element]][series, states, ,drop=FALSE], 
               H = x[[element]][series, series, ,drop=FALSE
            ], T = x[[element]][states, states, ,drop=FALSE], R = x[[element]][states, etas, ,drop=FALSE], 
               Q = x[[element]][etas, etas,,drop=FALSE ], a1 = x[[element]][states, 
            1,drop=FALSE], P1 = x[[element]][states, states,drop=FALSE], P1inf = x[[element]][states, states,drop=FALSE], )
        
    } else {
        switch(element, y = , u = x[[element]][times, series,drop=FALSE], Z = x[[element]][series, states, times,drop=FALSE], H = x[[element]][series, 
            series, times,drop=FALSE], T = x[[element]][states, states, times,drop=FALSE], R = x[[element]][states, etas, times,drop=FALSE], Q = x[[element]][etas, 
            etas, times,drop=FALSE], a1 = x[[element]][states, 1,drop=FALSE], P1 = x[[element]][states, states,drop=FALSE], P1inf = x[[element]][states, states,drop=FALSE], 
            )
    }
    
} 
