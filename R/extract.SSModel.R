#' Extract or Replace Parts of a State Space Model
#'
#' S3 methods for getting and setting parts of object of class
#' \code{SSModel}. These methods ensure that dimensions of system matrices are
#' not altered.
#'
#' If \code{element} is not one of 
#' \code{"y"}, \code{"Z"}, \code{"H"}, \code{"T"}, \code{"R"}, \code{"Q"}, 
#' \code{"a1"}, \code{"P1"}, \code{"P1inf"}, \code{"u"},
#' the default single bracket list extraction 
#' and assignments (\code{x[element]} and \code{x[element] <- value}) 
#' are used (and other arguments are ignored). 
#' 
#' If \code{element} is one of 
#' \code{"y"}, \code{"Z"}, \code{"H"}, \code{"T"}, \code{"R"}, \code{"Q"}, 
#' \code{"a1"}, \code{"P1"}, \code{"P1inf"}, \code{"u"} and if the arguments 
#' \code{states}, \code{etas}, \code{times} and \code{series} are 
#' all missing, the double bracket list 
#' extraction \code{x[[element]]} and modified double bracket list assignment
#' \code{x[[element]][] <- value} are used.
#' 
#' If neither of above holds, then for example in case of \code{element = Z} 
#' the extraction is of form \code{x$Z[series, states, times, drop]}.
#'   
#' @export
#' @rdname Extract.SSModel
#' @param x Object of class \code{SSModel}.
#' @param element Which element(s) is chosen.  Typical values are \code{"y"}, 
#' \code{"Z"}, \code{"H"}, \code{"T"}, \code{"R"}, \code{"Q"}, \code{"a1"}, 
#' \code{"P1"}, \code{"P1inf"}, and \code{"u"}. See details.
#' @param states Which states are chosen. Either a numeric vector containing the indices of the
#'   states, or a character vector defining the types of the states. Possible choices are
#'   \code{"all"},  \code{"level"}, \code{"slope"},
#'   \code{"trend"},  \code{"regression"}, \code{"arima"}, \code{"custom"},
#'   \code{"cycle"} or \code{"seasonal"}, where \code{"trend"} extracts all states relating to trend.
#'    These can be combined. Default is \code{"all"}.
#' @param etas Which disturbances eta are chosen. Used for elements \code{"R"} and \code{"Q"}.
#'   Either a numeric vector containing the indices of the etas, or a character vector defining the
#'   types of the etas. Possible choices are \code{"all"},  \code{"level"}, \code{"slope"},
#'   \code{"trend"},  \code{"regression"}, \code{"arima"}, \code{"custom"},
#'   \code{"cycle"} or \code{"seasonal"}, where \code{"trend"} extracts all etas relating to trend.
#' These can be combined. Default is \code{"all"}.
#' @param series Numeric. Which series are chosen. Used for elements
#' \code{"y"}, \code{"Z"}, and \code{"u"}.
#' @param times Numeric. Which time points are chosen.
#' @param drop Logical. If \code{TRUE} (default) the result is coerced to the lowest possible
#' dimension.
#' @param value A value to be assigned to x.
#' @param ... Ignored.
#' @return A selected subset of the chosen element or a value.
#' @examples
#' set.seed(1)
#' model <- SSModel(rnorm(10) ~ 1)
#' model["H"]
#' model["H"] <- 10
#' # H is still an array:
#' model["H"]
#' logLik(model)
#' model$H <- 1
#' # model["H"] throws an error as H is now scalar:
#' model$H
#' logLik(model, check.model = TRUE) #with check.model = FALSE R crashes!
`[<-.SSModel` <-  function(x, element, states, etas, series, times, ...,
  value) {
  
  if (missing(element)) {
    element <- 1:length(element)
  }
    
  if (length(element) != 1 || !is.character(element)) {
    class(x) <- "list"
    x[element] <- value
    class(x) <- "SSModel"
    return(x)
  }
  
  tmp <- try(element <- match.arg(arg = element,
    choices = c("y", "Z", "H", "T", "R", "Q", "a1", "P1",
      "P1inf", "u")), silent = TRUE)
  
  if (inherits(tmp, "try-error")) {
    class(x) <- "list"
    x[element] <- value
    class(x) <- "SSModel"
    return(x)
  }
  if (missing(states) && missing(etas) && missing(series) && missing(times)) {
    x[[element]][] <- value
    return(x)
  }
  if (!(element %in% c("y", "u", "Q", "H"))) {
    if (missing(states) || "all" %in% states) {
      states <- 1:attr(x, "m")
    } else {
      if (is.numeric(states)) {
        states <- as.integer(states)
        if (min(states) < 1 | max(states) > attr(x, "m"))
          stop("Vector states should contain the indices or names of the
                 states which are modified.")
      } else {
        states <- match.arg(arg = states,
          choices = c("all", "arima", "custom", "cycle",
            "seasonal", "trend", "level", "slope", "regression"),
          several.ok = TRUE)
        
        
        if ("trend" %in% states) {
          states <- c(states, "level", "slope")
        }
        states <- which(attr(x, "state_types") %in% states)
        
      }
    }
  }
  if (element %in% c("R", "Q")) {
    if (missing(etas) || "all" %in% etas) {
      etas <- 1:attr(x, "k")
    } else {
      if (is.numeric(etas)) {
        etas <- as.integer(etas)
        if (min(etas) < 1 | max(etas) > attr(x, "k"))
          stop("Vector etas should contain the indices or names of the etas
                 which are modified.")
      } else {
        etas <- match.arg(arg = etas,
          choices = c("all", "arima", "custom", "cycle",
            "seasonal", "trend", "level", "slope", "regression"),
          several.ok = TRUE)
        
        
        if ("trend" %in% etas) {
          etas <- c(etas, "level", "slope")
        }
        etas <- which(attr(x, "eta_types") %in% etas)
        
      }
    }
  }
  if (element %in% c("y", "u", "Z", "H")) {
    if (missing(series)) {
      series <- 1:attr(x, "p")
    } else if (!all(series %in% (1:attr(x, "p"))))
      stop("Argument series must have values between 1 to p, where p is the
             number of time series in model. ")
  }
  if (missing(times)) {
    switch(element, 
      y = x$y[, series] <- value, 
      u = x$u[, series] <- value,
      Z = x$Z[series, states, ] <- value,
      H = x$H[series, series, ] <- value,
      T = x$T[states, states, ] <- value,
      R = x$R[states, etas, ] <- value,
      Q = x$Q[etas, etas, ] <- value,
      a1 = x$a1[states, 1] <- value,
      P1 = x$P1[states, states] <- value,
      P1inf = x$P1inf[states, states] <- value)
  } else {
    switch(element, 
      y = x$y[times, series] <- value,
      u = x$u[times, series] <- value,
      Z = x$Z[series, states, times] <- value,
      H = x$H[series, series, times] <- value,
      T = x$T[states, states, times] <- value,
      R = x$R[states, etas, times] <- value,
      Q = x$Q[etas, etas, times] <- value,
      a1 = x$a1[states, 1] <- value,
      P1 = x$P1[states, states] <- value,
      P1inf = x$P1inf[states, states] <- value)
  }
  
  x
  
}
#' @export
#' @rdname Extract.SSModel
`[.SSModel` <- function(x, element, states, etas, series, times, drop = TRUE, ...) {
  
  if (missing(element)) {
    element <- 1:length(element)
  }
  
  if (length(element) != 1 || !is.character(element)) {
    class(x) <- "list"
    return(x[element])
  }
  
  tmp <- try(element <- match.arg(arg = element,
    choices = c("y", "Z", "H", "T", "R", "Q", "a1", "P1",
      "P1inf", "u")), silent = TRUE)
  
  if (inherits(tmp, "try-error")) {
    class(x) <- "list"
    return(x[element])
  } 
  
  if (missing(states) && missing(etas) && missing(series) && missing(times)) {
    return(x[[element]])
  }
  if (!(element %in% c("y", "u", "Q", "H"))) {
    if (missing(states) || "all" %in% states) {
      states <- 1:attr(x, "m")
    } else {
      if (is.numeric(states)) {
        states <- as.integer(states)
        if (min(states) < 1 | max(states) > attr(x, "m"))
          stop("Vector states should contain the indices or types of the
                 states which are modified.")
      } else {
        states <- match.arg(arg = states,
          choices = c("all", "arima", "custom", "cycle",
            "seasonal", "trend", "level", "slope","regression"),
          several.ok = TRUE)
        if ("trend" %in% states){
          states <- c(states, "level", "slope")
        }
        states <- which(attr(x, "state_types") %in% states)
      }
    }
  }
  if (element %in% c("R", "Q")) {
    if (missing(etas) || "all" %in% etas) {
      etas <- 1:attr(x, "k")
    } else {
      if (is.numeric(etas)) {
        etas <- as.integer(etas)
        if (min(etas) < 1 | max(etas) > attr(x, "k"))
          stop("Vector etas should contain the indices or types of the etas
                 which are modified.")
      } else {
        etas <- match.arg(arg = etas,
          choices = c("all", "arima", "custom", "cycle",
            "seasonal", "trend", "level", "slope", "regression"),
          several.ok = TRUE)
        
        
        if ("trend" %in% etas)
          etas <- c(etas, "level", "slope")
        etas <- which(attr(x, "eta_types") %in% etas)
        
      }
    }
  }
  if (element %in% c("y", "u", "Z", "H")) {
    if (missing(series)) {
      series <- 1:attr(x, "p")
    } else if (!all(series %in% (1:attr(x, "p"))))
      stop("Argument series must have values between 1 to p, where p is the number of time series in model. ")
    if ((element == "u" && identical(x$u, "Omitted")) || 
        (element == "H" && identical(x$H, "Omitted"))) {
      return("Omitted")
    }
  }
  if (missing(times)) {
    switch(element, 
      y = x$y[, series, drop = drop],
      u = x$u[, series, drop = drop],
      Z = x$Z[series, states, , drop = drop],
      H = x$H[series, series, , drop = drop],
      T = x$T[states, states, , drop = drop],
      R = x$R[states, etas, , drop = drop],
      Q = x$Q[etas, etas, , drop = drop],
      a1 = x$a1[states, 1, drop = drop],
      P1 = x$P1[states, states, drop = drop],
      P1inf = x$P1inf[states, states, drop = drop])
  } else {
    switch(element, 
      y = x$y[times, series, drop = drop],
      u = x$u[times, series, drop = drop],
      Z = x$Z[series, states, times, drop = drop],
      H = x$H[series, series, times, drop = drop],
      T = x$T[states, states, times, drop = drop],
      R = x$R[states, etas, times, drop = drop],
      Q = x$Q[etas, etas, times, drop = drop],
      a1 = x$a1[states, 1, drop = drop],
      P1 = x$P1[states, states, drop = drop],
      P1inf = x$P1inf[states, states, drop = drop])
  }
}



