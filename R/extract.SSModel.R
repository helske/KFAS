#' Extract or Replace Parts of a State Space Model
#'
#' S3 methods for getting and setting parts of object of class
#' \code{SSModel}. These methods ensure that dimensions of system matrices are
#' not altered.
#'
#' @export
#' @rdname Extract.SSModel
#' @param x Object of class \code{SSModel}.
#' @param element Which element is chosen. Possible choices are \code{"y"}, \code{"Z"}, \code{"H"},
#'   \code{"T"}, \code{"R"}, \code{"Q"}, \code{"a1"}, \code{"P1"}, \code{"P1inf"}, and \code{"u"}.
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
  
  element <- match.arg(arg = element,
    choices = c("y", "Z", "H", "T", "R", "Q", "a1", "P1",
      "P1inf", "u"))
  if (!(element %in% c("y", "u", "Q"))) {
    if (missing(states)) {
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
        
        if ("all" %in% states) {
          states <- 1:attr(x, "m")
        } else {
          if ("trend" %in% states)
            states <- c(states, "level", "slope")
          states <- which(attr(x, "state_types") %in% states)
        }
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
          stop("Vector etas should contain the indices or names of the etas
                 which are modified.")
      } else {
        etas <- match.arg(arg = etas,
          choices = c("all", "arima", "custom", "cycle",
            "seasonal", "trend", "level", "slope", "regression"),
          several.ok = TRUE)
        
        if ("all" %in% etas) {
          etas <- 1:attr(x, "k")
        } else {
          if ("trend" %in% etas)
            etas <- c(etas, "level", "slope")
          etas <- which(attr(x, "eta_types") %in% etas)
        }
      }
    }
  }
  if (element %in% c("y", "u", "Z")) {
    if (missing(series)) {
      series <- 1:attr(x, "p")
    } else if (!all(series %in% (1:attr(x, "p"))))
      stop("Argument series must have values between 1 to p, where p is the
             number of time series in model. ")
  }
  if (missing(times)) {
    switch(element, y = , u = x[[element]][, series] <- value,
      Z = x[[element]][series, states, ] <- value,
      H = x[[element]][series, series, ] <- value,
      T = x[[element]][states, states, ] <- value,
      R = x[[element]][states, etas, ] <- value,
      Q = x[[element]][etas, etas, ] <- value,
      a1 = x[[element]][states, 1] <- value,
      P1 = x[[element]][states, states] <- value,
      P1inf = x[[element]][states, states] <- value, )
  } else {
    switch(element, y = , u = x[[element]][times, series] <- value,
      Z = x[[element]][series, states, times] <- value,
      H = x[[element]][series, series, times] <- value,
      T = x[[element]][states, states, times] <- value,
      R = x[[element]][states, etas, times] <- value,
      Q = x[[element]][etas, etas, times] <- value,
      a1 = x[[element]][states, 1] <- value,
      P1 = x[[element]][states, states] <- value,
      P1inf = x[[element]][states, states] <- value, )
  }
  x
  
}
#' @export
#' @rdname Extract.SSModel
`[.SSModel` <-
  function(x, element, states, etas, series, times, drop=TRUE, ...) {
    
    if (is.numeric(element)) {
      unclass(x)[element]
    } else {
      
      element <- match.arg(arg = element,
        choices = c("y", "Z", "H", "T", "R", "Q", "a1", "P1",
          "P1inf", "u"))
      
      if (!(element %in% c("y", "u", "Q"))) {
        if (missing(states)) {
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
            
            if ("all" %in% states) {
              states <- 1:attr(x, "m")
            } else {
              if("trend" %in% states)
                states <- c(states, "level", "slope")
              states <- which(attr(x, "state_types") %in% states)
            }
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
              stop("Vector etas should contain the indices or types of the etas
                 which are modified.")
          } else {
            etas <- match.arg(arg = etas,
              choices = c("all", "arima", "custom", "cycle",
                "seasonal", "trend", "level", "slope", "regression"),
              several.ok = TRUE)
            
            if ("all" %in% etas) {
              etas <- 1:attr(x, "k")
            } else {
              if("trend" %in% etas)
                etas <- c(etas, "level", "slope")
              etas <- which(attr(x, "eta_types") %in% etas)
            }
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
        switch(element, y = ,
          u = x[[element]][, series, drop = drop],
          Z = x[[element]][series, states, , drop = drop],
          H = x[[element]][series, series, , drop = drop],
          T = x[[element]][states, states, , drop = drop],
          R = x[[element]][states, etas, , drop = drop],
          Q = x[[element]][etas, etas, , drop = drop],
          a1 = x[[element]][states, 1, drop = drop],
          P1 = x[[element]][states, states, drop = drop],
          P1inf = x[[element]][states, states, drop = drop],
        )
      } else {
        switch(element, y = ,
          u = x[[element]][times, series, drop = drop],
          Z = x[[element]][series, states, times, drop = drop],
          H = x[[element]][series, series, times, drop = drop],
          T = x[[element]][states, states, times, drop = drop],
          R = x[[element]][states, etas, times, drop = drop],
          Q = x[[element]][etas, etas, times, drop = drop],
          a1 = x[[element]][states, 1, drop = drop],
          P1 = x[[element]][states, states, drop = drop],
          P1inf = x[[element]][states, states, drop = drop], )
      }
    }
  }
