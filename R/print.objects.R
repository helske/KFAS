#' Print SSModel Object
#' @export
#' @param x SSModel object
#' @param ... Ignored.
print.SSModel <-  function(x, ...) {

  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
    sep = "")
  cat("State space model object of class SSModel\n\n")
  cat("Dimensions:\n")
  print(paste0("Number of time points: ", attr(x, "n")), quote = FALSE)
  print(paste0("Number of time series: ", attr(x, "p")), quote = FALSE)
  print(paste0("Number of disturbances: ", attr(x, "k")), quote = FALSE)
  print(paste0("Number of states: ", attr(x, "m")), quote = FALSE)
  cat("Names of the states:\n")
  print.default(format(rownames(x$a1)), quote = FALSE, print.gap = 2L)
  cat("Distributions of the time series:\n")
  if(length(distr <- unique(x$distribution)) > 1)
    distr <- x$distribution
  print.default(format(distr), quote = FALSE, print.gap = 2L)
  if (is.SSModel(x)) {
    cat("\nObject is a valid object of class SSModel.")
  } else {
    is.SSModel(x, return.logical = FALSE)
  }
}
#' Print Ouput of Kalman Filter and Smoother
#'
#' @export
#' @param x output object from function KFS.
#' @param type What to print. Possible values are \code{"state"} (default),
#'   \code{"signal"}, and  \code{"mean"}. Multiple choices are allowed.
#' @param digits minimum number of digits to be printed.
#' @param ... Ignored.
print.KFS <-
  function(x, type = "state", digits = max(3L, getOption("digits") - 3L), ...) {
    
    p <- x$dims$p
    m <- x$dims$m
    n <- x$dims$n
    type <- match.arg(type,choices = c("state", "signal", "mean"),
      several.ok = TRUE)

    pdiag <- 1 + 0:(p - 1) * (p + 1)
    mdiag <- 1 + 0:(m - 1) * (m + 1)
    namesx <- names(x)
    if("state" %in% type && any(c("a", "alphahat") %in% namesx)) {
      if (!("alphahat" %in% namesx)) {
        gaussian<-all(x$all_gaussian)
        print_this <- cbind(x$a[n + gaussian, ], sqrt(x$P[, , n + gaussian][mdiag]))
        colnames(print_this) <- c("Estimate", "Std. Error")
        if(gaussian){
          cat(paste0("Filtered values of states and standard errors at time n+1 = ",n+1,":\n"))
        } else cat(paste0("\n Filtered values of states and standard errors at time n = ",n,":\n"))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
      } else {
        print_this <- cbind(x$alphahat[n, ], sqrt(x$V[, , n][mdiag]))
        colnames(print_this) <- c("Estimate", "Std. Error")
        cat(paste0("Smoothed values of states and standard errors at time n = ",n,":\n"))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
      }
      cat("\n")
    }
    if("signal" %in% type && any(c("t", "thetahat") %in% namesx)){
      if (!("thetahat" %in% namesx)) {
        print_this <- cbind(x$t[n, ], sqrt(x$P_theta[, , n][pdiag]))
        colnames(print_this) <- c("Estimate", "Std. Error")
        cat(paste0("Filtered values of signal and standard errors at time n = ",n,":\n"))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
      } else {
        print_this <- cbind(x$thetahat[n,], sqrt(x$V_theta[, , n][pdiag]))
        colnames(print_this) <- c("Estimate", "Std. Error")
        cat(paste0("Smoothed values of signal and standard errors at time n = ",n,":\n"))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
      }
      cat("\n")
    }
    if("mean" %in% type && any(c("m", "muhat") %in% namesx)){
      if (!("muhat" %in% namesx)) {
        print_this <- cbind(x$m[n, ], sqrt(x$P_mu[, , n][pdiag]))
        colnames(print_this) <- c("Estimate", "Std. Error")
        cat(paste0("Filtered values of mean and standard errors at time n = ",n,":\n"))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
      } else {
        print_this <- cbind(x$muhat[n,], sqrt(x$V_mu[, , n][pdiag]))
        colnames(print_this) <- c("Estimate", "Std. Error")
        cat(paste0("Smoothed values of mean and standard errors at time n = ",n,":\n"))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
      }
      cat("\n")
    }

  }
