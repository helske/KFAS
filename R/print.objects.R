#' Print SSModel Object
#' @S3method print SSModel
#' @method print SSModel
#' @param x SSModel object
#' @param ... Ignored.

print.SSModel <- function(x, ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("State space model object of class SSModel\n\n")
    cat("Dimensions:\n")
    print(paste0("Number of time points: ", attr(x, "n")), quote = FALSE)
    print(paste0("Number of time series: ", attr(x, "p")), quote = FALSE)
    print(paste0("Number of disturbances: ", attr(x, "k")), quote = FALSE)
    print(paste0("Number of states: ", attr(x, "m")), quote = FALSE)
    cat("Names of the states:\n")
    print.default(format(rownames(x$a1)), quote = FALSE, print.gap = 2L)
    cat("Distributions of the time series:\n")
    print.default(format(x$distribution), quote = FALSE, print.gap = 2L)
    
    if (is.SSModel(x)) {
        cat("\nObject is a valid object of class SSModel.")
    } else {
        is.SSModel(x, return.logical = FALSE)
    }
}

#' Print Ouput of Kalman Filter and Smoother
#' @S3method print KFS
#' @method print KFS
#' @param x output object from function KFS.
#' @param digits minimum number of digits to be printed.
#' @param ... Ignored.

print.KFS <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    p <- attr(x$model, "p")
    m <- attr(x$model, "m")
    n <- attr(x$model, "n")
    pdiag <- 1 + 0:(p - 1) * (p + 1)
    mdiag <- 1 + 0:(m - 1) * (m + 1)
    if (!is.null(x$a) && is.null(x$alphahat)){
        
        cat("\n Filtered values of the first and last states and their standard errors:\n")
        
        print_this <- cbind(t(x$a[c(1, n + 1), , drop = FALSE]), sqrt(x$P[, , 1][mdiag]), sqrt(x$P[, , n + 1][mdiag]))
        colnames(print_this) <- c("a_1", paste0("a_", n + 1), "se_1", paste0("se_", n + 1))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
        
    }
    if (!is.null(x$theta) && is.null(x$thetahat)){
      
      cat("\n Filtered values of the first and last signals and their standard errors:\n")
      
      print_this <- cbind(t(x$theta[c(1, n), , drop = FALSE]), 
                          sqrt(x$P_theta[, , 1][pdiag]), 
                          sqrt(x$P_theta[, , n][pdiag]))
      colnames(print_this) <- c("theta_1", paste0("theta_", n), "se_1", paste0("se_", n))
      print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
      
    }
    if (!is.null(x$mu) && is.null(x$muhat)) {
      cat("\n Filtered values of the first and last mus and their standard errors:\n")
      
      print_this <- cbind(t(x$mu[c(1, n), ]), sqrt(x$P_mu[, , 1][pdiag]), sqrt(x$P_mu[, , n][pdiag]))
      colnames(print_this) <- c("mu_1", paste0("mu_", n), "se_1", paste0("se_", n))
      print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
      
    }
    
    if (!is.null(x$alphahat)) {
        cat("\n Smoothed values of the first and last states and their standard errors:\n")
        
        print_this <- cbind(t(x$alphahat[c(1, n), , drop = FALSE]), sqrt(x$V[, , 1][mdiag]), sqrt(x$V[, , n][mdiag]))
        colnames(print_this) <- c("alphahat_1", paste0("alphahat_", n), "se_1", paste0("se_", n))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
    }
    
    if (!is.null(x$thetahat)) {
        cat("\n Smoothed values of the first and last signals and their standard errors:\n")
        
        print_this <- cbind(t(x$thetahat[c(1, n), ]), sqrt(x$V_theta[, , 1][pdiag]), sqrt(x$V_theta[, , n][pdiag]))
        colnames(print_this) <- c("thetahat_1", paste0("thetahat_", n), "se_1", paste0("se_", n))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
        
    }
    if (!is.null(x$muhat)) {
        cat("\n Smoothed values of the first and last mus and their standard errors:\n")
        
        print_this <- cbind(t(x$muhat[c(1, n), ]), sqrt(x$V_mu[, , 1][pdiag]), sqrt(x$V_mu[, , n][pdiag]))
        colnames(print_this) <- c("mu_1", paste0("mu_", n), "se_1", paste0("se_", n))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
        
    }
    
    if (!is.null(x$epshat)) {
        cat("\n Smoothed values of the first and last epsilon disturbances and their standard errors:\n")
        
        print_this <- cbind(t(x$epshat[c(1, n), ]), sqrt(x$V_eps[, 1]), sqrt(x$V_eps[, n]))
        colnames(print_this) <- c("epshat_1", paste0("epshat_", n), "se_1", paste0("se_", n))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
    }
    if (!is.null(x$etahat)) {
        cat("\n Smoothed values of the first and last eta disturbances and their standard errors:\n")
        kdiag <- 1 + 0:(attr(x$model, "k") - 1) * (attr(x$model, "k") + 1)
        
        print_this <- cbind(t(x$etahat[c(1, n), ]), sqrt(x$V_eta[, , 1][kdiag]), sqrt(x$V_eta[, , n][kdiag]))
        colnames(print_this) <- c("etahat_1", paste0("etahat_", n), "se_1", paste0("se_", n))
        print.default(format(print_this, digits = digits), quote = FALSE, print.gap = 2)
    }
    
} 
