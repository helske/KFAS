#' Log-likelihood of the State Space Model.
#'
#' Function \code{logLik.SSmodel} computes the log-likelihood value of a state space model.
#'
#'
#' @export
#' @S3method logLik SSModel
#' @method logLik SSModel
#' @aliases logLik logLik.SSModel
#' @param object State space model of class \code{SSModel}.
#' @param nsim Number of independent samples used in estimating the
#' log-likelihood of the non-Gaussian state space model. Default is 0, which
#' gives good starting value for optimization. Only used for non-Gaussian model.
#' @param antithetics Logical. If TRUE, two antithetic variables are used in
#' simulations, one for location and another for scale. Default is TRUE. Only used for non-Gaussian model.
#' @param theta Initial values for conditional mode theta. Only used for non-Gaussian model.
#' @param check.model Logical. If TRUE, function \code{is.SSModel} is called before computing the likelihood. Default is FALSE.
#' @param transform How to transform the model in case of non-diagonal
#' covariance matrix \eqn{H}. Defaults to \code{'ldl'}. See function \code{\link{transformSSM}} for
#' details. 
#' @param maxiter The maximum number of iterations used in linearisation. Default is 50. Only used for non-Gaussian model.
#' @param seed The value is used as a seed via set.seed function. Only used for non-Gaussian model.
#' @param convtol Tolerance parameter for convergence checks for Gaussian approximation.
#'  Iterations are continued until 
#'  \eqn{tol>abs(dev_{old}-dev_{new})/(abs(dev_{new})+0.1))}.
#' @param ... Ignored.
#' @return \item{}{log-likelihood of the state space model.}
logLik.SSModel <- function(object, nsim = 0, antithetics = TRUE, theta, check.model = FALSE, 
                           transform = c("ldl","augment"), maxiter = 50, 
                           seed, convtol=1e-15,...) {
  if (check.model) {
    if (!is.SSModel(object, na.check = TRUE)) {
      return(-.Machine$double.xmax)
    }
  }
  
  p <- attr(object, "p")
  m <- attr(object, "m")
  k <- attr(object, "k")
  n <- attr(object, "n")
  ymiss <- array(is.na(object$y), dim = c(n, p))
  storage.mode(ymiss) <- "integer"
  tv <- array(0, dim = 5)
  tv[1] <- dim(object$Z)[3] > 1
  tv[2] <- any(object$distribution != "gaussian") || dim(object$H)[3] > 1
  tv[3] <- dim(object$T)[3] > 1
  tv[4] <- dim(object$R)[3] > 1
  tv[5] <- dim(object$Q)[3] > 1
  if (all(object$distribution == "gaussian")) {
      if(all(c(object$Q,object$H)==0) || all(c(object$R,object$H)==0)|| any(!is.finite(c(object$R,object$Q,object$H)==0)))
        return(-.Machine$double.xmax^0.75)
    kfout <- NULL
    if (p == 1) {
      kfout <- .Fortran(fglogliku, NAOK = TRUE, object$y, ymiss, as.integer(tv), object$Z, object$H, object$T, object$R, 
                        object$Q, object$a1, object$P1, object$P1inf, as.integer(m), as.integer(k),as.integer(n), lik = double(1), object$tol, as.integer(sum(object$P1inf)))
      
    } else {
      if (any(abs(apply(object$H, 3, "[", !diag(p))) > object$tol)) {
        object <- tryCatch(transformSSM(object, type = match.arg(arg = transform, choices = c("ldl","augment"))), error = function(e) e)
        if (!inherits(object, "SSModel")) {
          warning(object$message)
          return(-.Machine$double.xmax^0.75)
        }
        tv[1] <- dim(object$Z)[3] > 1
        tv[2] <- dim(object$H)[3] > 1
        tv[5] <- dim(object$Q)[3] > 1
      }
      
      kfout <- .Fortran(fgloglik, NAOK = TRUE, object$y, ymiss, as.integer(tv), object$Z, object$H, object$T, object$R, object$Q, 
                        object$a1, object$P1, object$P1inf, as.integer(p), as.integer(m), as.integer(k),as.integer(n), lik = double(1), object$tol, as.integer(sum(object$P1inf)))
    }
    logLik <- kfout$lik
  } else {
    if(all(c(object$Q,object$u)==0) || all(c(object$R,object$u)==0) || any(!is.finite(c(object$R,object$Q,object$u)==0)))
      return(-.Machine$double.xmax^0.75)
    if (missing(theta)) {
      theta <- sapply(1:p, function(i)
        switch(object$distribution[i], 
               gaussian = object$y[, i], 
               poisson = log(pmax(object$y[,i]/object$u[, i], 0.1, na.rm = TRUE)), 
               binomial = qlogis((ifelse(is.na(object$y[, i]), 0.5, object$y[, i]) + 0.5)/(object$u[, i] + 1)), 
               gamma = log(pmax(object$y[, i], 1, na.rm = TRUE)), 
               `negative binomial` = log(pmax(object$y[, i], 1/6, na.rm = TRUE))))
    } else theta <- array(theta, dim = c(n, p))
    
    if (nsim == 0) {
      nsim <- 1
      sim <- 0
      epsplus <- array(0, c(1, 1, 1))
      etaplus <- array(0, c(1, 1, 1))
      aplus1 <- array(0, dim = c(1, 1))
      c2 <- numeric(1)
      nnd <- 0
      nd <- which(diag(object$P1inf) == 0)
    } else {
      sim <- 1
      epsplus <- array(0, c(p, n, nsim))
      etaplus <- array(0, c(k, n, nsim))
      aplus1 <- array(0, dim = c(m, nsim))
      c2 <- numeric(nsim)
      
      x <- array(t(!ymiss), c(p, n, nsim))
      dfeps <- sum(x)/nsim
      
      x2 <- array(apply(object$Q, 3, diag) > object$tol, c(k, (n - 1) * tv[5] + 1))
      x2 <- array(x2, c(k, n, nsim))
      dfeta <- sum(x2)/nsim
      
      nde <- which(diag(object$P1) > object$tol)
      nnd <- length(nde)
      nd <- which(diag(object$P1inf) == 0)
      dfu <- dfeps + dfeta + nnd
      if (missing(seed)) 
        seed <- 123
      set.seed(seed)
      u <- rnorm(n = dfu * nsim, mean = 0, sd = 1)
      
      if (dfeps > 0) 
        epsplus[x] <- u[1:(dfeps * nsim)]
      if (dfeta > 0) 
        etaplus[x2] <- u[(dfeps * nsim + 1):(dfeps * nsim + dfeta * nsim)]
      if (nnd > 0) 
        aplus1[nde, ] <- u[(dfeps * nsim + dfeta * nsim + 1):(dfu * nsim)]
      
      
      if (antithetics) {
        for (i in 1:nsim) {
          u <- c(etaplus[, , i], epsplus[, , i], aplus1[, i])
          c2[i] <- t(u) %*% c(u)
        }
        q <- pchisq(c2, df = dfu)
        c2 <- sqrt(qchisq(1 - q, dfu)/c2)
      }
    }      
    nsim2 <- as.integer(max(sim * (3 * antithetics * nsim + nsim), 1))
    
    out <- .Fortran(fngloglik, NAOK = TRUE, object$y, ymiss, as.integer(tv), object$Z, object$T, object$R, object$Q, object$a1, 
                    object$P1, object$P1inf, as.integer(p), as.integer(m), as.integer(k),as.integer(n), lik = double(1), theta = theta, object$u, pmatch(x = object$distribution, table = c("gaussian", 
                                                                                                                                                                                            "poisson", "binomial", "gamma", "negative binomial"), duplicates.ok = TRUE), maxiter=as.integer(maxiter), as.integer(sum(object$P1inf)), 
                    convtol, as.integer(nnd), as.integer(nsim), epsplus, etaplus, aplus1, c2, object$tol, info = integer(1), as.integer(antithetics), 
                    as.integer(sim), nsim2, as.integer(nd), as.integer(length(nd)),diff=double(1))
    if (!is.finite(out$diff)){
      warning("Non-finite difference in approximation algoritm. Returning -Inf.")
      return(-.Machine$double.xmax^0.75)
    }
    if(out$maxiter==maxiter){
      warning(paste("Maximum number of iterations reached, 
                    the approximation algorithm did not converge. Latest difference was",out$diff))
    }
    
    
    
    # add the scaling factor from approximating model
    logLik <- out$lik +
      sum(sapply(1:p,function(i) 
        switch(object$distribution[i], 
               gaussian = 0, 
               poisson = sum(dpois(x = object$y[,i], lambda = exp(out$theta[, i])*object$u[,i], log = TRUE), na.rm = TRUE),
               binomial = sum(dbinom(x = object$y[, i], size = object$u[, i], prob = (exp(out$theta[, i])/(1 + exp(out$theta[, i]))), log = TRUE), na.rm = TRUE), 
               gamma = sum(dgamma(x = object$y[,i], shape = object$u[, i], scale = exp(out$theta[, i])/object$u[, i], log = TRUE), na.rm = TRUE), 
               `negative binomial` = sum(dnbinom(x = object$y[, i], size = object$u[, i], 
                                                 mu = exp(out$theta[, i]), log = TRUE), na.rm = TRUE))))
    
  }
  logLik
} 
