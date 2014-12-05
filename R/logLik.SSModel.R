#' Log-likelihood of the State Space Model.
#'
#' Function \code{logLik.SSmodel} computes the log-likelihood value of a state 
#' space model.
#' 
#' As KFAS is based on diffuse initialization, the likelihood is also diffuse, which coincides with 
#' restricted likelihood (REML) in appropriate (mixed) models. However, in typical REML estimation 
#' constant term \eqn{log|X'X|} is omitted from the log-likelihood formula. Similar term is also 
#' missing in diffuse likelihood formulations, but unlike in simpler linear models this term is not
#' necessarily constant. Therefore omitting this term can lead to suboptimal results in model 
#' estimation if there is unknown parameters in diffuse parts of Zt or Tt (Francke et al. 2010). 
#' See also Gurka (2006) for model comparison in mixed model settings with and without the additional (constant) term.
#'  
#' Note that for non-Gaussian models with importance sampling derivative-free 
#' optimization methods such as Nelder-Mead might be more reliable than methods 
#' which use finite difference approximations. This is due to noise caused by 
#' the relative stopping criterion used for finding approximating Gaussian model.
#' 
#' @references Francke, M. K., Koopman, S. J. and De Vos, A. F. (2010), 
#' Likelihood functions for state space models with diffuse initial conditions. 
#' Journal of Time Series Analysis, 31: 407-414.\cr
#' 
#' Gurka, M. J (2006),
#' Selecting the Best Linear Mixed Model Under REML.
#' The American Statistician, Vol. 60, Iss. 1.
#' 
#' @export
#' @aliases logLik logLik.SSModel
#' @param object State space model of class \code{SSModel}.
#' @param marginal Logical. Compute marginal instead of diffuse likelihood (see details). 
#' Default is FALSE.
#' @param nsim Number of independent samples used in estimating the
#' log-likelihood of the non-Gaussian state space model. Default is 0, which
#' gives good starting value for optimization. Only used for non-Gaussian model.
#' @param antithetics Logical. If TRUE, two antithetic variables are used in
#' simulations, one for location and another for scale. Default is TRUE. Only used for non-Gaussian model.
#' @param theta Initial values for conditional mode theta. Only used for non-Gaussian model.
#' @param check.model Logical. If TRUE, function \code{is.SSModel} is called before computing the likelihood. 
#' Default is TRUE. Function is slightly faster if no checking is done, but note that if you manually 
#' modify model object improperly, underlying Fortran code can crash R.
#' @param transform How to transform the model in case of non-diagonal
#' covariance matrix \eqn{H}. Defaults to \code{'ldl'}. See function \code{\link{transformSSM}} for
#' details. 
#' @param maxiter The maximum number of iterations used in linearisation. Default is 50. Only used for non-Gaussian model.
#' @param seed The value is used as a seed via set.seed function. Only used for non-Gaussian model.
#' @param convtol Tolerance parameter for convergence checks for Gaussian approximation.
#'  Iterations are continued until the scaled norm between three successive iterations is smaller than \code{convtol}.
#' @param stepmax Maximum stepsize used in Gaussian approximation. Only used for non-Gaussian models.
#' @param ... Ignored.
#' @return \item{}{log-likelihood of the state space model.}
logLik.SSModel <- 
  function(object, marginal=FALSE, nsim = 0, antithetics = TRUE, theta, check.model = TRUE, 
           transform = c("ldl", "augment"), maxiter = 50, seed, convtol = 1e-8, stepmax,...) {
    # Check that the model object is of proper form
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
    tv <- attr(object, "tv")
    if (all(object$distribution == "gaussian")) {
      if (all(c(object$Q, object$H) == 0) || all(c(object$R, object$H) == 0) || 
            any(!is.finite(c(object$R, object$Q, object$H) == 0))) 
        return(-.Machine$double.xmax^0.75)     
      if (p == 1) {
        out <- .Fortran(fglogliku, NAOK = TRUE, object$y, ymiss, tv, 
                        object$Z, object$H, object$T, object$R, object$Q, object$a1, object$P1, 
                        object$P1inf, m, k, n, lik = double(1), 
                        object$tol, as.integer(sum(object$P1inf)),marginal=as.integer(marginal))
      } else {
        if (any(abs(apply(object$H, 3, "[", !diag(p))) > object$tol)) {
          object <- tryCatch(transformSSM(object, type = match.arg(arg = transform, 
                                                                   choices = c("ldl", "augment"))), error = function(e) e)
          if (!inherits(object, "SSModel")) {
            warning(object$message)
            return(-.Machine$double.xmax^0.75)
          }
          tv <- attr(object, "tv") 
        } 
        
        out <- .Fortran(fgloglik, NAOK = TRUE, object$y, ymiss, tv, 
                        object$Z, object$H, object$T, object$R, object$Q, object$a1, object$P1, 
                        object$P1inf, p, m, k, n, 
                        lik = double(1), object$tol, as.integer(sum(object$P1inf)),marginal=as.integer(marginal))
        
      }
      if(out$marginal==-1){        
        warning("Computation of marginal likelihood failed, could not compute the additional term.")
        return(-.Machine$double.xmax^0.75)
      }      
    } else {
      if(maxiter<1)
        stop("Argument maxiter must a positive integer. ")
      if (all(c(object$Q, object$u) == 0) || all(c(object$R, object$u) == 0) || 
            any(!is.finite(c(object$R, object$Q, object$u) == 0))) 
        return(-.Machine$double.xmax^0.75)
      if (missing(theta)) {
        theta <- init_theta(object$y, object$u, object$distribution)
      } else theta <- array(theta, dim = c(n, p))
      if(missing(stepmax))
        stepmax<-max(1,sqrt(sum(theta^2)))
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
        x2 <- array(apply(object$Q, 3, diag) > object$tol, c(k, (n - 1) * tv[5] + 
                                                               1))
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
      out <- .Fortran(fngloglik, NAOK = TRUE, object$y, ymiss, tv, 
                      object$Z, object$T, object$R, object$Q, object$a1, object$P1, object$P1inf, 
                      p, m, k, n, lik = double(1), 
                      theta = theta, object$u, pmatch(x = object$distribution, table = c("gaussian", 
                                                                                         "poisson", "binomial", "gamma", "negative binomial"), duplicates.ok = TRUE), 
                      maxiter = as.integer(maxiter), as.integer(sum(object$P1inf)), convtol, 
                      as.integer(nnd), as.integer(nsim), epsplus, etaplus, aplus1, c2, object$tol, 
                      info = integer(1), as.integer(antithetics), as.integer(sim), nsim2, as.integer(nd), 
                      as.integer(length(nd)), diff = double(1),stepmax=as.double(stepmax),marginal=as.integer(marginal))
      if(out$info!=0){
        if (out$info==1)
          warning("Non-finite value of likelihood or linear predictor in approximation algorithm.")
        if(out$info==2)
          warning(paste("Maximum number of iterations reached, latest difference was", signif(out$diff,3)))
        if(out$info==5)
          warning("Computation of marginal likelihood failed, could not compute the additional term.")
        if (out$info == -2)
          warning("Couldn't compute LDL decomposition of Q.")  
        if (out$info == -3)
          warning("Couldn't compute LDL decomposition of P1.")
        if(out$info!=2)
          return(-.Machine$double.xmax^0.75)  
      }    
    } 
    out$lik
  }

