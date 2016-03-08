#' Log-likelihood of the State Space Model.
#'
#' Function \code{logLik.SSmodel} computes the log-likelihood value of a state
#' space model.
#'
#' As KFAS is based on diffuse initialization, the likelihood is also diffuse,
#' which coincides with restricted likelihood (REML) in appropriate (mixed)
#' models. However, in typical REML estimation constant term \eqn{log|X'X|} is
#' omitted from the log-likelihood formula. Similar term is also missing in
#' diffuse likelihood formulations, but unlike in simpler linear models this
#' term is not necessarily constant. Therefore omitting this term can lead to
#' suboptimal results in model estimation if there is unknown parameters in
#' diffuse parts of Zt or Tt (Francke et al. 2010). See also Gurka (2006) for
#' model comparison in mixed model settings with and without the additional
#' (constant) term (for BIC it could be better to use marginal likelihood instead
#' of diffuse likelihood and vice versa for AIC) and Casals et al. (2014).
#' The marginal likelihood can be computed by setting \code{marginal = TRUE}.
#'
#' Note that for non-Gaussian models with importance sampling derivative-free
#' optimization methods such as Nelder-Mead might be more reliable than methods
#' which use finite difference approximations. This is due to noise caused by
#' the relative stopping criterion used for finding approximating Gaussian
#' model.
#'
#' @references Francke, M. K., Koopman, S. J. and De Vos, A. F. (2010),
#'   Likelihood functions for state space models with diffuse initial
#'   conditions. Journal of Time Series Analysis, 31: 407--414.\cr
#'
#'   Gurka, M. J (2006), Selecting the Best Linear Mixed Model Under REML. The
#'   American Statistician, Vol. 60.\cr
#'
#'   Casals, J., Sotoca, S., Jerez, M. (2014), Minimally conditioned likelihood
#'   for a nonstationary state space model. Mathematics and Computers in
#'   Simulation, Vol. 100.
#'
#' @export
#' @importFrom stats logLik
#' @aliases logLik logLik.SSModel
#' @param object State space model of class \code{SSModel}.
#' @param marginal Logical. Compute marginal instead of diffuse likelihood (see
#'   details). Default is \code{FALSE}.
#' @param nsim Number of independent samples used in estimating the
#'   log-likelihood of the non-Gaussian state space model. Default is 0, which
#'   gives good starting value for optimization. Only used for non-Gaussian
#'   model.
#' @param antithetics Logical. If TRUE, two antithetic variables are used in
#'   simulations, one for location and another for scale. Default is TRUE. Only
#'   used for non-Gaussian model.
#' @param theta Initial values for conditional mode theta. Only used for
#'   non-Gaussian model.
#' @param check.model Logical. If TRUE, function \code{is.SSModel} is called
#'   before computing the likelihood. Default is \code{TRUE}.
#' @param transform How to transform the model in case of non-diagonal
#'   covariance matrix \eqn{H}. Defaults to \code{"ldl"}. See function
#'   \code{\link{transformSSM}} for details.
#' @param maxiter The maximum number of iterations used in linearisation.
#'   Default is 50. Only used for non-Gaussian model.
#' @param seed The value is used as a seed via \code{set.seed} function. Only used for
#'   non-Gaussian model.
#' @param convtol Tolerance parameter for convergence checks for Gaussian
#'   approximation.
#' @param ... Ignored.
#' @return Log-likelihood of the model.
logLik.SSModel <- function(object, marginal=FALSE, nsim = 0,
  antithetics = TRUE, theta, check.model = TRUE,
  transform = c("ldl", "augment"), maxiter = 50, seed, convtol = 1e-8, ...) {

  # Check that the model object is of proper form
  if (check.model) {
    if (!is.SSModel(object, na.check = TRUE)) {
      return(-.Machine$double.xmax ^ 0.75)
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
    # degenerate case
    if (all(c(object$Q, object$H) < .Machine$double.eps) || all(c(object$R, object$H) < .Machine$double.eps))
      return(-.Machine$double.xmax ^ 0.75)
    htol <- max(100, max(apply(object$H, 3, diag))) * .Machine$double.eps
    if (p > 1 && any(abs(apply(object$H, 3, "[", !diag(p))) > htol)) {
      object <-
        tryCatch(transformSSM(object, type = match.arg(arg = transform, choices = c("ldl", "augment"))),
          error = function(e) e)
      if (!inherits(object, "SSModel")) {
        warning(object$message)
        return(-.Machine$double.xmax ^ 0.75)
      }
      m <- attr(object, "m")
      k <- attr(object, "k")
      tv <- attr(object, "tv")
    }
    out <- .Fortran(fgloglik, NAOK = TRUE, t(object$y), t(ymiss), tv,
      aperm(object$Z,c(2,1,3)), object$H, object$T, object$R, object$Q, object$a1, object$P1,
      object$P1inf, p, m, k, n,
      lik = double(1), object$tol, as.integer(sum(object$P1inf)),marginal=as.integer(marginal))


    if (out$marginal == -1) {
      warning("Computation of marginal likelihood failed, could not compute the additional term.")
      return(-.Machine$double.xmax ^ 0.75)
    }
  } else {
    if (maxiter < 1)
      stop("Argument maxiter must a positive integer. ")
    if (all(c(object$Q, object$u) == 0) || all(c(object$R, object$u) == 0) ||
        any(!is.finite(c(object$R, object$Q, object$u) == 0)))
      return(-.Machine$double.xmax ^ 0.75)
    if (missing(theta) || is.null(theta)) {
      theta <- initTheta(object$y, object$u, object$distribution)
    } else theta <- array(theta, dim = c(n, p))
    if (nsim == 0) {
      nsim <- 1
      sim <- 0
      simtmp <- list(epsplus = array(0, c(1, 1, 1)), etaplus = array(0, c(1, 1, 1)),
        aplus1 = array(0, dim = c(1, 1)), c2 = numeric(1),
        nonzeroP1 = which(diag(object$P1) > object$tol),
        nNonzeroP1 = length(which(diag(object$P1) > object$tol)),
        zeroP1inf = which(diag(object$P1inf) > 0),
        nNonzeroP1inf = as.integer(sum(object$P1inf)))
    } else {
      sim <- 1
      if (missing(seed))
        seed <- 123
      set.seed(seed)
      simtmp <- simHelper(object, nsim, antithetics)
    }
    nsim2 <- as.integer(max(sim * (3 * antithetics * nsim + nsim), 1))
    out <- .Fortran(fngloglik, NAOK = TRUE, object$y, ymiss, tv,
      object$Z, object$T, object$R, object$Q, object$a1, object$P1, object$P1inf,
      p, m, k, n, lik = double(1),
      theta = theta, object$u,
      pmatch(x = object$distribution,
        table = c("gaussian", "poisson", "binomial", "gamma", "negative binomial"),
        duplicates.ok = TRUE),
      maxiter = as.integer(maxiter), simtmp$nNonzeroP1inf, convtol,
      simtmp$nNonzeroP1, as.integer(nsim), simtmp$epsplus, simtmp$etaplus,
      simtmp$aplus1, simtmp$c2, object$tol,
      info = integer(1), as.integer(antithetics), as.integer(sim), nsim2,
      simtmp$zeroP1inf, length(simtmp$zeroP1inf), diff = double(1),
      marginal = as.integer(marginal))

    if(out$info!=0){
      warning(switch(as.character(out$info),
        "-3" = "Couldn't compute LDL decomposition of P1.",
        "-2" =  "Couldn't compute LDL decomposition of Q.",
        "1" = "Gaussian approximation failed due to non-finite value in linear predictor.",
        "2" = "Gaussian approximation failed due to non-finite value of p(theta|y).",
        "3" = "Maximum number of iterations reached, the approximation did not converge.",
        "5" = "Computation of marginal likelihood failed, could not compute the additional term."
      ))
      if(out$info!=3) return(-.Machine$double.xmax^0.75)
    }
  }
  out$lik
}
