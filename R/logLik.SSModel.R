#' Log-likelihood of the State Space Model.
#'
#' Function \code{logLik.SSmodel} computes the log-likelihood value of a state
#' space model.
#'
#' As KFAS is based on diffuse initialization, the log-likelihood is also diffuse,
#' which coincides with restricted likelihood (REML) in an appropriate (mixed)
#' models. However, in typical REML estimation constant term \eqn{log|X'X|} is
#' omitted from the log-likelihood formula. Similar term is also missing in
#' diffuse log-likelihood formulations of state space models, but unlike in simpler
#' linear models this term is not necessarily constant. Therefore omitting this
#' term can lead to suboptimal results in model estimation if there is unknown
#' parameters in diffuse parts of Zt or Tt (Francke et al. 2011). Therefore
#' so called marginal log-likelihood (diffuse likelihood + extra term) is
#' recommended. See also Gurka (2006) for model comparison in mixed model
#' settings using REML with and without the additional (constant) term.
#' The marginal log-likelihood can be computed by setting \code{marginal = TRUE}.
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
#' @param expected Logical value defining the approximation of H_t in case of Gamma 
#' and negative binomial distribution. Default is \code{FALSE} which matches the 
#' algorithm of Durbin & Koopman (1997), whereas \code{TRUE} uses the expected value
#' of observations in the equations, leading to results which match with \code{glm} (where applicable).
#' The latter case was the default behaviour of KFAS before version 1.3.8.
#' Essentially this is the difference between observed and expected information.
#' @param ... Ignored.
#' @return Log-likelihood of the model.
#' @examples 

#' # Example of estimating AR model with covariates, and how to deal with possible
#' # non-stationarity in optimization.
#' 
#' set.seed(1)
#' x <- rnorm(100)
#' y <- 2 * x + arima.sim(n = 100, model = list(ar = c(0.5, -0.3)))
#' 
#' model<- SSModel(y ~ SSMarima(ar = c(0.5, -0.3),  Q = 1) + x, H = 0)
#' 
#' obj <- function(pars, model, estimate = TRUE) {
#'   #guard against stationarity
#'   armamod <- try(SSMarima(ar = artransform(pars[1:2]), Q = exp(pars[3])), silent = TRUE)
#'   if(class(armamod) == "try-error") {
#'     return(Inf)
#'   } else {
#'     # use advanced subsetting method for SSModels, see ?`[.SSModel`
#'     model["T", states = "arima"] <- armamod$T
#'     model["Q", eta = "arima"]  <- armamod$Q
#'     model["P1", states = "arima"]  <- armamod$P1
#'     if(estimate) {
#'       -logLik(model)
#'     } else {
#'       model
#'     }
#'   }
#' }
#' fit <- optim(p = c(0.5,-0.5,1), fn = obj, model = model, method ="BFGS")
#' 
#' model <- obj(fit$par, model, FALSE)
#' model$T
#' model$Q
#' coef(KFS(model), last = TRUE)
#' 
logLik.SSModel <- function(object, marginal=FALSE, nsim = 0,
  antithetics = TRUE, theta, check.model = TRUE,
  transform = c("ldl", "augment"), maxiter = 50, seed, convtol = 1e-8, 
  expected = FALSE, ...) {

  # Check that the model object is of proper form
  if (check.model) {
    if (!is.SSModel(object, na.check = TRUE)) {
      return(-.Machine$double.xmax ^ 0.75)
    }
  }
  if (!is.logical(expected))
    stop("Argument expected should be logical. ")
  expected <- as.integer(expected)
  p <- attr(object, "p")
  m <- attr(object, "m")
  k <- attr(object, "k")
  n <- attr(object, "n")
  ymiss <- array(is.na(object$y), dim = c(n, p))
  storage.mode(ymiss) <- "integer"
  tv <- attr(object, "tv")
  if (all(object$distribution == "gaussian")) {
    # degenerate case
    if (all(c(object$Q, object$H) < .Machine$double.eps^0.75) || all(c(object$R, object$H) < .Machine$double.eps^0.75))
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
      lik = double(1), object$tol, as.integer(sum(object$P1inf)))
    
    if (marginal) {
      logdetxx <- .Fortran(fmarginalxx, NAOK = TRUE, object$P1inf, object$Z, object$T, m, p, n, 
        as.integer(sum(object$P1inf)), tv, lik = double(1), info = integer(1))
      
      if (logdetxx$info == -1) {
        warning("Computation of marginal likelihood failed, could not compute the additional term.")
        return(-.Machine$double.xmax ^ 0.75)
      } else {
        out$lik <- out$lik + logdetxx$lik
      }
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
      diff = double(1),
      marginal = as.integer(marginal), expected)

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
