#' Simulation of a Gaussian State Space Model
#'
#' Function \code{simulateSMM} simulates states, signals, disturbances or missing observations of
#' the Gaussian state space model either conditional on the data (simulation smoother) or
#' unconditionally.
#'
#' Simulation smoother algorithm is based on article by J. Durbin and S.J. Koopman (2002).
#' The simulation filter (\code{filtered = TRUE}) is a straightforward modification
#' of the simulations smoother, where only filtering steps are performed.
#'
#' Function can use two antithetic variables, one for location and other for scale, so output
#' contains four blocks of simulated values which correlate which each other (ith block correlates
#' negatively with (i+1)th block, and positively with (i+2)th block etc.).
#'
#' Note that KFAS versions 1.2.0 and older, for unconditional simulation the initial
#' distribution of states was fixed so that \code{a1} was set to the smoothed estimates
#' of the first state and the initial variance was set to zero. Now original
#' \code{a1} and \code{P1} are used, and \code{P1inf} is ignored (i.e. diffuse states are
#' fixed to corresponding elements of \code{a1}).
#'
#' @export
#' @param object Gaussian state space object of class \code{SSModel}.
#' @param type What to simulate.
#' @param filtered Simulate from \eqn{p(\alpha_t|y_{t-1},...,y_1)}{p(\alpha[t]|y[t-1],...,y[1])}
#'   instead of \eqn{p(\alpha|y)}.
#' @param nsim Number of independent samples. Default is 1.
#' @param antithetics Use antithetic variables in simulation. Default is \code{FALSE}.
#' @param conditional Simulations are conditional to data. If \code{FALSE}, the
#' states having exact diffuse initial distribution (as defined by \code{P1inf}
#' are fixed to corresponding values of \code{a1}. See details.
#' @return An n x k x nsim array containing the simulated series, where k is number of observations,
#'   signals, states or disturbances.
#' @references Durbin J. and Koopman, S.J. (2002). A simple and efficient simulation smoother for
#'   state space time series analysis, Biometrika, Volume 89, Issue 3
#' @examples
#'
#' model <- SSModel(matrix(NA, 100, 1) ~ SSMtrend(1, 1, P1inf = 0), H = 1)
#'
#' set.seed(123)
#' sim <- simulateSSM(model, "obs", nsim = 2, antithetics = TRUE)
#' # first time points
#' sim[1,,]
#' # correlation structure between simulations with two antithetics
#' cor(sim[,1,])
#'
#' out_NA <- KFS(model, filtering = "none", smoothing = "state")
#' model["y"] <- sim[, 1, 1]
#' out_obs <- KFS(model, filtering = "none", smoothing = "state")
#'
#' set.seed(40216)
#' # simulate states from the p(alpha | y)
#' sim_conditional <- simulateSSM(model, nsim = 10, antithetics = TRUE)
#'
#' # mean of the simulated states is exactly correct due to antithetic variables
#' mean(sim_conditional[2, 1, ])
#' out_obs$alpha[2]
#' # for variances more simulations are needed
#' var(sim_conditional[2, 1, ])
#' out_obs$V[2]
#'
#' set.seed(40216)
#' # no data, simulations from p(alpha)
#' sim_unconditional <- simulateSSM(model, nsim = 10, antithetics = TRUE,
#'   conditional = FALSE)
#' mean(sim_unconditional[2, 1, ])
#' out_NA$alpha[2]
#' var(sim_unconditional[2, 1, ])
#' out_NA$V[2]
#'
#' ts.plot(cbind(sim_conditional[,1,1:5], sim_unconditional[,1,1:5]),
#'   col = rep(c(2,4), each = 5))
#' lines(out_obs$alpha, lwd=2)
#'
simulateSSM <- function(object,
  type = c("states", "signals", "disturbances", "observations", "epsilon", "eta"),
  filtered = FALSE, nsim = 1, antithetics = FALSE, conditional = TRUE) {

  # Check that the model object is of proper form
  is.SSModel(object, na.check = TRUE, return.logical = FALSE)
  sim.what <- match.arg(arg = type,
    choices = c("states", "signals", "disturbances", "observations", "epsilon", "eta"))
  if (any(object$distribution != "gaussian"))
    stop("Function is only for Gaussian models.")

  if (conditional && sim.what == "observations" && all(!is.na(object$y)))
    stop("There is no missing observations, nothing to simulate.")
  p <- attr(object, "p")
  if (sim.what == "observations") {
    object <- transformSSM(object, type = "augment")
  } else {
    htol <- max(100, max(apply(object$H, 3, diag))) * .Machine$double.eps
    if (p > 1 && any(abs(apply(object$H, 3, "[", !diag(p))) > htol))
      object <- transformSSM(object = object, type = "ldl")
  }
  m <- attr(object, "m")
  k <- attr(object, "k")
  n <- attr(object, "n")
  tv <- attr(object, "tv")
  ymiss <- is.na(object$y)
  storage.mode(ymiss) <- "integer"

  if (!conditional || all(is.na(object$y))) { # cannot simulate from N(a1,Inf)
    object$P1inf[] <- 0
  }
  simtmp <- simHelper(object, nsim, antithetics)

  sim.what <- which(c("epsilon", "eta", "disturbances", "states", "signals", "observations") == sim.what)
  simdim <- as.integer(switch(sim.what, p, k, p + k, m, p, p))

  if (!conditional || all(is.na(object$y))) {
    out <- .Fortran(fsimgaussianuncond, NAOK = TRUE, tv,
      object$Z, object$H, object$T, object$R, object$Q, object$a1, object$P1,
      object$P1inf, simtmp$nNonzeroP1, as.integer(nsim), simtmp$epsplus,
      simtmp$etaplus, simtmp$aplus1, p, n, m, k, info = as.integer(0),
      simtmp$nNonzeroP1inf, object$tol, simtmp$zeroP1inf, length(simtmp$zeroP1inf),
      sim = array({
        if (sim.what == 6) t(object$y) else 0
      }, c(simdim, n, 3 * nsim * antithetics + nsim)), simtmp$c2, sim.what,
      simdim, as.integer(antithetics))
  } else {
    if (!filtered) { #simulation smoother
      out <- .Fortran(fsimgaussian, NAOK = TRUE, ymiss, tv, object$y,
        object$Z, object$H, object$T, object$R, object$Q, object$a1, object$P1,
        object$P1inf, simtmp$nNonzeroP1, as.integer(nsim), simtmp$epsplus,
        simtmp$etaplus, simtmp$aplus1, p, n, m, k, info = as.integer(0),
        simtmp$nNonzeroP1inf, object$tol, simtmp$zeroP1inf, length(simtmp$zeroP1inf),
        sim = array({
          if (sim.what == 6) t(object$y) else 0
        }, c(simdim, n, 3 * nsim * antithetics + nsim)), simtmp$c2, sim.what,
        simdim, as.integer(antithetics))
    } else { # simulate from predictive distribution
      if (!(sim.what %in% (4:5)))
        stop("Only state and signal simulation filtering is supported.")
      out <- .Fortran(fsimfilter, NAOK = TRUE, ymiss, tv, object$y,
        object$Z, object$H, object$T, object$R, object$Q, object$a1, object$P1,
        object$P1inf, simtmp$nNonzeroP1, as.integer(nsim), simtmp$epsplus,
        simtmp$etaplus, simtmp$aplus1, p, n, m, k, info = as.integer(0),
        simtmp$nNonzeroP1inf, object$tol, simtmp$zeroP1inf, length(simtmp$zeroP1inf),
        sim = array(0, c(simdim, n, 3 * nsim * antithetics + nsim)),
        simtmp$c2, sim.what, simdim, as.integer(antithetics))
    }
  }
  if(out$info!=0){
    stop(switch(as.character(out$info),
      "-3" = "Couldn't compute LDL decomposition of P1.",
      "-2" =  "Couldn't compute LDL decomposition of Q."
    ))
  }
  rownames(out$sim) <- switch(sim.what, rep("eps", p), rep("eta", k),
    c(rep("eps", p), rep("eta", k)),
    rownames(object$a1), colnames(object$y), colnames(object$y))
  aperm(out$sim, c(2, 1, 3))
}
