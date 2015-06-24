#' Simulation of a Gaussian State Space Model
#' 
#' Function \code{simulateSMM} simulates states, signals, disturbances or missing observations of 
#' the Gaussian state space model.
#' 
#' 
#' Simulation smoother algorithm is based to article by J. Durbin and S.J. Koopman (2002).
#' 
#' Function can use two antithetic variables, one for location and other for scale, so output
#' contains four blocks of simulated values which correlate which each other (ith block correlates
#' negatively with (i+1)th block, and positively with (i+2)th block etc.).
#' 
#' @export
#' @param object Gaussian state space object of class \code{SSModel}.
#' @param type What to simulate.
#' @param filtered Simulate from \eqn{p(\alpha_t|y_{t-1},...,y_1)}{p(\alpha[t]|y[t-1],...,y[1])}
#'   instead of \eqn{p(\alpha|y)}.
#' @param nsim Number of independent samples. Default is 1.
#' @param antithetics Use antithetic variables in simulation. Default is FALSE.
#' @param conditional Simulations are conditional to data. If FALSE, the initial state
#'   \eqn{\alpha_1}{\alpha[1]} is set to \code{alphahat[1,]} computed by \code{KFS}, 
#'   and all the observations are removed from the model. Default is TRUE.
#' @return An n x k x nsim array containing the simulated series, where k is number of observations,
#'   signals, states or disturbances.
#' @references Durbin J. and Koopman, S.J. (2002). A simple and efficient simulation smoother for
#'   state space time series analysis, Biometrika, Volume 89, Issue 3
simulateSSM <- function(object, 
  type = c("states", "signals", "disturbances", "observations", "epsilon", "eta"), 
  filtered = FALSE, nsim = 1, antithetics = FALSE, conditional = TRUE) {
  
  # Check that the model object is of proper form
  is.SSModel(object, na.check = TRUE, return.logical = FALSE)
  sim.what <- match.arg(arg = type, 
    choices = c("states", "signals", "disturbances", "observations", "epsilon", "eta"))
  if (any(object$distribution != "gaussian")) 
    stop("Function is only for gaussian models.")
  if (!conditional) {
    out <- KFS(object, smoothing = "state")
    object$y[] <- NA
    object$a1[] <- out$alphahat[1, ]
    object$P1inf[] <- 0
    object$P1[] <- 0
  }
  if (sim.what == "observations" && sum(is.na(object$y)) == 0) 
    stop("There is no missing observations, nothing to simulate.")
  p <- attr(object, "p")
  if (sim.what == "observations") {
    object <- transformSSM(object, type = "augment")
  } else {
    if (p > 1 && any(abs(apply(object$H, 3, "[", !diag(p))) > object$tol)) 
      object <- transformSSM(object = object, type = "ldl")
  }
  m <- attr(object, "m")
  k <- attr(object, "k")
  n <- attr(object, "n")
  tv <- attr(object, "tv") 
  ymiss <- is.na(object$y)
  storage.mode(ymiss) <- "integer"    
  
  simtmp <- simHelper(object, ymiss, nsim, antithetics)
  
  sim.what <- which(c("epsilon", "eta", "disturbances", "states", "signals", "observations") == sim.what)
  simdim <- as.integer(switch(sim.what, p, k, p + k, m, p, p))
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
