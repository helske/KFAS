#' Importance Sampling of Exponential Family State Space Model
#'
#' Function \code{importanceSSM} simulates states or signals of the exponential
#' family state space model conditioned with the observations, returning the
#' simulated samples of the states/signals with the corresponding importance
#' weights.
#'
#' Function can use two antithetic variables, one for location and other for
#' scale, so output contains four blocks of simulated values which correlate
#' which each other (ith block correlates negatively with (i+1)th block, and
#' positively with (i+2)th block etc.).
#'
#' @export
#' @param model Exponential family state space model of class \code{SSModel}.
#' @param type What to simulate, \code{"states"} or \code{"signals"}. Default is
#'   \code{"states"}
#' @param filtered Simulate from \eqn{p(\alpha_t|y_{t-1},...,y_1)} instead of
#'   \eqn{p(\alpha|y)}. Note that for large models this can be very slow. Default is FALSE.
#' @param nsim Number of independent samples. Default is 1000.
#' @param save.model Return the original model with the samples. Default is
#'   FALSE.
#' @param theta Initial values for the conditional mode theta.
#' @param antithetics Logical. If TRUE, two antithetic variables are used in
#'   simulations, one for location and another for scale. Default is FALSE.
#' @param maxiter Maximum number of iterations used in linearisation. Default is
#'   50.
#' @return A list containing elements
#' \item{samples}{Simulated samples. }
#' \item{weights}{Importance weights. }
#' \item{model}{Original model in case of \code{save.model==TRUE}.}
#' @examples
#' data("sexratio")
#' model <- SSModel(Male ~ SSMtrend(1, Q = list(NA)), u = sexratio[,"Total"], data = sexratio,
#'                 distribution = "binomial")
#' fit <- fitSSM(model, inits = -15, method = "BFGS")
#' fit$model$Q #1.107652e-06

#' # Computing confidence intervals for sex ratio
#' # Uses importance sampling on response scale (1000 samples with antithetics)
#' set.seed(1)
#' imp <- importanceSSM(fit$model, nsim = 250, antithetics = TRUE)
#' sexratio.smooth <- numeric(length(model$y))
#' sexratio.ci <- matrix(0, length(model$y), 2)
#' w <- imp$w/sum(imp$w)
#' for(i in 1:length(model$y)){
#'   sexr <- exp(imp$sample[i,1,])
#'   sexratio.smooth[i]<-sum(sexr*w)
#'   oo <- order(sexr)
#'   sexratio.ci[i,] <- c(sexr[oo][which.min(abs(cumsum(w[oo]) - 0.05))],
#'                    sexr[oo][which.min(abs(cumsum(w[oo]) - 0.95))])
#' }
#'
#' \dontrun{
#' # Filtered estimates
#' impf <- importanceSSM(fit$model, nsim = 250, antithetics = TRUE,filtered=TRUE)
#' sexratio.filter <- rep(NA,length(model$y))
#' sexratio.fci <- matrix(NA, length(model$y), 2)
#' w <- impf$w/rowSums(impf$w)
#' for(i in 2:length(model$y)){
#'   sexr <- exp(impf$sample[i,1,])
#'   sexratio.filter[i] <- sum(sexr*w[i,])
#'   oo<-order(sexr)
#'   sexratio.fci[i,] <- c(sexr[oo][which.min(abs(cumsum(w[i,oo]) - 0.05))],
#'                     sexr[oo][which.min(abs(cumsum(w[i,oo]) - 0.95))])
#' }
#'
#' ts.plot(cbind(sexratio.smooth,sexratio.ci,sexratio.filter,sexratio.fci),
#'         col=c(1,1,1,2,2,2),lty=c(1,2,2,1,2,2))
#' }
importanceSSM <-  function(model, type = c("states", "signals"),
  filtered = FALSE,  nsim = 1000, save.model = FALSE, theta,
  antithetics = FALSE, maxiter = 50) {

  if (all(model$distribution == "gaussian")) {
    stop("Model is completely Gaussian, use simulateSSM instead. ")
  }
  if(maxiter<1)
    stop("Argument maxiter must a positive integer. ")
  if(nsim<1)
    stop("Argument nsim must a positive integer. ")
  # Check that the model object is of proper form
  is.SSModel(model, na.check = TRUE, return.logical = FALSE)


  p <- attr(model, "p")
  m <- attr(model, "m")
  k <- attr(model, "k")
  n <- attr(model, "n")
  tv <- attr(model, "tv")
  ymiss <- is.na(model$y)
  storage.mode(ymiss) <- "integer"

  # initial values for linear predictor theta
  if (missing(theta)) {
    theta <- initTheta(model$y, model$u, model$distribution)
  } else theta <- array(theta, dim = c(n, p))

  # generate standard normal variables for importance sampling
  simtmp <- simHelper(model, nsim, antithetics)
  sim.what <- which(c("epsilon", "eta", "disturbances", "states", "signals", "observations") ==
      match.arg(arg = type, choices = c("states", "signals")))
  simdim <- as.integer(switch(sim.what, p, k, p + k, m, p, p))
  if (!filtered) {
    out <- .Fortran(fisample, NAOK = TRUE, model$y, ymiss, tv, model$Z,
      model$T, model$R, model$Q, model$a1, model$P1, model$P1inf, model$u,
      dist = pmatch(x = model$distribution,
        table = c("gaussian", "poisson",  "binomial", "gamma", "negative binomial"),
        duplicates.ok = TRUE),
      p, n, m, k, theta, maxiter = as.integer(maxiter),
      simtmp$nNonzeroP1inf, 1e-08, simtmp$nNonzeroP1, as.integer(nsim),
      simtmp$epsplus, simtmp$etaplus, simtmp$aplus1, simtmp$c2, model$tol,
      info = integer(1), as.integer(antithetics),
      w = numeric(3 * nsim * antithetics + nsim),
      sim = array(0, c(simdim,  n, 3 * nsim * antithetics + nsim)), simtmp$zeroP1inf,
      length(simtmp$zeroP1inf), sim.what, simdim)
  } else {
    out <- .Fortran(fisamplefilter, NAOK = TRUE, model$y, ymiss, as.integer(tv),
      model$Z, model$T, model$R, model$Q, model$a1, model$P1, model$P1inf,model$u,
      dist = pmatch(x = model$distribution,
        table = c("gaussian",  "poisson", "binomial", "gamma", "negative binomial"),
        duplicates.ok = TRUE),
      p, n, m, k, theta, maxiter = as.integer(maxiter),
      simtmp$nNonzeroP1inf, 1e-08, simtmp$nNonzeroP1, as.integer(nsim),
      simtmp$epsplus, simtmp$etaplus, simtmp$aplus1, simtmp$c2,
      model$tol, info = integer(1), as.integer(antithetics),
      w = array(0, c(n, 3 * nsim * antithetics + nsim)),
      sim = array(0, c(simdim, n, 3 * nsim * antithetics + nsim)), simtmp$zeroP1inf,
      length(simtmp$zeroP1inf), sim.what, simdim)
  }
  if(out$info!=0){
    switch(as.character(out$info),
      "-3" = stop("Couldn't compute LDL decomposition of P1."),
      "-2" =  stop("Couldn't compute LDL decomposition of Q."),
      "1" = stop("Gaussian approximation failed due to non-finite value in linear predictor."),
      "2" = stop("Gaussian approximation failed due to non-finite value of p(theta|y)."),
      "3" = warning("Maximum number of iterations reached, the approximation did not converge.")
    )
  }

  out <- list(samples = aperm(out$sim, c(2, 1, 3)), weights = out$w)
  if (save.model)
    out$model <- model
  out
}
