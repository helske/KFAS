#' Importance Sampling of Exponential Family State Space Model
#'
#' Importance Sampling of Exponential Family State Space Model.
#'
#' Function \code{importanceSSM} simulates states or signals of the exponential family state space 
#' model conditioned with the observations, returning the simulated samples of the states/signals 
#' with the corresponding importance weights.
#'
#' Function can use two antithetic variables, one for location and other for
#' scale, so output contains four blocks of simulated values which correlate
#' which each other (ith block correlates negatively with (i+1)th block, and
#' positively with (i+2)th block etc.).
#' 
#' @export
#' @param model Exponential family state space model of class \code{SSModel}.
#' @param type What to simulate, \code{'states'} or \code{'signals'}. Default is \code{'states'}
#'@param filtered Simulate from \eqn{p(\alpha_t|y_{t-1},...,y_1)} instead of \eqn{p(\alpha|y)}.
#' @param nsim Number of independent samples. Default is 1000.
#' @param save.model Return the original model with the samples. Default is FALSE.
#' @param theta Initial values for conditional mode theta.
#' @param antithetics Logical. If TRUE, two antithetic variables are used in
#' simulations, one for location and another for scale. Default is FALSE.
#' @param maxiter Maximum number of iterations used in linearisation. Default is 50.
#' @param stepmax Maximum stepsize used in Gaussian approximation. Only used for non-Gaussian models.
#' @return A list containing elements \code{samples}, \code{weights} and \code{model} (if \code{save.model==TRUE}).
importanceSSM <- 
  function(model, type = c("states", "signals"), filtered = FALSE, 
    nsim = 1000, save.model = FALSE, theta, antithetics = FALSE, maxiter = 50,stepmax) {
    if(maxiter<1)
      stop("Argument maxiter must a positive integer. ")
    # Check that the model object is of proper form
    is.SSModel(model, na.check = TRUE, return.logical = FALSE)
    sim.what <- match.arg(arg = type, choices = c("states", "signals"))
    p <- attr(model, "p")
    m <- attr(model, "m")
    k <- attr(model, "k")
    n <- attr(model, "n")
    tv <- attr(model, "tv")
    ymiss <- is.na(model$y)
    storage.mode(ymiss) <- "integer"
    # initial values for linear predictor theta
    if (missing(theta)) {
        theta <- init_theta(model$y, model$u, model$distribution)
    } else theta <- array(theta, dim = c(n, p))
    if(missing(stepmax))
      stepmax<-max(1,sqrt(sum(theta^2)))
    epsplus <- array(0, c(p, n, nsim))
    etaplus <- array(0, c(k, n, nsim))
    aplus1 <- array(0, dim = c(m, nsim))
    x <- array(t(!ymiss), c(p, n, nsim))
    dfeps <- sum(x)/nsim
    x2 <- array(abs(apply(model$Q, 3, diag)) > model$tol, c(k, (n - 1) * tv[5] + 
        1))
    x2 <- array(x2, c(k, n, nsim))
    dfeta <- sum(x2)/nsim
    nonzero_P1 <- which(diag(model$P1) > model$tol)
    N_nonzero_P1 <- length(nonzero_P1)
    zero_P1inf <- which(diag(model$P1inf) == 0)
    dfu <- dfeps + dfeta + N_nonzero_P1
    u <- rnorm(dfu * nsim, mean = 0, sd = 1)
    if (dfeps > 0) 
        epsplus[x] <- u[1:(dfeps * nsim)]
    if (dfeta > 0) 
        etaplus[x2] <- u[(dfeps * nsim + 1):(dfeps * nsim + dfeta * nsim)]
    if (N_nonzero_P1 > 0) 
        aplus1[nonzero_P1, ] <- u[(dfeps * nsim + dfeta * nsim + 1):(dfu * nsim)]
    c2 <- numeric(nsim)
    if (antithetics) {
        for (i in 1:nsim) {
            u <- c(etaplus[, , i], epsplus[, , i], aplus1[, i])
            c2[i] <- t(u) %*% c(u)
        }
        q <- pchisq(c2, df = dfu)
        c2 <- sqrt(qchisq(1 - q, dfu)/c2)
    }
    sim.what <- which(c("epsilon", "eta", "disturbances", "states", "signals", "observations") == 
        sim.what)
    simdim <- as.integer(switch(sim.what, p, k, p + k, m, p, p))
    if (!filtered) {
        out <- .Fortran(fisample, NAOK = TRUE, model$y, ymiss, tv, model$Z, 
            model$T, model$R, model$Q, model$a1, model$P1, model$P1inf, model$u, 
            dist = pmatch(x = model$distribution, table = c("gaussian", "poisson", 
                "binomial", "gamma", "negative binomial"), duplicates.ok = TRUE), 
            p, n, m, k, theta, maxiter = as.integer(maxiter), 
            as.integer(sum(model$P1inf)), 1e-08, as.integer(N_nonzero_P1), as.integer(nsim), 
            epsplus, etaplus, aplus1, c2, model$tol, info = integer(1), as.integer(antithetics), 
            w = numeric(3 * nsim * antithetics + nsim), sim = array(0, c(simdim, 
                n, 3 * nsim * antithetics + nsim)), as.integer(zero_P1inf), as.integer(length(zero_P1inf)), 
            sim.what, simdim,stepmax=as.double(stepmax))
    } else {
        # warning('Filtered samples of non-Gaussian models is at a prototype stage,
        # results can be meaningless.')
        out <- .Fortran(fisamplefilter, NAOK = TRUE, model$y, ymiss, as.integer(tv), 
            model$Z, model$T, model$R, model$Q, model$a1, model$P1, model$P1inf, 
            model$u, dist = pmatch(x = model$distribution, table = c("gaussian", 
                "poisson", "binomial", "gamma", "negative binomial"), duplicates.ok = TRUE), 
            p, n, m, k, theta, maxiter = as.integer(maxiter), 
            as.integer(sum(model$P1inf)), 1e-08, as.integer(N_nonzero_P1), as.integer(nsim), 
            epsplus, etaplus, aplus1, c2, model$tol, info = integer(1), as.integer(antithetics), 
            w = array(0, c(n, 3 * nsim * antithetics + nsim)), sim = array(0, c(simdim, 
                n, 3 * nsim * antithetics + nsim)), as.integer(zero_P1inf), as.integer(length(zero_P1inf)), 
            sim.what, simdim,stepmax=as.double(stepmax))
    }
    if (maxiter == out$maxiter) 
        warning("Maximum number of iterations reached, the linearization did not converge.")
    out <- list(samples = aperm(out$sim, c(2, 1, 3)), weights = out$w)
    if (save.model) 
        out$model <- model
    # class(out) <- 'importanceSSM'
    out
} 
