#' Kalman Filter and Smoother with Exact Diffuse Initialization for Exponential
#' Family State Space Models
#'
#' Performs Kalman filtering and smoothing with exact diffuse initialization
#' using univariate approach for exponential family state space models.
#'
#' Notice that in case of multivariate Gaussian observations, \code{v}, \code{F},
#' \code{Finf}, \code{K} and \code{Kinf} are usually not the same as those
#' calculated in usual multivariate Kalman filter. As filtering is done one
#' observation element at the time, the elements of the prediction error
#' \eqn{v_t}{v[t]} are uncorrelated, and \code{F}, \code{Finf}, \code{K} and
#' \code{Kinf} contain only the diagonal elemens of the corresponding covariance
#' matrices. The usual multivariate versions of \code{F} and \code{v} can be
#' obtained from the output of \code{KFS} using the function
#' \code{\link{mvInnovations}}.
#'
#' In rare cases of a diffuse initialization phase with highly correlated
#' states, cumulative rounding errors in computing \code{Finf} and \code{Pinf}
#' can sometimes cause the diffuse phase end too early,
#' or the backward smoothing gives negative variances. Changing the tolerance
#' parameter \code{tol} of the model (see \code{\link{SSModel}}) to smaller (or
#' larger) can sometimes help. Another option is to redefine the prior state
#' variances more informative.
#'
#' Fon non-Gaussian models the components corresponding to diffuse filtering
#' (\code{Finf}, \code{Pinf}, \code{d}, \code{Kinf}) are not returned even
#' when \code{filtering} is used. Results based on approximating Gaussian model
#' can be obtained by running \code{KFS} using the output of \code{approxSSM}.
#'
#' In case of non-Gaussian models with \code{nsim = 0}, the smoothed estimates
#' relate to the conditional mode of \eqn{p(\alpha|y)}. When using importance
#' sampling (\code{nsim>0}), results correspond to the conditional mean of
#'  \eqn{p(\alpha|y)}.
#'
#' @export
#' @importFrom stats start frequency tsp<- tsp ts
#' @param model Object of class \code{SSModel}.
#' @param filtering Types of filtering. Possible choices are \code{"state"},
#'   \code{"signal"}, \code{"mean"}, and \code{"none"}. Default is
#'   \code{"state"} for Gaussian and \code{"none"} for non-Gaussian models.
#'   Multiple values are allowed. For Gaussian models, the signal is the mean.
#'   Note that filtering for non-Gaussian models with importance sampling can be
#'   very slow with large models.
#' @param smoothing Types of smoothing. Possible choices are \code{"state"},
#'   \code{"signal"}, \code{"mean"}, \code{"disturbance"}, and \code{"none"}. Default is \code{"state"} and \code{"mean"}. For
#'   non-Gaussian models, option \code{"disturbance"} is not supported, and for
#'   Gaussian models option \code{"mean"} is identical to \code{"signal"}. Multiple values are
#'   allowed.
#' @param simplify If \code{FALSE} and the model is completely Gaussian, \code{KFS} returns some
#' generally not so interesting variables from filtering and smoothing. Default
#' is \code{TRUE}.
#' @param transform How to transform the model in case of non-diagonal
#'   covariance matrix \code{H}. Defaults to \code{"ldl"}. See function
#'   \code{\link{transformSSM}} for details.
#' @param nsim The number of independent samples used in importance sampling.
#' Only used for non-Gaussian models. Default is 0, which computes the
#' approximating Gaussian model by \code{\link{approxSSM}} and performs the
#' usual Gaussian filtering/smoothing so that the smoothed state estimates
#' equals to the conditional mode of \eqn{p(\alpha_t|y)}{p(\alpha[t]|y)}.
#' In case of \code{nsim = 0}, the mean estimates and their variances are computed using
#' the Delta method (ignoring the covariance terms).
#' @param theta Initial values for conditional mode theta. Only used for
#'   non-Gaussian models.
#' @param maxiter The maximum number of iterations used in Gaussian
#'   approximation. Default is 50. Only used for non-Gaussian models.
#' @param convtol Tolerance parameter for convergence checks for Gaussian
#'   approximation. Only used for non-Gaussian models.
#' @return What \code{KFS} returns depends on the arguments \code{filtering},
#'   \code{smoothing} and \code{simplify}, and whether the model is Gaussian or
#'   not:
#'
#'   \item{model}{Original state space model. }
#'
#'   \item{KFS_transform}{How the non-diagonal \code{H} was handled. }
#'
#'   \item{logLik}{Value of the log-likelihood function. Only returned for fully
#'   Gaussian models. }
#'
#'   \item{a}{One-step-ahead predictions of states, \eqn{a_t = E(\alpha_t | y_{t-1},
#'   \ldots, y_{1})}{a[t] = E(\alpha[t] | y[t-1], \ldots, y[1])}. }
#'
#'   \item{P}{Non-diffuse parts of the error covariance matrix of predicted states,
#'   \eqn{P_t = Var(\alpha_t | y_{t-1}, \ldots, y_{1})
#'   }{P[t] = Var(\alpha[t] | y[t-1], \ldots, y[1])}. }
#'
#'   \item{Pinf}{Diffuse part of the error covariance matrix of predicted states.
#'   Only returned for Gaussian models. }
#'
#'   \item{t}{One-step-ahead predictions of signals, \eqn{E(Z_t\alpha_t | y_{t-1},
#'   \ldots, y_{1})}{E(Z[t]\alpha[t] | y[t-1], \ldots, y[1])}. }
#'
#'   \item{P_theta}{Non-diffuse part of \eqn{Var(Z_t\alpha_t | y_{t-1}, \ldots,
#'   y_{1})}{Var(Z[t]\alpha[t] | y[t-1], \ldots, y[1])}. }
#'
#'   \item{m}{One-step-ahead predictions \eqn{f(\theta_t) | y_{t-1}, \ldots,
#'   y_{1})}{f(\theta[t]) | y[t-1], \ldots, y[1])}, where \eqn{f} is the
#'   inverse link function. In case of Poisson distribution these predictions are
#'   multiplied with exposure \eqn{u_t}{u[t]}.  }
#'
#'   \item{P_mu}{Non-diffuse part of \eqn{Var(f(\theta_t) |
#'   y_{t-1}, \ldots, y_{1})}{Var(f(\theta[t]) | y[t-1], \ldots, y[1])}.
#'   In case of Poisson distribution this is \eqn{Var(u_t f(\theta_t) | y_{t-1},
#'   \ldots, y_{1})}{Var(t[t]f(\theta[t]) | y[t-1], \ldots, y[1])}.
#'   If \code{nsim = 0}, only diagonal elements (variances) are computed, using the
#'   Delta method. }
#'
#'   \item{alphahat}{Smoothed estimates of states, \eqn{E(\alpha_t | y_1, \ldots,
#'   y_n)}{E(\alpha[t] | y[1], \ldots, y[n])}. }
#'
#'   \item{V}{Error covariance matrices of smoothed states, \eqn{Var(\alpha_t | y_1,
#'   \ldots, y_n)}{Var(\alpha[t] | y[1], \ldots, y[n])}. }
#'
#'   \item{thetahat}{Smoothed estimates of signals, \eqn{E(Z_t\alpha_t | y_1,
#'   \ldots, y_n)}{E(Z[t]\alpha[t] | y[1], \ldots, y[n])}. }
#'
#'   \item{V_theta}{Error covariance matrices of smoothed signals
#'   \eqn{Var(Z[t]\alpha_t | y_1, \ldots, y_n).}{Var(Z[t]\alpha[t]
#'   | y[1], \ldot , y[n])}. }
#'
#'   \item{muhat}{Smoothed estimates of \eqn{f(\theta_t) | y_1, \ldots,
#'   y_n)}{f(\theta[t]) | y[1], \ldots, y[n])}, where \eqn{f} is the inverse
#'   link function, or in Poisson case \eqn{u_t f(\theta_t) | y_1, \ldots,
#'   y_n)}{u[t]f(\theta[t]) | y[1], \ldots, y[n])}, where \eqn{u} is the exposure term. }
#'
#'
#'   \item{V_mu}{Error covariances \eqn{Cov(f(\theta_t)| y_1, \ldots,
#'   y_n)}{Cov(f(\theta[t]) | y[1], \ldots, y[n])} (or the covariances of
#'   \eqn{u_t f(\theta_t)}{u[t]f(\theta[t])} given the data in case of Poisson
#'   distribution). If \code{nsim = 0}, only diagonal elements (variances) are
#'   computed, using the Delta method.  }
#'
#'   \item{etahat}{Smoothed disturbance terms \eqn{E(\eta_t | y_1, \ldots,
#'   y_n)}{E(\eta[t] | y[1], \ldots, y[n])}. Only for Gaussian models. }
#'
#'   \item{V_eta}{Error covariances \eqn{Var(\eta_t | y_1, \ldots, y_n)}{Var(\eta[t]
#'   | y[1], \ldots, y[n])}. }
#'
#'   \item{epshat}{Smoothed disturbance terms \eqn{E(\epsilon_{t,i} | y_1,
#'   \ldots, y_n)}{E(\epsilon[t,i] | y[1], \ldots, y[n])}. Note that due to
#'   the possible diagonalization these are on transformed scale.
#'   Only for Gaussian models. }
#'
#'   \item{V_eps}{Diagonal elements of \eqn{Var(\epsilon_{t} | y_1, \ldots,
#'   y_n)}{Var(\epsilon[t] | y[1], \ldots, y[n])}. Note that due to the
#'   diagonalization the off-diagonal elements are zero.
#'   Only for Gaussian models.  }
#'
#'   \item{iterations}{The number of iterations used in linearization of
#'   non-Gaussian model. }
#'
#'   \item{v}{Prediction errors \eqn{v_{t,i} = y_{t,i} - Z_{i,t}a_{t,i},
#'   i = 1, \ldots,p}{v[t,i] = y[t,i] - Z[i,t]a[t,i], i = 1, \ldots, p}, where
#'   \deqn{a_{t,i} = E(\alpha_t | y_{t,i-1}, \ldots, y_{t,1}, \ldots,
#'   y_{1,1})}{a[t,i] = E(\alpha[t] | y[t,i-1], \ldots, y[t,1], \ldots, y[1,1])}.
#'   Only returned for Gaussian models.  }
#'
#'   \item{F}{Prediction error variances \eqn{Var(v_{t,i})}{Var(v[t,i])}. Only
#'   returned for Gaussian models. }
#'
#'   \item{Finf}{Diffuse part of prediction error variances. Only returned for Gaussian
#'   models. }
#'
#'   \item{d}{The last time index of diffuse phase, i.e. the non-diffuse
#'   phase began at time \eqn{d+1}. Only returned for Gaussian models.  }
#'
#'
#'   \item{j}{The last observation index \eqn{i} of \eqn{y_{i,t}}{y[i,t]} of the
#'   diffuse phase. Only returned for Gaussian models.  }
#'
#'   In addition, if argument \code{simplify = FALSE}, list contains following
#'   components:
#'
#'   \item{K}{Covariances \eqn{Cov(\alpha_{t,i}, y_{t,i} | y_{t,i-1}, \ldots,
#'   y_{t,1}, y_{t-1}, \ldots , y_{1}), \quad i = 1, \ldots, p}{Cov(\alpha[t,i],
#'   y[t,i] | y[t,i-1], \ldots, y[t,1], y[t-1], \ldots, y[1]), i = 1, \ldots, p}.
#'   }
#'
#'   \item{Kinf}{Diffuse part of \eqn{K_t}{K[t]}.  }
#'
#'   \item{r}{Weighted sums of innovations \eqn{v_{t+1}, \ldots, v_{n}}{v[t+1],
#'   \ldots, v[n]}.  Notice that in literature \eqn{t} in \eqn{r_t}{r[t]} goes from
#'   \eqn{0, \ldots, n}. Here \eqn{t = 1, \ldots, n + 1}. Same applies to all \eqn{r} and
#'   \eqn{N} variables.  }
#'
#'   \item{r0, r1}{Diffuse phase decomposition of \eqn{r_t}{r[t]}.  }
#'
#'   \item{N}{Covariances \eqn{Var(r_t)}{Var(r[t])}.  }
#'
#'   \item{N0, N1, N2}{Diffuse phase decomposition of \eqn{N_t}{N[t]}.   }
#'
#' @references Koopman, S.J. and Durbin J. (2000).  Fast filtering and
#' smoothing for non-stationary time series models, Journal of American
#' Statistical Assosiation, 92, 1630-38.  \cr
#'
#' Koopman, S.J. and Durbin J. (2001).  Time Series Analysis by State Space
#' Methods. Oxford: Oxford University Press.  \cr
#'
#' Koopman, S.J. and Durbin J. (2003).  Filtering and smoothing of state vector
#' for diffuse state space models, Journal of Time Series Analysis, Vol. 24,
#' No. 1.  \cr
#'
KFS <-  function(model, filtering, smoothing, simplify = TRUE,
  transform = c("ldl", "augment"), nsim = 0, theta, maxiter = 50,
  convtol = 1e-08) {

  # Check that the model object is of proper form
  is.SSModel(model, na.check = TRUE, return.logical = FALSE)
  if (missing(filtering)) {
    if (all(model$distribution == "gaussian")) {
      filtering <- "state"
    } else filtering <- "none"
  } else {
    filtering <- match.arg(arg = filtering,
      choices = c("state", "signal", "mean", "none"),
      several.ok = TRUE)
    if ("signal" %in% filtering && all(model$distribution == "gaussian")) {
      filtering[filtering == "signal"] <- "mean"
      filtering <- unique(filtering)
    }
  }
  if (missing(smoothing)) {
    smoothing <- c("state", "mean")
  } else {
    smoothing <- match.arg(arg = smoothing,
      choices = c("state", "signal", "disturbance",
        "mean", "none"), several.ok = TRUE)

    if ("signal" %in% smoothing && all(model$distribution == "gaussian")) {
      smoothing[smoothing == "signal"] <- "mean"
      smoothing <- unique(smoothing)
    }
    if ("disturbance" %in% smoothing && !all(model$distribution == "gaussian")) {
      warning("disturbance smoothing is not supported for non-gaussian models.")
      smoothing <- smoothing[smoothing != "disturbance"]
    }
  }
  p <- attr(model, "p")
  m <- attr(model, "m")
  k <- attr(model, "k")
  n <- attr(model, "n")
  tv <- attr(model, "tv")
  ymiss <- is.na(model$y)
  storage.mode(ymiss) <- "integer"
  out <- list(model = model)

  # non-Gaussian case
  if (any(model$distribution != "gaussian")) {

    if(maxiter<1)
      stop("Argument maxiter must a positive integer. ")

    # initial values for theta
    if (missing(theta)) {
      theta <- initTheta(model$y, model$u, model$distribution)
    } else theta <- array(theta, dim = c(n, p))


    if (nsim > 0) {
      simtmp <- simHelper(model, nsim, TRUE)

      if (!("none" %in% filtering)) {
        filterout <- .Fortran(fngfilter, NAOK = TRUE, model$y, ymiss, tv,
          model$Z, model$T, model$R, model$Q, model$a1, model$P1, model$P1inf,
          model$u, theta,
          pmatch(x = model$distribution,
            table = c("gaussian", "poisson", "binomial", "gamma", "negative binomial"),
            duplicates.ok = TRUE),
          p, n, m, k, simtmp$nNonzeroP1inf,
          simtmp$nNonzeroP1, as.integer(nsim), simtmp$epsplus,
          simtmp$etaplus, simtmp$aplus1, simtmp$c2, model$tol, info = integer(1),
          maxiter = as.integer(maxiter),
          convtol = convtol, simtmp$zeroP1inf, length(simtmp$zeroP1inf),
          a = array(0, ("state" %in% filtering) * c(m - 1, n - 1) + 1),
          P = array(0, ("state" %in% filtering) * c(m - 1, m - 1, n - 1) + 1),
          theta = array(0, ("signal" %in% filtering) * c(p - 1, n - 1) + 1),
          P_theta = array(0, ("signal" %in% filtering) * c(p - 1, p - 1, n - 1) + 1),
          mu = array(0, ("mean" %in% filtering) * c(p - 1, n - 1) + 1),
          P_mu = array(0, ("mean" %in% filtering) * c(p - 1, p - 1, n - 1) + 1),
          as.integer("state" %in%  filtering), as.integer("signal" %in% filtering),
          as.integer("mean" %in%  filtering))
        if(filterout$info!=0){
          switch(as.character(filterout$info),
            "-3" = stop("Couldn't compute LDL decomposition of P1."),
            "-2" =  stop("Couldn't compute LDL decomposition of Q."),
            "1" =  stop("Gaussian approximation failed due to non-finite value in linear predictor."),
            "2" = stop("Gaussian approximation failed due to non-finite value of p(theta|y)."),
            "3" = warning("Maximum number of iterations reached, the approximation did not converge.")
          )
        }


        if ("state" %in% filtering) {
          out <- c(out, list(a = ts(t(filterout$a), start = start(model$y),
            frequency = frequency(model$y)), P = filterout$P))
          colnames(out$a) <- rownames(model$a1)
        }
        if ("signal" %in% filtering) {
          out <- c(out, list(t = ts(t(filterout$theta), start = start(model$y),
            frequency = frequency(model$y)), P_theta = filterout$P_theta))
          colnames(out$t) <- colnames(model$y)
        }
        if ("mean" %in% filtering) {
          out <- c(out, list(m = ts(t(filterout$mu), start = start(model$y),
            frequency = frequency(model$y)), P_mu = filterout$P_mu))
          colnames(out$m) <- colnames(model$y)
        }
        out <- c(out, iterations = filterout$maxiter)
      }
      if (!("none" %in% smoothing)) {
        smoothout <-
          .Fortran(fngsmooth, NAOK = TRUE, model$y, ymiss, tv,
            model$Z, model$T, model$R, model$Q, model$a1, model$P1, model$P1inf,
            model$u, theta,
            pmatch(x = model$distribution,
              table = c("gaussian", "poisson", "binomial", "gamma", "negative binomial"),
              duplicates.ok = TRUE),
            p, n, m, k, simtmp$nNonzeroP1inf, simtmp$nNonzeroP1,
            as.integer(nsim), simtmp$epsplus, simtmp$etaplus, simtmp$aplus1, simtmp$c2,
            model$tol, info = integer(1),
            maxiter = as.integer(maxiter),  convtol = convtol,
            simtmp$zeroP1inf, length(simtmp$zeroP1inf),
            alphahat = array(0, ("state" %in% smoothing) * c(m - 1, n - 1) + 1),
            V = array(0, ("state" %in% smoothing) * c(m - 1, m - 1, n - 1) + 1),
            thetahat = array(0, ("signal" %in% smoothing) * c(p - 1, n - 1) + 1),
            V_theta = array(0, ("signal" %in% smoothing) * c(p - 1, p - 1, n - 1) + 1),
            muhat = array(0, ("mean" %in% smoothing) * c(p - 1, n - 1) + 1),
            V_mu = array(0, ("mean" %in% smoothing) * c(p - 1, p - 1, n - 1) + 1),
            as.integer("state" %in%  smoothing), as.integer("signal" %in% smoothing),
            as.integer("mean" %in% smoothing))
        if(smoothout$info!=0){
          switch(as.character(smoothout$info),
            "-3" = stop("Couldn't compute LDL decomposition of P1."),
            "-2" =  stop("Couldn't compute LDL decomposition of Q."),
            "1" =  stop("Gaussian approximation failed due to non-finite value in linear predictor."),
            "2" = stop("Gaussian approximation failed due to non-finite value of p(theta|y)."),
            "3" = warning("Maximum number of iterations reached, the approximation did not converge.")
          )
        }
        if ("state" %in% smoothing) {
          out <- c(out, list(alphahat = ts(t(smoothout$alphahat), start = start(model$y),
            frequency = frequency(model$y)), V = smoothout$V))
          colnames(out$alphahat) <- rownames(model$a1)
        }
        if ("signal" %in% smoothing) {
          out <- c(out, list(thetahat = ts(t(smoothout$thetahat), start = start(model$y),
            frequency = frequency(model$y)), V_theta = smoothout$V_theta))
          colnames(out$thetahat) <- colnames(model$y)
        }
        if ("mean" %in% smoothing) {
          out <- c(out, list(muhat = ts(t(smoothout$muhat), start = start(model$y),
            frequency = frequency(model$y)), V_mu = smoothout$V_mu))
          colnames(out$muhat) <- colnames(model$y)
        }
        if ("none" %in% filtering)
          out <- c(out, iterations = smoothout$maxiter)
      }

      out$call <- match.call(expand.dots = FALSE)
      class(out) <- "KFS"
      return(out)
    } else {
      # Approximating model

      app <- .Fortran(fapprox, NAOK = TRUE, model$y, ymiss, tv,
        model$Z, model$T, model$R, Htilde = array(0, c(p, p, n)), model$Q,
        model$a1, model$P1, model$P1inf, p, n, m,
        k, theta = theta, model$u, ytilde = array(0, dim = c(n, p)),
        pmatch(x = model$distribution,
          table = c("gaussian", "poisson", "binomial", "gamma", "negative binomial"),
          duplicates.ok = TRUE),
        maxiter = as.integer(maxiter), model$tol, as.integer(sum(model$P1inf)),
        convtol, diff = double(1),lik=double(1), info=integer(1))

      if(app$info!=0){
        switch(as.character(app$info),
          "-3" = stop("Couldn't compute LDL decomposition of P1."),
          "-2" =  stop("Couldn't compute LDL decomposition of Q."),
          "1" =  stop("Gaussian approximation failed due to non-finite value in linear predictor."),
          "2" = stop("Gaussian approximation failed due to non-finite value of p(theta|y)."),
          "3" = warning("Maximum number of iterations reached, the approximation did not converge.")
        )
      }


      tsp(app$ytilde) <- tsp(model$y)
      model$y <- app$ytilde
      model$H <- app$Htilde

    }
  }

  if (all(model$distribution == "gaussian")) {
    transform <- match.arg(arg = transform, choices = c("ldl", "augment"))
    # Deal with the possible non-diagonality of H
    htol <- max(100, max(apply(model$H, 3, diag))) * .Machine$double.eps
    if (any(abs(apply(model$H, 3, "[", !diag(p))) > htol)) {
      model <- transformSSM(model, type = transform)
      KFS_transform <- transform
      tv <- attr(model, "tv")
      m <- attr(model, "m")
      k <- attr(model, "k")
    } else KFS_transform <- "none"
  } else KFS_transform <- "none"

  filtersignal <- ("signal" %in% filtering) || ("mean" %in% filtering)

  filterout <- .Fortran(fkfilter, NAOK = TRUE, model$y, ymiss, tv,
    model$Z, model$H, model$T, model$R, model$Q, model$a1, P1 = model$P1, model$P1inf,
    p, n, m, k, d = integer(1), j = integer(1), a = array(0, dim = c(m, n + 1)),
    P = array(0, dim = c(m, m, n + 1)), v = array(0, dim = c(p, n)),
    F = array(0, dim = c(p, n)), K = array(0, dim = c(m, p, n)),
    Pinf = array(0, dim = c(m, m, n + 1)), Finf = array(0, dim = c(p, n)),
    Kinf = array(0, dim = c(m, p, n)), lik = double(1), model$tol,
    as.integer(sum(model$P1inf)), theta = array(0, c(filtersignal * n, p)),
    P_theta = array(0, c(p, p, filtersignal * n)), as.integer(filtersignal))

  if (filterout$d == n & filterout$j == p)
    warning("Model is degenerate, diffuse phase did not end.")
  if (!all(is.finite(filterout$P))) {
    stop("Non-finite values on covariance matrix P. ")
  }
  filterout$Pinf <- filterout$Pinf[1:m, 1:m, 1:filterout$d, drop = FALSE]
  if (filterout$d > 0 & m > 1 & min(apply(filterout$Pinf, 3, diag)) < -.Machine$double.eps^0.5)
    warning(paste0("Possible error in diffuse filtering: Negative variances in Pinf, ",
      "check the model or try changing the tolerance parameter tol or P1/P1inf of the model."))
  if (sum(filterout$Finf > 0) != sum(diag(model$P1inf)))
    warning(paste0("Possible error in diffuse filtering: Number of nonzero elements in ",
      "Finf is not equal to the number of diffuse states. \n",
      "Either model is degenerate or numerical errors occured. ",
      "Check the model or try changing the tolerance parameter tol or P1/P1inf of the model."))

  if (filterout$d > 0) {
    filterout$Finf <- filterout$Finf[, 1:filterout$d, drop = FALSE]
    filterout$Kinf <- filterout$Kinf[, , 1:filterout$d, drop = FALSE]
  } else {
    filterout$Finf <- filterout$Kinf <- NA
  }
  out$KFS_transform <- KFS_transform
  out$d <- filterout$d
  out$j <- filterout$j

  if (all(model$distribution == "gaussian"))
    out$logLik <- filterout$lik

  if (!("none" %in% filtering)) {

    if (all(model$distribution == "gaussian")) {
      filterout$v[as.logical(t(ymiss))] <- NA
      if ("state" %in% filtering) {
        rownames(filterout$a) <- rownames(model$a1)
        out <- c(out, list(a = ts(t(filterout$a), start = start(model$y),
          frequency = frequency(model$y)), P = filterout$P, Pinf = filterout$Pinf))
      }
      if ("mean" %in% filtering) {
        colnames(filterout$theta) <- colnames(model$y)
        out <- c(out, list(m = ts(filterout$theta, start = start(model$y),
          frequency = frequency(model$y)), P_mu = filterout$P_theta))
      }
      out <- c(out, list(v = ts(t(filterout$v), start = start(model$y), frequency = frequency(model$y)),
        F = filterout$F, Finf = filterout$Finf))
      if (!simplify)
        out <- c(out, list(K = filterout$K, Kinf = filterout$Kinf))
    } else {
      if ("state" %in% filtering) {
        rownames(filterout$a) <- rownames(model$a1)
        out <- c(out, list(a = ts(t(filterout$a), start = start(model$y),
          frequency = frequency(model$y)), P = filterout$P))
      }
      if ("signal" %in% filtering) {
        colnames(filterout$theta) <- colnames(model$y)
        out <- c(out, list(t = ts(filterout$theta, start = start(model$y),
          frequency = frequency(model$y)), P_theta = filterout$P_theta))
      }
      if ("mean" %in% filtering) {
        out$m <- out$model$y
        out$P_mu <- array(0, c(p, p, n))
        for (i in 1:p) {
          out$m[, i] <- switch(model$distribution[i],
            gaussian = filterout$theta[, i],
            poisson = exp(filterout$theta[, i]) * model$u[, i],
            binomial = exp(filterout$theta[, i]) / (1 + exp(filterout$theta[, i])),
            gamma = exp(filterout$theta[, i]),
            `negative binomial` = exp(filterout$theta[, i]))
          out$P_mu[i, i, ] <- switch(model$distribution[i],
            gaussian = filterout$P_theta[i, i, ],
            poisson = filterout$P_theta[i, i, ] * out$m[,i]^2,
            binomial = filterout$P_theta[i, i, ] *
              (exp(filterout$theta[, i]) / (1 + exp(filterout$theta[, i]))^2)^2,
            gamma = filterout$P_theta[i, i, ] * out$m[,i]^2,
            `negative binomial` = filterout$P_theta[i, i, ] *
              out$m[,i]^2)
        }
      }
    }
  }


  if (!("none" %in% smoothing)) {
    transformedZ <- KFS_transform == "ldl" && ("signal" %in% smoothing || "mean" %in% smoothing)
    if (transformedZ){
      originalZ <- out$model$Z
    } else originalZ <- double(1)
    smoothout <- .Fortran(fgsmoothall, NAOK = TRUE, ymiss, tv, model$Z,
      model$H, model$T, model$R, model$Q, p, n, m,
      k, filterout$d, filterout$j, filterout$a, filterout$P, filterout$v,
      filterout$F, filterout$K, r = array(0, dim = c(m, n + 1)),
      r0 = array(0, dim = c(m, filterout$d + 1)),
      r1 = array(0, dim = c(m, filterout$d + 1)),
      N = array(0, dim = c(m, m, n + 1)),
      N0 = array(0, dim = c(m, m, filterout$d + 1)),
      N1 = array(0, dim = c(m, m, filterout$d + 1)),
      N2 = array(0, dim = c(m, m, filterout$d + 1)), filterout$Pinf, filterout$Kinf,
      filterout$Finf, alphahat = array(0, dim = c(m, n)*("state" %in% smoothing)),
      V = array(0, dim = c(m, m, n)),
      epshat = array(0, dim = c(p, n)*("disturbance" %in% smoothing)),
      V_eps = array(0, dim = c(p, n)*("disturbance" %in% smoothing)),
      etahat = array(0, dim = c(k, n)*("disturbance" %in% smoothing)),
      V_eta = array(0, dim = c(k, k, n)*("disturbance" %in% smoothing)),
      thetahat = array(0, dim = c(p, n)*("signal" %in% smoothing || "mean" %in% smoothing)),
      V_theta = array(0, dim = c(p, p, n)*("signal" %in% smoothing || "mean" %in% smoothing)),
      as.integer(transformedZ),
      originalZ, as.integer(dim(out$model$Z)[3] > 1),
      as.integer(KFS_transform != "augment"),
      as.integer("state" %in% smoothing), as.integer("disturbance" %in% smoothing),
      as.integer(("signal" %in% smoothing || "mean" %in% smoothing)))

    if (m > 1 & min(apply(smoothout$V, 3, diag)) < -.Machine$double.eps^0.5)
      warning(paste0("Possible error in smoothing: Negative variances in V, ",
        "check the model or try changing the tolerance parameter tol or P1/P1inf of the model."))
    if ("state" %in% smoothing) {
      out$alphahat <- ts(t(smoothout$alphahat), start = start(model$y), frequency = frequency(model$y))
      colnames(out$alphahat) <- rownames(model$a1)
      out$V <- smoothout$V
    }
    if ("disturbance" %in% smoothing) {
      out$etahat <- ts(t(smoothout$etahat), start = start(model$y), frequency = frequency(model$y))
      colnames(out$etahat) <- rownames(model$Q[, , 1])
      out$V_eta <- smoothout$V_eta
      if (KFS_transform != "augment") {
        out$epshat <- ts(t(smoothout$epshat), start = start(model$y), frequency = frequency(model$y))
        colnames(out$epshat) <- rownames(model$H[, , 1])
        out$V_eps <- smoothout$V_eps
      }
    }
    if ("signal" %in% smoothing) {
      out$thetahat <- ts(t(smoothout$thetahat), start = start(model$y), frequency = frequency(model$y))
      colnames(out$thetahat) <- rownames(model$H[, , 1])
      out$V_theta <- smoothout$V_theta
    }
    if ("mean" %in% smoothing) {
      out$muhat <- out$model$y
      out$V_mu <- array(0, c(p, p, n))
      for (i in 1:p) {
        out$muhat[, i] <- switch(model$distribution[i],
          gaussian = smoothout$thetahat[i, ],
          poisson = exp(smoothout$thetahat[i, ]) * model$u[, i],
          binomial = exp(smoothout$thetahat[i, ]) / (1 + exp(smoothout$thetahat[i, ])),
          gamma = exp(smoothout$thetahat[i, ]),
          `negative binomial` = exp(smoothout$thetahat[i, ]))
        out$V_mu[i, i, ] <- switch(model$distribution[i],
          gaussian = smoothout$V_theta[i, i, ],
          poisson = smoothout$V_theta[i, i, ] * out$muhat[, i]^2,
          binomial = smoothout$V_theta[i, i, ] *
            (exp(smoothout$thetahat[i, ]) / (1 + exp(smoothout$thetahat[i, ]))^2)^2,
          gamma = smoothout$V_theta[i, i, ] * out$muhat[, i]^2,
          `negative binomial` = smoothout$V_theta[i, i, ] * out$muhat[, i]^2)
      }
    }
    if (!simplify && all(model$distribution == "gaussian"))
      out <- c(out, list(r = smoothout$r, r0 = smoothout$r0, r1 = smoothout$r1,
        N = smoothout$N, N0 = smoothout$N0, N1 = smoothout$N1, N2 = smoothout$N2))
  }
  out$call <- match.call(expand.dots = FALSE)
  class(out) <- "KFS"
  out
}
