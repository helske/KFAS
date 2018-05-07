#' State Space Model Predictions
#'
#' Function \code{predict.SSModel} predicts the future observations of a state
#' space model of class \code{\link{SSModel}}.
#'
#' For non-Gaussian models, the results depend whether importance sampling is
#' used (\code{nsim>0}). without simulations, the confidence intervals are based
#' on the Gaussian approximation of \eqn{p(\alpha | y)}. Confidence intervals in
#' response scale are computed in linear predictor scale, and then transformed
#' to response scale. The prediction intervals are not supported. With
#' importance sampling, the confidence intervals are computed as the empirical
#' quantiles from the weighted sample, whereas the prediction intervals contain
#' additional step of simulating the response variables from the sampling
#' distribution \eqn{p(y|\theta^i)}.
#'
#' Predictions take account the uncertainty in state estimation
#' (given the prior distribution for the initial states), but not the uncertainty
#' of estimating the parameters in the system matrices (i.e. \eqn{Z}, \eqn{Q} etc.).
#' Thus the obtained confidence/prediction intervals can underestimate the true
#' uncertainty for short time series and/or complex models.
#'
#' If no simulations are used, the standard errors in response scale are
#' computed using the Delta method.
#'
#' @export
#' @importFrom stats predict end qnorm
#' @aliases predict predict.SSModel
#' @param object Object of class \code{SSModel}.
#' @param newdata A compatible \code{SSModel} object to be added in the end of
#'   the old object for which the predictions are required. If omitted,
#'   predictions are either for the past data points, or if argument
#'   \code{n.ahead} is given, \code{n.ahead} time steps ahead.
#' @param n.ahead Number of steps ahead at which to predict. Only used if
#'   \code{newdata} is omitted. Note that when using \code{n.ahead}, object
#'   cannot contain time varying system matrices.
#' @param interval Type of interval calculation.
#' @param level Confidence level for intervals.
#' @param type Scale of the prediction, \code{"response"} or \code{"link"}.
#' @param states Which states are used in computing the predictions. Either a
#'   numeric vector containing the indices of the corresponding states, or a
#'   character vector defining the types of the corresponding states. Possible choices are
#'   \code{"all"},  \code{"level"}, \code{"slope"},
#'   \code{"trend"},  \code{"regression"}, \code{"arima"}, \code{"custom"},
#'   \code{"cycle"} or \code{"seasonal"}, where \code{"trend"} extracts all states relating to trend.
#'    These can be combined. Default is \code{"all"}.
#' @param nsim Number of independent samples used in importance sampling. Used
#'   only for non-Gaussian models.
#' @param se.fit If TRUE, standard errors of fitted values are computed. Default is FALSE.
#' @param prob if TRUE (default), the predictions in binomial case are
#'   probabilities instead of counts.
#' @param maxiter The maximum number of iterations used in approximation Default
#'   is 50. Only used for non-Gaussian model.
#' @param filtered If \code{TRUE}, compute predictions based on filtered
#' (one-step-ahead) estimates. Default is FALSE i.e. predictions are based on
#' all available observations given by user. For diffuse phase,
#' interval bounds and standard errors of fitted values are set to \code{-Inf}/\code{Inf}
#' (If the interest is in the first time points it might be useful to use
#' non-exact diffuse initialization.).
#' @param \dots Ignored.
#' @return A matrix or list of matrices containing the predictions, and
#'   optionally standard errors.
#' @examples
#'
#' set.seed(1)
#' x <- runif(n=100,min=1,max=3)
#' y <- rpois(n=100,lambda=exp(x-1))
#' model <- SSModel(y~x,distribution="poisson")
#' xnew <- seq(0.5,3.5,by=0.1)
#' newdata <- SSModel(rep(NA,length(xnew))~xnew,distribution="poisson")
#' pred <- predict(model,newdata=newdata,interval="prediction",level=0.9,nsim=100)
#' plot(x=x,y=y,pch=19,ylim=c(0,25),xlim=c(0.5,3.5))
#' matlines(x=xnew,y=pred,col=c(2,2,2),lty=c(1,2,2),type="l")
#'
#' model <- SSModel(Nile~SSMtrend(1,Q=1469),H=15099)
#' pred <- predict(model,n.ahead=10,interval="prediction",level=0.9)
#' pred
predict.SSModel <- function(object, newdata, n.ahead,
  interval = c("none", "confidence", "prediction"), level = 0.95,
  type = c("response", "link"), states = NULL, se.fit = FALSE,  nsim = 0,
  prob = TRUE, maxiter = 50, filtered = FALSE, ...) {

  interval <- match.arg(interval)
  type <- match.arg(type)
  # Check that the model object is of proper form
  is.SSModel(object, na.check = TRUE, return.logical = FALSE)
  m <- attr(object, "m")
  p <- attr(object, "p")
  k <- attr(object, "k")
  if (missing(states)) {
    states <- as.integer(1:m)
  } else {
    if (is.numeric(states)) {
      states <- as.integer(states)
      if (min(states) < 1 | max(states) > m)
        stop("Vector states should contain the indices or types of the states which are combined.")
    } else {
      states <- match.arg(arg = states, choices = c("all", "arima", "custom", "level","slope",
        "cycle", "seasonal", "trend", "regression"),
        several.ok = TRUE)
      if ("all" %in% states) {
        states <- as.integer(1:attr(object, "m"))
      } else {
        if("trend" %in% states)
          states <- c(states, "level", "slope")
        states <- which(attr(object, "state_types") %in% states)
      }
    }
  }
  gaussianmodel <- all(object$distribution == "gaussian")
  if (!missing(newdata) && !is.null(newdata)) {
    #expand model with newdata object
    # Check that the model object is of proper form
    is.SSModel(newdata, na.check = TRUE, return.logical = FALSE)
    if (p != attr(newdata, "p"))
      stop("Different number of time series for 'object' and 'newdata'.")
    if (m != attr(newdata, "m"))
      stop("Different number of states for 'object' and 'newdata'.")
    if (k != attr(newdata, "k"))
      stop("Different number of disturbance terms for 'object' and 'newdata'.")
    if (!identical(object$distribution, newdata$distribution))
      stop("Different distributions for 'object' and 'newdata'")
    no <- attr(object, "n")
    nn <- attr(newdata, "n")
    n <- attr(object, "n") <- no + nn
    timespan <- (no + 1):n
    object$y <- ts(rbind(object$y, newdata$y),
      start = start(object$y), frequency = frequency(object$y))
    endtime <- end(object$y)
    tvo <- attr(object, "tv")
    tvn <- attr(newdata, "tv")
    same <- function(x, y) isTRUE(all.equal(x, y, tolerance = 0,
      check.attributes = FALSE))
    if (tvo[1] || tvn[1] || !same(object$Z, newdata$Z)) {
      object$Z <- array(data = c(array(object$Z, dim = c(m, p, no)),
        array(newdata$Z, dim = c(m, p, nn))), dim = c(p, m, n))
      attr(object, "tv")[1] <- 1L
    }
    if (gaussianmodel && (tvo[2] || tvn[2] || !same(object$H, newdata$H))) {
      object$H <- array(data = c(array(object$H, dim = c(p, p, no)),
        array(newdata$H, dim = c(p, p, nn))), dim = c(p, p, n))
      attr(object, "tv")[2] <- 1L
    } else if(!gaussianmodel) object$u <- rbind(object$u, matrix(newdata$u, nn, p))

    if (tvo[3] || tvn[3] || !same(object$T, newdata$T)) {
      object$T <- array(data = c(array(object$T, dim = c(m, m, no)),
        array(newdata$T, dim = c(m, m, nn))), dim = c(m, m, n))
      attr(object, "tv")[3] <- 1L
    }
    if (tvo[4] || tvn[4] || !same(object$R, newdata$R)) {
      object$R <- array(data = c(array(object$R, dim = c(m, k, no)),
        array(newdata$R, dim = c(m, k, nn))), dim = c(m, k, n))
      attr(object, "tv")[4] <- 1L
    }
    if (tvo[5] || tvn[5] || !same(object$Q, newdata$Q)) {
      object$Q <- array(data = c(array(data = object$Q, dim = c(k, k, no)),
        array(data = newdata$Q, dim = c(k, k, nn))), dim = c(k, k, n))
      attr(object, "tv")[5] <- 1L
    }
  } else {
    if (!missing(n.ahead) && !is.null(n.ahead)) {
      tv <- attr(object, "tv")
      if(ifelse(gaussianmodel,any(tv), any(c(apply(object$u, 2, function(x) length(unique(x)) > 1))) || any(tv[-2])))
        stop("Model contains time varying system matrices, cannot use argument 'n.ahead'. Use 'newdata' instead.")
      timespan <- attr(object, "n") + 1:n.ahead
      n <- attr(object, "n") <- attr(object, "n") + as.integer(n.ahead)
      endtime<-end(object$y) + c(0, n.ahead)
      object$y <- window(object$y, end = endtime, extend = TRUE)
      if (any(object$distribution != "gaussian"))
        object$u <- rbind(object$u, matrix(object$u[1, ], nrow = n.ahead,
          ncol = ncol(object$u), byrow = TRUE))
    } else {
      timespan <- 1:attr(object, "n")
      endtime <- end(object$y)
    }
  }
  if (!gaussianmodel && interval == "prediction") {
    if (type == "link")
      stop("Prediction intervals can only be computed at response scale.")
    if (nsim < 1)
      stop("Cannot compute prediction intervals for non-gaussian models without importance sampling.")
  }
  pred <- vector("list", length = p)
  if (gaussianmodel) { #Gaussian case
    if (identical(states, as.integer(1:m))) {
      if (filtered) {
        out <- KFS(model = object, filtering = "mean", smoothing = "none")
        names(out)[6:7] <- c("muhat", "V_mu")
        if (out$d > 0) {
        out$V_mu[,,1:out$d] <- Inf #diffuse phase
        }
      } else {
        out <- KFS(model = object, filtering = "none", smoothing = "mean")
      }
    } else {
      if (filtered) {
        out <- KFS(model = object, filtering = "state", smoothing = "none")
        d <- out$d
        out <- signal(out,  states = states, filtered = TRUE)
        names(out) <- c("muhat", "V_mu")
        if (d > 0) {
        out$V_mu[,,1:d] <- Inf #diffuse phase
        }
      } else {
        out <- signal(KFS(model = object, filtering = "none", smoothing = "state"),
          states = states)
        names(out) <- c("muhat", "V_mu")
      }

    }
    for (i in 1:p) {
      pred[[i]] <- cbind(fit = out$muhat[timespan, i],
        switch(interval, none = NULL,
          confidence = out$muhat[timespan, i] +
            qnorm((1 - level)/2) * sqrt(out$V_mu[i, i, timespan]) %o% c(1, -1),
          prediction = out$muhat[timespan, i] + qnorm((1 - level)/2) *
            sqrt(out$V_mu[i, i, timespan] +
                object$H[i, i, if (dim(object$H)[3] > 1) timespan else 1]) %o% c(1, -1)),
        se.fit = if (se.fit)
          sqrt(out$V_mu[i, i, timespan]))
      if (interval != "none")
        colnames(pred[[i]])[2:3] <- c("lwr", "upr")
    }
  } else {
    if (nsim < 1) { #Using approximating model
      if (identical(states, as.integer(1:m))) {
        if (filtered) {
          out <- KFS(model = object, filtering = "signal", smoothing = "none",
            maxiter = maxiter)
          names(out)[5:6] <- c("thetahat", "V_theta")
          if (out$d > 0) {
          out$V_theta[,,1:out$d] <- Inf #diffuse phase
          }
        } else {
          out <- KFS(model = object, smoothing = "signal", maxiter = maxiter)
        }
      } else {
        if (filtered) {
          out <- KFS(model = object, filtering = "state", smoothing = "none",
            maxiter = maxiter)
          d <- out$d
          out <- signal(out, states = states, filtered = TRUE)
          names(out) <- c("thetahat", "V_theta")
          if (d > 0) {
          out$V_theta[,,1:d] <- Inf #diffuse phase
          }
        } else {
          out <- signal(KFS(model = object, smoothing = "state", maxiter = maxiter),
            states = states)
          names(out) <- c("thetahat", "V_theta")
        }

      }

      for (i in 1:p) {
        pred[[i]] <- cbind(fit = out$thetahat[timespan, i] +
            (if (object$distribution[i] == "poisson")  log(object$u[timespan, i]) else 0),
          switch(interval, none = NULL,
            out$thetahat[timespan, i] +
              (if (object$distribution[i] == "poisson") log(object$u[timespan, i]) else 0) +
              qnorm((1 - level)/2) * sqrt(out$V_theta[i, i, timespan]) %o%
              c(1, -1)), se.fit = if (se.fit)
                sqrt(out$V_theta[i, i, timespan]))
        if (interval == "confidence")
          colnames(pred[[i]])[2:3] <- c("lwr", "upr")
      }
      if (type == "response") {
        if (se.fit) {
          tmp <- which(colnames(pred[[1]]) == "se.fit")
          for (i in 1:p) {
            pred[[i]][, "se.fit"] <- switch(object$distribution[i],
              gaussian = pred[[i]][,"se.fit"],
              poisson = pred[[i]][, "se.fit"] * exp(pred[[i]][, 1]),
              binomial = pred[[i]][, "se.fit"] *
                (if (!prob) object$u[timespan, i] else 1) *
                exp(pred[[i]][, 1])/(1 + exp(pred[[i]][, 1]))^2,
              gamma = pred[[i]][, "se.fit"] * exp(pred[[i]][, 1]),
              `negative binomial` = pred[[i]][, "se.fit"] * exp(pred[[i]][, 1]))
            pred[[i]][, -tmp] <- switch(object$distribution[i],
              gaussian = pred[[i]][, -tmp],
              poisson = exp(pred[[i]][, -tmp]),
              binomial = (if (!prob) object$u[timespan, i] else 1) *
                exp(pred[[i]][, -tmp])/(1 + exp(pred[[i]][, -tmp])),
              gamma = exp(pred[[i]][, -tmp]),
              `negative binomial` = exp(pred[[i]][, -tmp]))
          }
        } else {
          for (i in 1:p) pred[[i]] <- switch(object$distribution[i],
            gaussian = pred[[i]],
            poisson = exp(pred[[i]]),
            binomial = (if (!prob) object$u[timespan, i] else 1) *
              exp(pred[[i]])/(1 + exp(pred[[i]])),
            gamma = exp(pred[[i]]),
            `negative binomial` = exp(pred[[i]]))
        }
      }

    } else {# with importance sampling
      if (filtered) {
        d <- KFS(approxSSM(object, maxiter = maxiter), smoothing = "none")$d
      }

      if (interval == "none") {
        imp <- importanceSSM(object,
          ifelse(identical(states, as.integer(1:m)), "signal", "states"),
          nsim = nsim, antithetics = TRUE, maxiter = maxiter, filtered = filtered)
        nsim <- as.integer(4 * nsim)
        if (!identical(states, as.integer(1:m))) {
          imp$samples <- .Fortran(fzalpha, as.integer(dim(object$Z)[3] > 1),
            object$Z, imp$samples, signal = array(0, c(n, p, nsim)),
            p, m, n, nsim, length(states), states)$signal
        }

        w <- imp$weights/sum(imp$weights)
        if (type == "response") {
          for (i in 1:p) {
            imp$samples[timespan, i, ] <- switch(object$distribution[i],
              gaussian = imp$samples[timespan, i, ],
              poisson = object$u[timespan, i] * exp(imp$samples[timespan, i, ]),
              binomial = (if (!prob) object$u[timespan, i] else 1) *
                exp(imp$samples[timespan, i, ])/(1 + exp(imp$samples[timespan, i, ])),
              gamma = exp(imp$samples[timespan, i, ]),
              `negative binomial` = exp(imp$samples[timespan, i, ]))
          }
        } else {
          for (i in 1:p) if (object$distribution[i] == "poisson")
            imp$samples[timespan, i, ] <- imp$samples[timespan, i, ] + log(object$u[timespan,
              i])
        }
        varmean <- .Fortran(fvarmeanw, imp$samples[timespan, , , drop = FALSE], w,
          p, length(timespan),
          nsim, mean = array(0, c(length(timespan), p)),
          var = array(0, c(length(timespan), p)), as.integer(se.fit))

        if (se.fit) {
          if (filtered && d > 0) {
            varmean$var[1:d, ] <- Inf #diffuse phase
          }
          pred <- lapply(1:p, function(j) cbind(fit = varmean$mean[, j],
            se.fit = sqrt(varmean$var[, j])))

        } else {
          pred <- lapply(1:p, function(j) varmean$mean[, j])
        }
      } else {
        pred <- interval(object, interval = interval, level = level, type = type,
          states = states, nsim = nsim, se.fit = se.fit, timespan = timespan,
          prob = prob, maxiter = maxiter, filtered = filtered)
        if (filtered && d > 0) {
          for (i in 1:p) {
            pred[[i]][1:d, "lwr"] <- -Inf
            pred[[i]][1:d, "upr"] <- Inf #diffuse phase
            if (se.fit) {
              pred[[i]][1:d, "se.fit"] <- Inf #diffuse phase
            }
          }

        }
      }
    }
  }
  names(pred) <- colnames(object$y)
  pred <- lapply(pred, ts,end = endtime, frequency = frequency(object$y))
  if (p == 1) pred <- pred[[1]]
  pred
}
