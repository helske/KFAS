# Function for computing the prediction and confidence intervals for non-Gaussian
# models.
# Used by predict.SSModel method.
#' @importFrom stats quantile rnorm rbinom rpois rnbinom rgamma
interval <- function(model, interval = c("confidence", "prediction"), level,
  type = c("response",  "link"), states = NULL, nsim, se.fit = TRUE,
  timespan, prob = TRUE, maxiter = 50, filtered = FALSE) {

  interval <- match.arg(interval)

  type <- match.arg(type)
  if (type == "link" && interval == "prediction")
    stop("Prediction intervals can only be computed at response scale.")
  m <- attr(model, "m")
  n <- as.integer(length(timespan))
  p <- attr(model, "p")
  # Generate sample via importance sampling
  imp <- importanceSSM(model, ifelse(identical(states, as.integer(1:m)), "signal", "states"),
    nsim = nsim, antithetics = TRUE, maxiter = maxiter, filtered = filtered)
  nsim <- as.integer(4 * nsim)
  w <- imp$weights/sum(imp$weights)
  imp$samples <- imp$samples[timespan, , , drop = FALSE]
  # use only selected states
  if (!identical(states, 1:attr(model, "m")))
    imp$samples <- .Fortran(fzalpha, as.integer(dim(model$Z)[3] > 1),
      model$Z[, ,if (dim(model$Z)[3] > 1) timespan else 1, drop = FALSE],
      imp$samples, signal = array(0, c(n, p, nsim)),
      p, m, n, nsim, length(states), states)$signal

  for (j in 1:p) {
    if (model$distribution[j] == "poisson")
      imp$samples[, j, ] <- imp$samples[, j, ] + log(model$u[timespan, j])
  }
  # compute intervals using weighted sample
  if (type == "link") {
    int <- lapply(1:p, function(j) sapply(1:n, function(i) {
      or <- order(imp$samples[i, j, ])
      c(imp$samples[i, j, or][which.max(cumsum(w[or]) >= (1 - level)/2)],
        imp$samples[i, j, or][which.max(cumsum(w[or]) >= 1 - (1 - level)/2)])
    }))
  } else {
    for (j in 1:p) {
      imp$samples[, j, ] <-
        switch(model$distribution[j],
          gaussian = imp$samples[, j, ],
          poisson = exp(imp$samples[, j, ]),
          binomial = (if (!prob) c(model$u[timespan, j]) else 1) *
            exp(imp$samples[, j, ])/(1 + exp(imp$samples[, j, ])),
          gamma = exp(imp$samples[, j, ]),
          `negative binomial` = exp(imp$samples[, j, ]))
    }
    if (interval == "confidence") {
      int <- lapply(1:p, function(j) {
        sapply(1:n, function(i) {
          or <- order(imp$samples[i, j, ])
          c(imp$samples[i, j, or][which.max(cumsum(w[or]) >= (1 - level)/2)],
            imp$samples[i, j, or][which.max(cumsum(w[or]) >= 1 - (1 - level)/2)])
        })
      })
    } else {
      # sample from observational density
      int <- lapply(1:p, function(j) {
        sapply(1:n, function(i) {
          sample_mu <- sample(imp$samples[i, j, ], size = nsim, replace = TRUE,
            prob = w)
          quantile(switch(model$distribution[j],
            gaussian = rnorm(n = nsim, mean = sample_mu, sd = model$u[timespan[i], j]),
            poisson = rpois(n = nsim, lambda = sample_mu),
            binomial = rbinom(n = nsim, size = (if (!prob) model$u[timespan[i], j] else 1),
              prob = sample_mu/(if (!prob) model$u[timespan[i], j] else 1)),
            gamma = rgamma(n = nsim, shape = model$u[timespan[i], j],
              scale = sample_mu/model$u[timespan[i], j]),
            `negative binomial` = rnbinom(n = nsim, size = model$u[timespan[i],j], mu = sample_mu)),
            prob = c((1 - level)/2, 1 - (1 - level)/2))
        })
      })
    }
  }
  varmean <- .Fortran(fvarmeanw, imp$samples, w, p, n,
    nsim, mean = array(0, c(n, p)), var = array(0, c(n, p)), as.integer(se.fit))
  if (se.fit) {
    pred <- lapply(1:p, function(j) cbind(fit = cbind(varmean$mean[, j],
      lwr = int[[j]][1, ], upr = int[[j]][2, ]),
      se.fit = sqrt(varmean$var[, j])))
  } else {
    pred <- lapply(1:p, function(j) cbind(fit = varmean$mean[, j],
      lwr = int[[j]][1, ], upr = int[[j]][2, ]))
  }
  pred
}
