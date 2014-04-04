# Function for computing the prediction and confidence intervals for non-Gaussian models
# Used by predict.SSModel method

interval <- function(model, interval = c("confidence", "prediction"), level, type = c("response", "link"), 
                     states=NULL,nsim, se.fit = TRUE, timespan, prob = TRUE,maxiter=50) {
  
  interval <- match.arg(interval)
  type <- match.arg(type)
  if (type == "link" && interval == "prediction") 
    stop("Prediction intervals can only be computed at response scale.")
  
  m <- attr(model, "m")
  n <- as.integer(length(timespan))
  p <- attr(model, "p")
  n2 <- as.integer(length(timespan))
  # Generate sample via importance sampling
  imp <- importanceSSM(model, ifelse(identical(states, as.integer(1:m)), "signal", "states"), 
                       nsim = nsim, antithetics = TRUE,maxiter=maxiter)
  nsim <- as.integer(4 * nsim)
  w <- imp$weights/sum(imp$weights)
  if (!identical(states, 1:attr(model, "m"))) # use only selected states
    imp$samples <- .Fortran(fzalpha, as.integer(dim(model$Z)[3] > 1), 
                            model$Z[, , if (dim(model$Z)[3] > 1) timespan else 1, drop = FALSE], 
                            imp$samples[timespan, , drop = FALSE], 
                            signal = array(0, c(n2, p, nsim)), as.integer(p), n2, as.integer(m), 
                            nsim, as.integer(length(states)), states)$signal
  
  for(j in 1:p)
    if(model$distribution[j]=="poisson")
      imp$samples[timespan, j, ]<-imp$samples[timespan, j, ]+log(model$u[timespan,j])
  
  # compute intervals using weighted sample
  if (type == "link") {
    int <- lapply(1:p, function(j) sapply(timespan, function(i) {
      or <- order(imp$samples[i, j, ])
      c(imp$samples[i, j, or][which.max(cumsum(w[or]) >= (1 - level)/2)], imp$samples[i, j, or][which.max(cumsum(w[or]) >= 
                                                                                                            1 - (1 - level)/2)])
    }))
  } else {
    for (j in 1:p) {
      imp$samples[timespan, j, ] <- switch(model$distribution[j], gaussian = imp$samples[timespan, j, ], poisson = exp(imp$samples[timespan, j, ]), binomial = (if (!prob) model$u[timespan, j] else 1) * exp(imp$samples[timespan, 
                                                                                                                                                                                                                          j, ])/(1 + exp(imp$samples[timespan, j, ])), gamma = exp(imp$samples[timespan, j, ]), `negative binomial` = exp(imp$samples[timespan, 
                                                                                                                                                                                                                                                                                                                                                      j, ]))
    }
    
    if (interval == "confidence") {
      int <- lapply(1:p, function(j) {
        sapply(timespan, function(i) {
          or <- order(imp$samples[i, j, ])
          c(imp$samples[i, j, or][which.max(cumsum(w[or]) >= (1 - level)/2)], imp$samples[i, j, or][which.max(cumsum(w[or]) >= 
                                                                                                                1 - (1 - level)/2)])
        })
      })
    } else {
      # sample from observational density
      int <- lapply(1:p, function(j) {
        sapply(timespan, function(i) {
          sample_mu <- sample(imp$samples[i, j, ], size = nsim, replace = TRUE, prob = w)
          q<-quantile(switch(model$distribution[j], gaussian = rnorm(n = nsim, mean = sample_mu, sd = model$u[i, j]), poisson = rpois(n = nsim, 
                                                                                                                                      lambda = sample_mu), binomial = rbinom(n = nsim, size = (if (!prob) model$u[i, j] else 1), prob = sample_mu/(if (!prob) model$u[i, 
                                                                                                                                                                                                                                                                      j] else 1)), gamma = rgamma(n = nsim, shape = model$u[i, j], scale = sample_mu/model$u[i, j]), `negative binomial` = rnbinom(n = nsim, 
                                                                                                                                                                                                                                                                                                                                                                                                   size = model$u[i, j], mu = sample_mu)), 
                      prob = c((1 - level)/2, 1 - (1 - level)/2),
                      type=switch(model$distribution[j], gaussian=,gamma=7, poisson=,binomial=,`negative binomial`=1))                  
        })
      })
      
    }
    
  }
  varmean <- .Fortran(fvarmeanw, imp$samples[timespan, , ], w, as.integer(p), n2, nsim, mean = array(0, c(n2, p)), var = array(0, c(n2, p)), 
                      as.integer(se.fit))
  if (se.fit) {
    pred <- lapply(1:p, function(j) cbind(fit = cbind(varmean$mean[, j], lwr = int[[j]][1, ], upr = int[[j]][2, ]), se.fit = sqrt(varmean$var[, 
                                                                                                                                              j])))
  } else {
    pred <- lapply(1:p, function(j) cbind(varmean$mean[, j], lwr = int[[j]][1, ], upr = int[[j]][2, ]))
  }
  pred
} 
