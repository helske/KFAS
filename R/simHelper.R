
# simulate error terms from standard normal distribution
# transformation to correct distribution is done in Fortran

simHelper <- function(model, ymiss, nsim, antithetics) {
  
  epsplus <- array(0, c(attr(model, "p"), attr(model, "n"), nsim))
  etaplus <- array(0, c(attr(model, "k"), attr(model, "n"), nsim))
  aplus1 <- array(0, dim = c(attr(model, "m"), nsim))
  if(any(model$distribution != "gaussian")){
    x <- array(t(!ymiss), c(attr(model, "p"), attr(model, "n"), nsim))
  } else {
    x <- array(abs(apply(model$H, 3, diag)) > model$tol, c(attr(model, "p"), attr(model, "n"))) & (!t(ymiss))
    x <- array(x, c(attr(model, "p"), attr(model, "n"), nsim))  
  }  
  
  dfeps <- sum(x)/nsim
  x2 <- array(abs(apply(model$Q, 3, diag)) > model$tol, c(attr(model, "k"), (attr(model, "n") - 1) * attr(model, "tv")[5] + 1))
  x2 <- array(x2, c(attr(model, "k"), attr(model, "n"), nsim))
  dfeta <- sum(x2)/nsim
  nonzeroP1 <- which(diag(model$P1) > model$tol)
  nNonzeroP1 <- length(nonzeroP1)   
  dfu <- dfeps + dfeta + nNonzeroP1
  u <- rnorm(dfu * nsim, mean = 0, sd = 1)
  if (dfeps > 0) 
    epsplus[x] <- u[1:(dfeps * nsim)]
  if (dfeta > 0) 
    etaplus[x2] <- u[(dfeps * nsim + 1):(dfeps * nsim + dfeta * nsim)]
  if (nNonzeroP1 > 0) 
    aplus1[nonzeroP1, ] <- u[(dfeps * nsim + dfeta * nsim + 1):(dfu * nsim)]
  c2 <- numeric(nsim)
  # for second antithetic
  if (antithetics) {
    for (i in 1:nsim) {
      u <- c(etaplus[, , i], epsplus[, , i], aplus1[, i])
      c2[i] <- t(u) %*% c(u)
    }
    q <- pchisq(c2, df = dfu)
    c2 <- sqrt(qchisq(1 - q, dfu)/c2)
  }
  list(epsplus = epsplus, etaplus = etaplus, aplus1 = aplus1, c2 = c2, 
       nonzeroP1 = as.integer(nonzeroP1), nNonzeroP1 = nNonzeroP1,
       zeroP1inf = which(diag(model$P1inf) > 0), nNonzeroP1inf = as.integer(sum(model$P1inf)))
}