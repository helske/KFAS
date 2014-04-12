#' Linear Gaussian Approximation for Exponential Family State Space Model
#'
#' Function \code{approxSMM} computes the linear Gaussian approximation of a
#' state space model where observations follow an exponential family distribution.
#'
#' The linear Gaussian approximating model is defined by
#' \deqn{\tilde y_t = Z_t \alpha_t + \epsilon_t, \quad \epsilon_t \sim N(0,\tilde H_t),}{ytilde[t] = Z[t]\alpha[t] + \epsilon[t], \epsilon[t] ~ N(0,Htilde[t]),}
#' \deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t, \quad \eta_t \sim N(0,Q_t),}{\alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t], \eta[t] ~ N(0,Q[t]),}
#' and \eqn{\alpha_1 \sim N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])}, 
#' where \eqn{\tilde y}{ytilde} and \eqn{\tilde H}{Htilde} are chosen in a way that the linear
#' Gaussian approximating model has the same conditional mode of \eqn{\theta=Z\alpha} 
#' given the observations \eqn{y} as the original non-gaussian model. 
#' Models also have a same curvature at the mode.
#'
#' The approximation of the exponential family state space model is based on iterative weighted 
#' least squares method, see McCullagh and Nelder (1983) p.31 and Durbin Koopman (2012) p. 243.
#'  
#' @seealso Importance sampling of non-Gaussian state space models \code{\link{importanceSSM}}, 
#' construct a \code{SSModel} object \code{\link{SSModel}}, and examples in \code{\link{KFAS}}.
#' @export
#' @param model A non-Gaussian state space model object of class \code{SSModel}.
#' @param theta Initial values for conditional mode theta.
#' @param maxiter The maximum number of iterations used in approximation Default is 50.
#' @param tol Tolerance parameter for convergence checks.
#'  Iterations are continued until 
#'  \eqn{tol>abs(dev_{old}-dev_{new})/(abs(dev_{new})+0.1))}.
#' @return An object which contains the approximating Gaussian state space model with following additional components:
#' \item{thetahat}{Mode of \eqn{p(\theta|y)}. }
#' \item{iterations}{Number of iterations used. }
approxSSM <- function(model, theta, maxiter = 50, tol = 1e-15) {
  
  # Check that the model object is of proper form
  is.SSModel(model, na.check = TRUE, return.logical = FALSE)
  if (all(model$distribution == "gaussian")) 
    stop("Model is completely Gaussian, nothing to approximate.")
  p <- attr(model, "p")
  m <- attr(model, "m")
  k <- attr(model, "k")
  n <- attr(model, "n")
  
  tv <- array(0, dim = 5)
  tv[1] <- dim(model$Z)[3] > 1
  tv[2] <- 1
  tv[3] <- dim(model$T)[3] > 1
  tv[4] <- dim(model$R)[3] > 1
  tv[5] <- dim(model$Q)[3] > 1
  
  ymiss <- is.na(model$y)
  storage.mode(ymiss) <- "integer"
  if(is.null(maxiter)) maxiter<-50
  
  # initial values for linear predictor theta
  if (missing(theta) || is.null(theta)) {
    theta <- sapply(1:p, function(i) 
      switch(model$distribution[i], 
             gaussian = model$y[, i], 
             poisson = log(pmax(model$y[, i]/model$u[, i], 0.1, na.rm = TRUE)), 
             binomial = qlogis((ifelse(is.na(model$y[, i]), 0.5, model$y[, i]) + 0.5)/(model$u[, i] + 1)), 
             gamma = log(pmax(model$y[, i], 1, na.rm = TRUE)), 
             `negative binomial` = log(pmax(model$y[, i], 1/6, na.rm = TRUE))))
  } else theta <- array(theta, dim = c(n, p))
 
  # call Fortran subroutine for model approximation
  out <- .Fortran(fapprox, NAOK = TRUE, model$y, ymiss, as.integer(tv), model$Z, 
                  model$T, model$R, Htilde = array(0, c(p, p, n)), model$Q, 
                  model$a1, model$P1, model$P1inf, as.integer(p), as.integer(n),
                  as.integer(m),as.integer(k), theta = theta, model$u, 
                  ytilde = array(0, dim = c(n, p)), 
                  pmatch(x = model$distribution, 
                         table = c("gaussian", "poisson", "binomial", "gamma", "negative binomial"), 
                         duplicates.ok = TRUE), maxiter = as.integer(maxiter), 
                  model$tol, as.integer(sum(model$P1inf)), as.double(tol),diff=double(1))
  
  if (!is.finite(out$diff)){
    stop("Non-finite difference in approximation algoritm.")
  }
  if(out$maxiter==maxiter){
    warning(paste("Maximum number of iterations reached, 
                  the approximation algorithm did not converge. Latest difference was",out$diff))
  }
  
  model$distribution <- rep("gaussian", p)
  model$y[] <- out$ytilde
  model$y[as.logical(ymiss)] <- NA
  model$H <- out$Htilde
  model$thetahat <- out$theta
  model$iterations <- out$maxiter
  model$difference <- out$diff
  class(model) <- c("approxSSM", "SSModel")
  invisible(model)
} 
