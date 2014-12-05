#' Linear Gaussian Approximation for Exponential Family State Space Model
#'
#' Function \code{approxSMM} performs the linear Gaussian approximation of a
#' state space model where observations follow an exponential 
#' family distribution.
#' 
#' This function is rarely needed itself, it is mainly available for 
#' illustrative and debugging purposes. The underlying Fortran code is used by 
#' other functions of KFAS for non-Gaussian modelling.
#' 
#' The linear Gaussian approximating model is defined by
#' \deqn{\tilde y_t = Z_t \alpha_t + \epsilon_t, 
#' \quad \epsilon_t \sim N(0,\tilde H_t),}{
#' ytilde[t] = Z[t]\alpha[t] + \epsilon[t], \epsilon[t] ~ N(0,Htilde[t]),}
#' \deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t, 
#' \quad \eta_t \sim N(0,Q_t),}{
#' \alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t], \eta[t] ~ N(0,Q[t]),}
#' 
#' and \eqn{\alpha_1 \sim N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])},
#'  
#' where \eqn{\tilde y}{ytilde} and \eqn{\tilde H}{Htilde} are chosen in a way 
#' that the linear Gaussian approximating model has the same conditional mode 
#' of \eqn{\theta=Z\alpha} given the observations \eqn{y} as the original 
#' non-gaussian model. Models also have a same curvature at the mode.
#'
#' The approximation of the exponential family state space model is based on 
#' iterative weighted least squares method, see McCullagh and Nelder (1983) 
#' p.31 and Durbin Koopman (2012) p. 243.
#' #'  
#' @seealso 
#' \code{\link{importanceSSM}},
#' \code{\link{SSModel}},
#' \code{\link{KFAS}}.
#' @export
#' @param model A non-Gaussian state space model object of class \code{SSModel}.
#' @param theta Initial values for conditional mode theta.
#' @param maxiter The maximum number of iterations used in approximation.
#'  Default is 50.
#' @param tol Tolerance parameter for convergence checks.
#'  Iterations are continued until the scaled norm between three successive iterations is smaller than \code{tol}.
#' @param stepmax Maximum stepsize used in Gaussian approximation. Only used for non-Gaussian models.
#' @return An object which contains the approximating Gaussian state space model
#'  with following additional components:
#' \item{thetahat}{Mode of \eqn{p(\theta|y)}. }
#' \item{iterations}{Number of iterations used. }
#' @references
#' Koopman, S.J. and Durbin J. (2012).  Time Series Analysis by State Space
#' Methods. Second edition. Oxford: Oxford University Press.
approxSSM <- 
  function(model, theta, maxiter = 50, tol = 1e-08,stepmax) {
    if (is.null(maxiter)) 
      maxiter <- 50
    if (is.null(tol)) 
      tol <- 1e-08
    if(maxiter<1)
      stop("Argument maxiter must a positive integer. ")    
    # Check that the model object is of proper form
    is.SSModel(model, na.check = TRUE, return.logical = FALSE)
    if (all(model$distribution == "gaussian")) 
      stop("Model is completely Gaussian, nothing to approximate.")
    
    
    p <- attr(model, "p")
    m <- attr(model, "m")
    k <- attr(model, "k")
    n <- attr(model, "n")
    tv <- attr(model, "tv")
    ymiss <- is.na(model$y)
    storage.mode(ymiss) <- "integer"    
    
    # Initial values for linear predictor theta
    if (missing(theta) || is.null(theta)) {
      theta <- init_theta(model$y, model$u, model$distribution)
    } else theta <- array(theta, dim = c(n, p))
    
    dist <- pmatch(x = model$distribution, duplicates.ok = TRUE,
                   table = c("gaussian", "poisson", 
                             "binomial", "gamma", "negative binomial"))
    
    if(missing(stepmax))
      stepmax<-max(1,sqrt(sum(theta^2)))
    
    # Call Fortran subroutine for model approximation
    out <- 
      .Fortran(fapprox, NAOK = TRUE, model$y, ymiss, tv, model$Z, 
               model$T, model$R, Htilde = array(0, c(p, p, n)), model$Q, 
               model$a1, model$P1, model$P1inf, p, n, 
               m, k, theta = theta, model$u, 
               ytilde = array(0, dim = c(n, p)), dist, 
               maxiter = as.integer(maxiter), model$tol, 
               as.integer(sum(model$P1inf)), 
               as.double(tol), diff = double(1),double(1),as.double(stepmax),info=integer(1))
    if(out$info!=0){
      if (out$info==1) {
        stop("Non-finite value of likelihood or linear predictor in approximation algorithm.")
      }
      if(out$info==2){
        warning(paste("Maximum number of iterations reached, lLatest difference was", signif(out$diff,3)))
      }   
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
