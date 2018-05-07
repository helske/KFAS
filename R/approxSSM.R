#' Linear Gaussian Approximation for Exponential Family State Space Model
#'
#' Function \code{approxSMM} performs a linear Gaussian approximation of an
#' exponential family state space model.
#'
#' This function is rarely needed itself, it is mainly available for
#' illustrative and debugging purposes. The underlying Fortran code is used by
#' other functions of KFAS for non-Gaussian state space modelling.
#'
#' The linear Gaussian approximating model is defined by
#' \deqn{\tilde y_t = Z_t \alpha_t + \epsilon_t, \quad \epsilon_t \sim N(0,\tilde H_t),}{
#' ytilde[t] = Z[t]\alpha[t] + \epsilon[t], \epsilon[t] ~ N(0,Htilde[t]),}
#' \deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t, \quad \eta_t \sim N(0,Q_t),}{
#' \alpha[t+1] = T[t]\alpha[t] + R[t]\eta[t], \eta[t] ~ N(0,Q[t]),}
#'
#' and \eqn{\alpha_1 \sim N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])},
#'
#' where \eqn{\tilde y}{ytilde} and \eqn{\tilde H}{Htilde} are chosen in a way
#' that the linear Gaussian approximating model has the same conditional mode of
#' \eqn{\theta=Z\alpha} given the observations \eqn{y} as the original
#' non-Gaussian model. Models also have a same curvature at the mode.
#'
#' The approximation of the exponential family state space model is based on
#' iterative weighted least squares method, see McCullagh and Nelder (1983) p.31
#' and Durbin Koopman (2012) p. 243.
#'
#' @seealso \code{\link{importanceSSM}}, \code{\link{SSModel}},
#'   \code{\link{KFS}}, \code{\link{KFAS}}.
#' @export
#' @param model A non-Gaussian state space model object of class \code{SSModel}.
#' @param theta Initial values for conditional mode theta.
#' @param maxiter The maximum number of iterations used in approximation.
#'   Default is 50.
#' @param tol Tolerance parameter for convergence checks.
#' @return An object of class \code{SSModel} which contains the approximating Gaussian state space model
#'   with following additional components:
#'   \item{thetahat}{Mode of \eqn{p(\theta|y)}. }
#'   \item{iterations}{Number of iterations used. }
#'   \item{difference}{Relative difference in the last step of approximation
#'   algorithm. }
#' @references \itemize{ \item McCullagh, P. and Nelder, J. A. (1983).
#'   Generalized linear models. Chapman and Hall.
#'   \item Koopman, S.J. and Durbin, J. (2012).
#'   Time Series Analysis by State Space Methods. Second edition. Oxford University Press. }
#' @examples
#'
#' # A Gamma example modified from ?glm (with log-link)
#' clotting <- data.frame(
#'   u = c(5,10,15,20,30,40,60,80,100),
#'   lot1 = c(118,58,42,35,27,25,21,19,18),
#'   lot2 = c(69,35,26,21,18,16,13,12,12))
#'
#' glmfit1 <- glm(lot1 ~ log(u), data = clotting, family = Gamma(link = "log"))
#' glmfit2 <- glm(lot2 ~ log(u), data = clotting, family = Gamma(link = "log"))
#'
#' # Model lot1 and lot2 together (they are still assumed independent)
#' # Note that Gamma distribution is parameterized by 1/dispersion i.e. shape parameter
#' model <- SSModel(cbind(lot1, lot2) ~ log(u),
#'                 u = 1/c(summary(glmfit1)$dispersion, summary(glmfit2)$dispersion),
#'                 data = clotting, distribution = "gamma")
#' approxmodel <- approxSSM(model)
#'
#' # Conditional modes of linear predictor:
#' approxmodel$thetahat
#' cbind(glmfit1$linear.predictor, glmfit2$linear.predictor)
#'
#' KFS(approxmodel)
#' summary(glmfit1)
#' summary(glmfit2)
#'
#' # approxSSM uses modified step-halving for more robust convergence than glm:
#' y <- rep (0:1, c(15, 10))
#' suppressWarnings(glm(formula = y ~ 1, family = binomial(link = "logit"), start = 2))
#' model <- SSModel(y~1, dist = "binomial")
#' KFS(model, theta = 2)
#' KFS(model, theta = 7)
approxSSM <- function(model, theta, maxiter = 50, tol = 1e-08) {

  if (maxiter < 1) {
    stop("Argument maxiter must a positive integer. ")
  }
  if (tol < 0) {
    stop("Argument tol must be non-negative. ")
  }
  # Check that the model object is of proper form
  is.SSModel(model, na.check = TRUE, return.logical = FALSE)
  if (all(model$distribution == "gaussian")) {
    stop("Model is completely Gaussian, nothing to approximate. ")
  }

  p <- attr(model, "p")
  m <- attr(model, "m")
  k <- attr(model, "k")
  n <- attr(model, "n")
  tv <- attr(model, "tv")
  ymiss <- is.na(model$y)
  storage.mode(ymiss) <- "integer"

  # Initial values for linear predictor theta
  if (missing(theta) || is.null(theta)) {
    theta <- initTheta(model$y, model$u, model$distribution)
  } else theta <- array(theta, dim = c(n, p))

  dist <- pmatch(x = model$distribution, duplicates.ok = TRUE,
    table = c("gaussian", "poisson",
      "binomial", "gamma", "negative binomial"))

  # Call Fortran subroutine for model approximation
  out <-
    .Fortran(fapprox, NAOK = TRUE, model$y, ymiss, tv, model$Z, model$T,
      model$R, Htilde = array(0, c(p, p, n)), model$Q, model$a1,
      model$P1, model$P1inf, p, n, m, k, theta = theta, model$u,
      ytilde = array(0, dim = c(n, p)), dist,
      maxiter = as.integer(maxiter), model$tol,
      as.integer(sum(model$P1inf)), as.double(tol), diff = double(1),
      double(1),info = integer(1))
  if (out$info != 0) {
    warning(switch(as.character(out$info),
      "-3" = "Couldn't compute LDL decomposition of P1.",
      "-2" =  "Couldn't compute LDL decomposition of Q.",
      "1" = paste0("Gaussian approximation failed due to ",
        "non-finite value in linear predictor."),
      "2" = paste0("Gaussian approximation failed due to ",
        "non-finite value of p(theta|y)."),
      "3" = paste0("Maximum number of iterations reached, latest ",
        "difference was ", signif(out$diff, 3))
    ))
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
