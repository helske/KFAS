#' Maximum Likelihood Estimation of a State Space Model
#'
#' Function \code{fitSSM} finds the maximum likelihood estimates for unknown
#' parameters of an arbitary state space model, given the user-defined model
#' updating function.
#'
#' This function is simple wrapper around \code{optim}. For optimal performance in
#' complicated problems, it is more efficient to use problem specific codes with
#' calls to \code{logLik} method directly.
#'
#' After updating the model, if  \code{checkfn(model) == TRUE}, the
#' log-likelihood is computed using a call \code{-logLik(model,check.model = FALSE)}.
#' Otherwise objective function returns \code{.Machine$double.xmax}.
#' If \code{checkfn} is missing, default \code{checkfn} is used which checks
#' that system matrices contain only finite values (components \code{y} and
#' \code{Z} are not tested), and checks that \code{Q} and \code{H} contain only
#' values smaller than 1e7 (larger values can cause numerical instability).'
#'
#' Note that for non-Gaussian models derivative-free optimization methods such as
#' Nelder-Mead might be more reliable than methods which use finite difference
#' approximations. This is due to noise caused by the relative stopping criterion
#' used for finding approximating Gaussian model. In most cases this does not
#' seem to cause any problems though.
#'
#' @export
#' @param inits Initial values for \code{optim}
#' @param model Model object of class \code{SSModel}.
#' @param updatefn User defined function which updates the model given the
#'   parameters. Must be of form \code{updatefn(pars, model, ...)}.
#'   If not supplied, a default function is used,
#'   which estimates the values marked as NA in unconstrained time-invariant
#'   covariance matrices Q and H. Note that the default \code{updatefn} function
#'   cannot be used with trigonometric seasonal components as its covariance
#'   structure is of form \eqn{\sigma}{\sigma}I.
#' @param checkfn Optional function for checking the validity of the model. See
#'   details.
#' @param update_args Optional list containing additional arguments to \code{updatefn}.
#' @param ... Further arguments for functions \code{optim} and
#'  \code{logLik.SSModel}, such as \code{nsim = 1000} and \code{method = "BFGS"}.
#'@return A list with elements
#'\item{optim.out}{Output from function \code{optim}. }
#'\item{model}{Model with estimated parameters. }
#' @examples
#'
#' # Example function for updating covariance matrices H and Q
#' # (also used as a default function in fitSSM)
#'
#' updatefn <- function(pars, model){
#'   if(any(is.na(model$Q))){
#'     Q <- as.matrix(model$Q[,,1])
#'     naQd  <- which(is.na(diag(Q)))
#'     naQnd <- which(upper.tri(Q[naQd,naQd]) & is.na(Q[naQd,naQd]))
#'     Q[naQd,naQd][lower.tri(Q[naQd,naQd])] <- 0
#'     diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
#'     Q[naQd,naQd][naQnd] <- pars[length(naQd)+1:length(naQnd)]
#'     model$Q[naQd,naQd,1] <- crossprod(Q[naQd,naQd])
#'   }
#'  if(!identical(model$H,'Omitted') && any(is.na(model$H))){#'
#'    H<-as.matrix(model$H[,,1])
#'    naHd  <- which(is.na(diag(H)))
#'    naHnd <- which(upper.tri(H[naHd,naHd]) & is.na(H[naHd,naHd]))
#'    H[naHd,naHd][lower.tri(H[naHd,naHd])] <- 0
#'    diag(H)[naHd] <-
#'      exp(0.5 * pars[length(naQd)+length(naQnd)+1:length(naHd)])
#'    H[naHd,naHd][naHnd] <-
#'      pars[length(naQd)+length(naQnd)+length(naHd)+1:length(naHnd)]
#'    model$H[naHd,naHd,1] <- crossprod(H[naHd,naHd])
#'    }
#'  model
#'}
#'
#' # Example function for checking the validity of covariance matrices.
#'
#' checkfn <- function(model){
#'   #test positive semidefiniteness of H and Q
#'   inherits(try(ldl(model$H[,,1]),TRUE),'try-error') ||
#'   inherits(try(ldl(model$Q[,,1]),TRUE),'try-error')
#' }
#'
#'
#' model <- SSModel(Nile ~ SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
#'
#' #function for updating the model
#' update_model <- function(pars, model) {
#'   model["H"] <- pars[1]
#'   model["Q"] <- pars[2]
#'   model
#' }
#'
#' #check that variances are non-negative
#' check_model <- function(model) {
#'   (model["H"] > 0 && model["Q"] > 0)
#' }
#'
#' fit <- fitSSM(inits = rep(var(Nile)/5, 2), model = model,
#'                  updatefn = update_model, checkfn = check_model)

fitSSM <- function(model, inits, updatefn, checkfn, update_args = NULL, ...) {

  is_gaussian <- all(model$distribution == "gaussian")

  if (missing(updatefn)) {
    # use default updating function, only for time invariant covariance matrices
    estH <- is_gaussian && any(is.na(model$H))
    estQ <- any(is.na(model$Q))
    if ((dim(model$H)[3] > 1 && estH || (dim(model$Q)[3] > 1) && estQ))
      stop("No model updating function supplied, but cannot use default
             function as the covariance matrices are time varying.")
    updatefn <- function(pars, model) {
      if (estQ) {
        Q <- as.matrix(model$Q[, , 1])
        naQd <- which(is.na(diag(Q)))
        naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
        Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
        diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
        Q[naQd, naQd][naQnd] <- pars[length(naQd) + 1:length(naQnd)]
        model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
      } else naQnd <- naQd <- NULL
      if (estH) {
        H <- as.matrix(model$H[, , 1])
        naHd <- which(is.na(diag(H)))
        naHnd <- which(upper.tri(H[naHd, naHd]) & is.na(H[naHd, naHd]))
        H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
        diag(H)[naHd] <- exp(0.5 * pars[length(naQd) + length(naQnd) +
            1:length(naHd)])
        H[naHd, naHd][naHnd] <- pars[length(naQd) + length(naQnd) +
            length(naHd) + 1:length(naHnd)]
        model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
      }
      model
    }
  }
  # Check that the model object is of proper form
  is.SSModel(do.call(updatefn, args = c(list(pars = inits, model = model), update_args)),
    na.check = TRUE, return.logical = FALSE)

  # initial values for theta can be computed beforehand
  if (!is_gaussian && is.null(list(...)$theta)) {
    theta <- initTheta(model$y, model$u, model$distribution)

  } else theta <- NULL

  if (missing(checkfn)) {
    # check for nonfinite values and overly large variances/covariances
    if (is_gaussian) {
      checkfn <- function(model){
        all(sapply(c("H", "T", "R", "Q", "a1", "P1", "P1inf"),
          function(x) {
            all(is.finite(model[[x]]))
          })) && max(model$Q) <= 1e7 && max(model$H) <= 1e7
      }
    } else {
      checkfn <- function(model){
        all(sapply(c("u", "T", "R", "Q", "a1", "P1", "P1inf"),
          function(x) {
            all(is.finite(model[[x]]))
          })) &&
          max(model$Q) <= 1e7
      }
    }

  }
  likfn <- function(pars, model, ...) {
    model <- do.call(updatefn, args = c(list(pars = pars, model = model), update_args))
    if (checkfn(model)) {
      return(-logLik(object = model, check.model = FALSE, theta = theta, ...))
    } else return(.Machine$double.xmax ^ 0.75)
  }

  out <- NULL
  out$optim.out <- optim(par = inits, fn = likfn, model = model, ...)
  out$model <- do.call(updatefn, args = c(list(pars = out$optim.out$par, model = model), update_args))
  # check that the obtained model is of proper form
  is.SSModel(out$model, na.check = TRUE, return.logical = FALSE)
  out
}
