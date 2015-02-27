#'Maximum Likelihood Estimation of a State Space Model
#'
#'Function \code{fitSSM} finds the maximum likelihood estimates for unknown 
#'parameters of an arbitary state space model, given the user-defined model 
#'updating function.
#'
#'This function is simple wrapper around optim. For maximum performance in 
#'complicated problems, user is suggested to use problem specific codes with 
#'calls to \code{logLik} method directly.
#'
#'After updating the model, if \code{checkfn} is non-missing and 
#'\code{checkfn(model)==TRUE}, \code{-logLik(model,check.model=FALSE)} is 
#'computed, otherwise objective function returns \code{.Machine$double.xmax}. If
#'\code{checkfn} is missing, checking is done via \code{is.SSModel} by 
#'\code{-logLik(model,check.model=TRUE)}
#'
#'Note that for non-Gaussian models derivative-free optimization methods such as
#'Nelder-Mead might be more reliable than methods which use finite difference
#'approximations. This is due to noise caused by the relative stopping criterion
#'used for finding approximating Gaussian model. In most cases this does not
#'seem to cause any problems though.
#'
#'@export
#'@param inits Initial values for \code{optim}
#'@param model Model object of class \code{SSModel}.
#'@param updatefn User defined function which updates the model given the 
#'  parameters. Must be of form \code{updatefn(pars, model,...)}, i.e. must 
#'  contain ellipsis \code{...}. If not supplied, a default function is used, 
#'  which estimates the values marked as NA in time invariant covariance 
#'  matrices Q and H. Note that the default function cannot be used with 
#'  trigonometric seasonal components as its covariance structure is of form 
#'  \eqn{\sigma}{\sigma}I.
#'@param checkfn Optional function for checking the validity of the model. See 
#'  details.
#'@param ... Further arguments for functions \code{optim}, \code{updatefn} and 
#'  \code{logLik.SSModel}, such as \code{method='BFGS'}.
#'@return A list with elements \item{optim.out}{Output from function 
#'  \code{optim}. } \item{model}{Model with estimated parameters. }
#' @examples
#' 
#' # Example function for updating covariance matrices H and Q 
#' # (also used as a default function in fitSSM)
#' 
#' updatefn <- function(pars,model,...){  
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
fitSSM <- 
  function(model, inits, updatefn, checkfn,...) {
   
    if (missing(updatefn)) {
      # use default updating function, only for time invariant covariance matrices
      estH <- !identical(model$H, "Omitted") && any(is.na(model$H))
      estQ <- any(is.na(model$Q))
      if ((dim(model$H)[3]>1 && estH || (dim(model$Q)[3] > 1) && estQ)) 
        stop("No model updating function supplied, but cannot use default 
             function as the covariance matrices are time varying.")
      updatefn <- function(pars, model, ...) {
        if(estQ){
          Q <- as.matrix(model$Q[, , 1])
          naQd <- which(is.na(diag(Q)))
          naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
          Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
          diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
          Q[naQd, naQd][naQnd] <- pars[length(naQd) + 1:length(naQnd)]
          model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
        } else naQnd<-naQd<-NULL
        if (estH){
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
    is.SSModel(updatefn(inits, model, ...), na.check = TRUE, 
               return.logical = FALSE)
    
    # use is.SSModel as default checking function
    if (missing(checkfn)) {
      likfn <- function(pars, model, ...) {
        model <- updatefn(pars, model, ...)
         -logLik(object = model, check.model = TRUE, ...)       
      }
    } else {
        likfn <- function(pars, model, ...) {
          model <- updatefn(pars, model, ...)
          if (checkfn(model)) {
            return(-logLik(object = model, check.model = FALSE, ...))
          } else return(.Machine$double.xmax^0.75)
        }
    }
    out <- NULL
    out$optim.out <- optim(par = inits, fn = likfn, model = model, ...)
    out$model <- updatefn(out$opt$par, model, ...)
    out
  } 
