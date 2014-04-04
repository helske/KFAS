#' Maximum Likelihood Estimation of a State Space Model
#'
#' Function \code{fitSSM} finds the maximum likelihood estimates for
#' unknown parameters of an arbitary state space model, given the user-defined model updating function.
#'
#'
#' @export
#' @param inits Initial values for \code{optim}
#' @param model Model object of class \code{SSModel}.
#' @param updatefn User defined function which updates the model given the parameters. 
#' Must be of form \code{updatefn(pars, model,...)}, i.e. must contain ellipsis \code{...}. 
#' If not supplied, a default function is used, which estimates the values marked as NA in time invariant covariance matrices Q and H.
#' @param checkfn Optional function for model checking. 
#' If supplied, after updating the model, if \code{checkfn(model)} returns TRUE, -log-likelihood is computed, 
#' otherwise \code{.Machine$double.xmax} is returned. See examples.
#' If not supplied, check.model=TRUE is used for checking possible NA or Inf values, see ?logLik.SSModel.
#' @param ... Further arguments for functions \code{optim}, \code{updatefn} and \code{logLik.SSModel}, such as \code{method='BFGS'}.
#' @return A list with elements 
#' \item{optim.out}{Output from function \code{optim}. }
#' \item{model}{Model with estimated parameters. }
#' @examples
#' 
#' # Example function for updating covariance matrices H and Q 
#' # (also used as a default function in fitSSM)
#' 
#' updatefn <- function(pars,model,...){       
#' Q<-as.matrix(model$Q[,,1])         
#' naQd  <- which(is.na(diag(Q))) 
#' naQnd <- which(upper.tri(Q[naQd,naQd]) & is.na(Q[naQd,naQd]))  
#' Q[naQd,naQd][lower.tri(Q[naQd,naQd])] <- 0
#' diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
#' Q[naQd,naQd][naQnd] <- pars[length(naQd)+1:length(naQnd)]
#' model$Q[naQd,naQd,1] <- crossprod(Q[naQd,naQd])
#' if(!identical(model$H,'Omitted')){
#'    H<-as.matrix(model$H[,,1])  
#'    naHd  <- which(is.na(diag(H)))
#'    naHnd <- which(upper.tri(H[naHd,naHd]) & is.na(H[naHd,naHd]))
#'    H[naHd,naHd][lower.tri(H[naHd,naHd])] <- 0         
#'    diag(H)[naHd] <- exp(0.5 * pars[length(naQd)+length(naQnd)+1:length(naHd)])
#'    H[naHd,naHd][naHnd] <- pars[length(naQd)+length(naQnd)+length(naHd)+1:length(naHnd)]
#'    model$H[naHd,naHd,1] <- crossprod(H[naHd,naHd])
#'  }            
#'  
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


fitSSM <- function(model, inits, updatefn, checkfn, ...) {
    
    # use default updating function, only for time invariant covariance matrices
    if (missing(updatefn)) {
        if (max(dim(model$H)[3], dim(model$Q)[3]) > 1) 
            stop("No model updating function supplied, but cannot use default function as covariance matrices are time varying.")
        updatefn <- function(pars, model, ...) {
            Q <- as.matrix(model$Q[, , 1])
            naQd <- which(is.na(diag(Q)))
            naQnd <- which(upper.tri(Q[naQd, naQd]) & is.na(Q[naQd, naQd]))
            Q[naQd, naQd][lower.tri(Q[naQd, naQd])] <- 0
            diag(Q)[naQd] <- exp(0.5 * pars[1:length(naQd)])
            Q[naQd, naQd][naQnd] <- pars[length(naQd) + 1:length(naQnd)]
            model$Q[naQd, naQd, 1] <- crossprod(Q[naQd, naQd])
            if (!identical(model$H, "Omitted")) {
                H <- as.matrix(model$H[, , 1])
                naHd <- which(is.na(diag(H)))
                naHnd <- which(upper.tri(H[naHd, naHd]) & is.na(H[naHd, naHd]))
                H[naHd, naHd][lower.tri(H[naHd, naHd])] <- 0
                diag(H)[naHd] <- exp(0.5 * pars[length(naQd) + length(naQnd) + 1:length(naHd)])
                H[naHd, naHd][naHnd] <- pars[length(naQd) + length(naQnd) + length(naHd) + 1:length(naHnd)]
                model$H[naHd, naHd, 1] <- crossprod(H[naHd, naHd])
            }
            model
        }
        
    }
    
    is.SSModel(updatefn(inits, model, ...), na.check = TRUE, return.logical = FALSE)
    
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
                return(-logLik(object = model, ...))
            } else return(.Machine$double.xmax)
        }
        
    }
    
    
    out <- NULL
    out$optim.out <- optim(par = inits, fn = likfn, model = model, ...)
    out$model <- updatefn(out$opt$par, model, ...)
    out
} 
