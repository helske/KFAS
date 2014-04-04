#' Extract Standardized Residuals from KFS output
#' @S3method rstandard KFS
#' @method rstandard KFS
#' @details For object of class KFS, several types of standardized residuals can be computed:
#' \itemize{
#' 
#' \item 'recursive': One-step ahead prediction residuals defined as 
#' \deqn{v_{t,i})/\sqrt{F_{i,t}},}
#' with residuals being undefined in diffuse phase. Only supported for fully Gaussian models.
#' 
#' 
#' \item 'pearson':  Standardized Pearson residuals 
#' \deqn{(y_{t,i}-\theta_{t,i})/\sqrt{V(\mu)_{t,i}\phi_i\sqrt{1-h_{t,i}}}, \quad i=1,\ldots,p,t=1,\ldots,n,}{(y[t,i]-\theta[t,i])/(V(\mu)[t,i]\phi[i](1-h[t,i]))^0.5, i=1,\ldots,p, t=1,\ldots,n,}
#'                   where \eqn{V(\mu_{t,i})}{V(\mu[t,i])} is the variance function of the model, 
#'                   \eqn{\phi_i}{\phi[i]} is the dispersion parameter and \eqn{h_{t,i}}{h[t,i]} 
#'                   is the hat value.  For gaussian models, these coincide with the smoothed 
#'                   \eqn{\epsilon} disturbance residuals.
#'
#' \item 'state':  Residuals based on the smoothed disturbance terms \eqn{\eta} are defined as
#' \deqn{L^{-1}_t \hat \eta_t, \quad t=1,\ldots,n,}{L^{-1}[t] \eta[t], t=1,\ldots,n,} where 
#' \eqn{L_t}{L[t]} is the lower triangular matrix from Cholesky decomposition of \eqn{V_{\eta,t}}{V[\eta,t]}.
#' 
#' \item 'deviance': Deviance residuals.
#' }
#'
#' @param model KFS object
#' @param type Type of residuals. See details.
#' @param ... Ignored.

rstandard.KFS <- function(model, type = c("recursive","deviance", "pearson", "state"), ...) {
  
  type <- match.arg(type)
  
  if ((type == "recursive" || type == "state") && any(model$model$distribution != "gaussian")) 
    stop("Recursive and state residuals are only supported for fully gaussian models.")
  
  recursive <- function(object) {
    if (is.null(object[["v",exact=TRUE]])) 
      stop("KFS object needs to contain innovations from state filtering. ")
    series <- object$v/sqrt(t(object$F))
    series[1:(object$d - 1), ] <- NA
    series[object$d, 1:object$j] <- NA
    series
  }
  
  
    
  pearson <- function(object) {
    if (all(object$model$distribution == "gaussian")) {
      w <- matrix(apply(object$model$H,3,diag), attr(object$model, "n"), attr(object$model, "p"),byrow=TRUE)    
    } else {
      w <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
      for (i in 1:attr(object$model, "p")) 
        w[, i] <- switch(object$model$distribution[i], 
                         gaussian = object$model$u[,i], 
                         poisson = 1, 
                         binomial = 1, 
                         gamma = 1/object$model$u[,i], 
                         `negative binomial` =1)
    }
    
    series <- object$model$y
    if(sum(bins<-object$model$distribution=="binomial")>0)
      series[, bins]<-
      series[, bins]/
      object$model$u[, bins]
    
    ((series - fitted(object))/sqrt(variance(object)))/sqrt(w*(1 - hatvalues(object)))
  }
  
  deviance <- function(object) {
    if (all(object$model$distribution == "gaussian")) { 
      w <- matrix(apply(object$model$H,3,diag), attr(object$model, "n"), attr(object$model, "p"),byrow=TRUE)
    } else {
      w <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
      for (i in 1:attr(object$model, "p")) 
        w[, i] <- switch(object$model$distribution[i], 
                         gaussian = object$model$u[,i], poisson = 1, 
                         binomial = 1, 
                         gamma = 1/object$model$u[,i], 
                         `negative binomial` =1)    
    }
    series <- object$model$y
    if(sum(bins<-object$model$distribution=="binomial")>0)
      series[, bins]<-
      series[, bins]/
      object$model$u[, bins]
    
    for (i in 1:attr(object$model, "p")) 
      series[, i] <- ifelse(series[, i]>object$muhat[,i],1,-1)*
      sqrt(switch(object$model$distribution[i], 
                  gaussian = (series[,i] - object$muhat[, i])^2, 
                  poisson =
                    2*(series[,i]*log(ifelse(series[,i] == 0, 1, series[,i]/object$muhat[, i])) - 
                         series[,i] + object$muhat[, i]), 
                  binomial = 2*object$model$u[,i]*(series[,i]*log(ifelse(series[,i] == 0, 1, 
                                                                         series[,i]/object$muhat[, i])) + 
                                                     (1 - series[,i])* 
                                                     log(ifelse(series[, i] == 1 | object$muhat[, i]==1, 1,
                                                                (1 - series[,i])/(1 - object$muhat[, i])))), 
                  gamma =
                    -2*(log(ifelse(object$model$y[, i] == 0, 1, 
                                   object$model$y[,  i]/object$muhat[, i])) - 
                          (object$model$y[, i] - object$muhat[, i])/object$muhat[, i]),
                  `negative binomial` = 
                    2*(object$model$y[, i] * log(pmax(1, object$model$y[, i])/object$muhat[, i]) - 
                         (object$model$y[, i] + object$model$u[, i]) * 
                         log((object$model$y[, i] + object$model$u[, i])/(object$muhat[, i] + 
                                                                            object$model$u[, i])))))   
    
    series/sqrt(w*(1 - hatvalues(object)))
  }
  
  state <- function(object) {
    if (is.null(object$etahat)) {
      stop("KFS object needs to contain smoothed estimates of state disturbances eta.")
    } else {
      k <- attr(object$model, "k")
      n <- attr(object$model, "n")
      if (dim(object$model$Q)[3] == 1) {
        z <- which(object$model$Q[, , 1][1 + 0:(k - 1) * (k + 1)] > 0)
        eta <- array(0, c(n, length(z)))
        for (i in 1:(n - 1)) {
          if (!isTRUE(all.equal(object$etahat[i, z], rep(0, length(z))))) {
            x <- try(chol(solve(object$model$Q[z, z, 1] - object$V_eta[z, z, i])) %*% object$etahat[i, z], TRUE)
            if (inherits(x, "try-error")) {
              warning(paste("Could not compute the standardized smoothed state residuals, V_eta[,,", i, "] is not invertible", 
                            sep = ""))
              break
            } else eta[i, ] <- x
          }
        }
        
      } else {
        z <- NULL
        for (i in 1:k) if (sum(object$model$Q[i, i, ]) > 0) 
          z <- c(z, i)
        zlength <- length(z)
        eta <- array(NA, c(n, zlength))
        if (zlength > 1) {
          for (i in 1:(n - 1)) {
            if (!isTRUE(all.equal(object$etahat[i, z][z2], rep(0, length(z2))))) {
              z2 <- which(object$V_eta[z, z, i][1 + 0:(zlength - 1) * (zlength + 1)] > 0)
              x <- try(chol(solve(object$model$Q[z, z, i][z2, z2] - object$V_eta[z, z, i][z2, z2])) %*% object$etahat[i, 
                                                                                                                      z][z2], TRUE)
              if (inherits(x, "try-error")) {
                warning(paste("Could not compute the standardized smoothed state residuals, V_eta[,,", i, "] is not invertible", 
                              sep = ""))
                break
              } else eta[i, z2] <- x
            }
          }
        } else {
          for (i in 1:n) {
            if (!isTRUE(all.equal(object$etahat[i, z], rep(0, length(z))))) 
              eta[i, 1] <- object$etahat[i, z]/sqrt(object$model$Q[z, z, i] - object$V_eta[z, z, i])
          }
        }
      }
      eta[n, ] <- 0
      eta
      
    }
    
  }
  variance<-function(object){
    vars<-object$model$y
    for(i in 1:length(object$model$distribution))
      vars[,i]<-switch(object$model$distribution[i], 
                       gaussian = 1, 
                       poisson = object$muhat[, i],
                       binomial = object$muhat[, i] * (1 - object$muhat[, i])/object$model$u[,i], 
                       gamma = object$muhat[, i]^2, 
                       `negative binomial` = object$muhat[, i] + object$muhat[, i]^2/object$model$u[, i])
    vars
  }
  
  
  return(ts(drop(do.call(type, list(model))),start=start(model$model$y),frequency=frequency(model$model$y),names=colnames(model$model$y)))
} 
