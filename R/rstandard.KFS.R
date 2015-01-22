#' Extract Standardized Residuals from KFS output
#' @export
#' @details For object of class KFS with fully Gaussian observations, several 
#'   types of standardized residuals can be computed:
#'   
#'   \item 'recursive': For Gaussian models the one-step ahead prediction
#'   residuals defined as 
#'   \deqn{v_{t,i})/\sqrt{F_{i,t}},}{v[t,i])/\sqrt{F[i,t]},} with residuals 
#'   being undefined in diffuse phase. For non-Gaussian models recursive
#'   residuals are obtained as 
#'   \deqn{L^{-1}_t(y_{t}-\mu_{t}),}{L^(-1)[t](y[t]-\mu[t]),} where 
#'   \eqn{L_t}{L[t]} is the lower triangular matrix from Cholesky decomposition 
#'   of \eqn{V(y_t)+V(\mu_t)}, and both the variance function V(y_t) and
#'   V(\mu_t) are based on the filtered means \mu_t$. Computing these for large
#'   non-Gaussian models can be time consuming as filtering is needed.
#'   
#'   \item 'observation':  Residuals based on the smoothed observation 
#'   disturbance terms \eqn{\epsilon} are defined as \deqn{L^{-1}_t \hat 
#'   \epsilon_t, \quad t=1,\ldots,n,}{L^{-1}[t] \epsilon[t], t=1,\ldots,n,} 
#'   where \eqn{L_t}{L[t]} is the lower triangular matrix from Cholesky 
#'   decomposition of \eqn{V_{\epsilon,t}}{V[\epsilon,t]}.
#'   
#'   
#'   \item 'state':  Residuals based on the smoothed state disturbance terms 
#'   \eqn{\eta} are defined as \deqn{L^{-1}_t \hat \eta_t, \quad 
#'   t=1,\ldots,n,}{L^{-1}[t] \eta[t], t=1,\ldots,n,} where \eqn{L_t}{L[t]} is 
#'   the lower triangular matrix from Cholesky decomposition of 
#'   \eqn{V_{\eta,t}}{V[\eta,t]}.
#'   
#'   
#'   \item 'pearson':  Standardized Pearson residuals 
#'   \deqn{L^{-1}_t(y_{t}-\theta_{i}), \quad 
#'   t=1,\ldots,n,}{L^(-1)[t](y[t]-\theta[t]), t=1,\ldots,n,}, where 
#'   \eqn{L_t}{L[t]} is the lower triangular matrix from Cholesky decomposition 
#'   of \eqn{V(y_t)-V(\mu_t)}. For gaussian models, these coincide with the 
#'   standardized smoothed \eqn{\epsilon} disturbance residuals, and for 
#'   generalized linear models these coincide with the standardized Pearson 
#'   residuals (hence the name).
#'   
#'   

#' \item 'deviance': Deviance residuals. 
#' Deprecated. This option was meant to be used only for the GLM comparisons, as their generalization to other models is lacking,
#' but these will be completely removed in future in order to avoid misleading results in non-GLM settings.
#' }
#'
#' @param model KFS object
#' @param type Type of residuals. See details.
#' @param ... Ignored.
rstandard.KFS <- 
  function(model, type = c("recursive", "observation","state", "deviance", "pearson"), 
           ...) {
    type <- match.arg(type)
    if(type=="deviance")
      .Deprecated(msg="Argument type=\"deviance\" is deprecated.")
    
    if (type == "state" && any(model$model$distribution != "gaussian")) 
      stop("State residuals are only supported for fully gaussian models.")
    
    recursive <- function(object) {
      if(any(object$model$distribution !=  "gaussian") && is.null(object[["m", exact = TRUE]]))
        stop("KFS object does not contain filtered means. ")
      if (all(object$model$distribution ==  "gaussian") && is.null(object[["v", exact = TRUE]])) 
        stop("KFS object does not contain prediction errors v. ")
      if(all(object$model$distribution ==  "gaussian")){
        series <- object$v/sqrt(t(object$F))
        if (object$d > 0) {
          series[1:(object$d - 1), ] <- NA
          series[object$d, 1:object$j] <- NA
        } 
      }else {
        
        dj<-KFS(approxSSM(object$model),filtering="state",smoothing="none")[c("d","j")]
        
        p <- attr(object$model, "p")  
        series<-object$model$y
        
        variance <- function(object) {
          vars <- object$model$y
          for (i in 1:length(object$model$distribution)){ 
            vars[, i] <- switch(object$model$distribution[i], 
                                gaussian = object$u[,i], 
                                poisson = object$m[, i], 
                                binomial = object$m[, i] * (1 - object$m[, i])/object$model$u[, i], 
                                gamma = object$m[, i]^2/object$model$u[,i], 
                                `negative binomial` = object$m[, i] + object$m[, i]^2/object$model$u[, i])
          }
          vars
        }
        vars<-variance(object)     
        if (sum(bins <- object$model$distribution == "binomial") > 0) 
          series[, bins] <- series[, bins]/object$model$u[, bins]
        series<- series-object[["m", exact = TRUE]]     
        for(t in (dj$d+1):attr(object$model, "n") ){
          yobs<-which(!is.na(series[t,]))
          if(length(yobs)>0){
            tmp <- chol(diag(vars[t,yobs],p)+object$P_mu[yobs,yobs,t])
            series[t,yobs] <- series[t,yobs]%*%solve(tmp)
          }
        }
        if(dj$d > 0){
          series[1:(dj$d - 1), ] <- NA
          series[dj$d, 1:dj$j] <- NA                           
        }
      }
      series
    }
    
    pearson <- function(object) {     
      n <- attr(object$model, "n")
      
      if(is.null(object$muhat))
        stop("KFS object needs to contain smoothed means. ")
      
      if (all(object$model$distribution == "gaussian")) {              
        
        tv<-dim(object$model$H)[3] > 1
        res<-object$model$y-object$muhat 
        for(t in 1:n){
          yobs<-which(!is.na(res[t,]))
          if(length(yobs)>0){
            tmp <- chol(object$model$H[yobs,yobs,(t-1)*tv+1]-object$V_mu[yobs,yobs,t])
            res[t,yobs] <- res[t,yobs]%*%solve(tmp)
          }
        }
      }else {    
        p <- attr(object$model, "p")  
        series<-object$model$y
        vars<-variance(object)     
        if (sum(bins <- object$model$distribution == "binomial") > 0) 
          series[, bins] <- series[, bins]/object$model$u[, bins]
        res<-(series-object$muhat)       
        for(t in 1:n){
          yobs<-which(!is.na(res[t,]))
          if(length(yobs)>0){
            tmp <- chol(diag(vars[t,yobs],p)-object$V_mu[yobs,yobs,t])
            res[t,yobs] <- res[t,yobs]%*%solve(tmp)
          }
        }        
      }
      res
    }
    deviance <- function(object) {
      if (all(object$model$distribution == "gaussian")) {
        w <- matrix(apply(object$model$H, 3, diag), attr(object$model, "n"), 
                    attr(object$model, "p"), byrow = TRUE)
      } else {
        w <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
        for (i in 1:attr(object$model, "p")) w[, i] <- switch(object$model$distribution[i], 
                                                              gaussian = object$model$u[, i], poisson = 1, binomial = 1, gamma = 1/object$model$u[, 
                                                                                                                                                  i], `negative binomial` = 1)
      }
      series <- object$model$y
      if (sum(bins <- object$model$distribution == "binomial") > 0) 
        series[, bins] <- series[, bins]/object$model$u[, bins]
      for (i in 1:attr(object$model, "p")) series[, i] <- ifelse(series[, i] > 
                                                                   object$muhat[, i], 1, -1) * sqrt(switch(object$model$distribution[i], 
                                                                                                           gaussian = (series[, i] - object$muhat[, i])^2, poisson = 2 * (series[, 
                                                                                                                                                                                 i] * log(ifelse(series[, i] == 0, 1, series[, i]/object$muhat[, i])) - 
                                                                                                                                                                            series[, i] + object$muhat[, i]), binomial = 2 * object$model$u[, 
                                                                                                                                                                                                                                            i] * (series[, i] * log(ifelse(series[, i] == 0, 1, series[, i]/object$muhat[, 
                                                                                                                                                                                                                                                                                                                         i])) + (1 - series[, i]) * log(ifelse(series[, i] == 1 | object$muhat[, 
                                                                                                                                                                                                                                                                                                                                                                                               i] == 1, 1, (1 - series[, i])/(1 - object$muhat[, i])))), gamma = -2 * 
                                                                                                             (log(ifelse(object$model$y[, i] == 0, 1, object$model$y[, i]/object$muhat[, 
                                                                                                                                                                                       i])) - (object$model$y[, i] - object$muhat[, i])/object$muhat[, 
                                                                                                                                                                                                                                                     i]), `negative binomial` = 2 * (object$model$y[, i] * log(pmax(1, 
                                                                                                                                                                                                                                                                                                                    object$model$y[, i])/object$muhat[, i]) - (object$model$y[, i] + 
                                                                                                                                                                                                                                                                                                                                                                 object$model$u[, i]) * log((object$model$y[, i] + object$model$u[, 
                                                                                                                                                                                                                                                                                                                                                                                                                                  i])/(object$muhat[, i] + object$model$u[, i])))))
      series/sqrt(w * (1 - hatvalues(object)))
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
              x <- try(chol(solve(object$model$Q[z, z, 1] - object$V_eta[z,z, i])) %*% 
                         object$etahat[i, z], TRUE)
              if (inherits(x, "try-error")) {
                warning(paste("Could not compute the standardized smoothed state residuals, V_eta[,,",i, "] is not invertible", sep = ""))
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
                x <- try(chol(solve(object$model$Q[z, z, i][z2, z2] - 
                                      object$V_eta[z, z, i][z2, z2])) %*% object$etahat[i, z][z2], TRUE)
                if (inherits(x, "try-error")) {
                  warning(paste("Could not compute the standardized smoothed state residuals, V_eta[,,", i, "] is not invertible", sep = ""))
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
    variance <- function(object) {
      vars <- object$model$y
      for (i in 1:length(object$model$distribution)){ 
        vars[, i] <- switch(object$model$distribution[i], 
                            gaussian = object$u[,i], 
                            poisson = object$muhat[, i], 
                            binomial = object$muhat[, i] * (1 - object$muhat[, i])/object$model$u[, i], 
                            gamma = object$muhat[, i]^2/object$model$u[,i], 
                            `negative binomial` = object$muhat[, i] + object$muhat[, i]^2/object$model$u[, i])
      }
      vars
    }
    x<-do.call(type, list(model))
    return(ts(drop(x), start = start(model$model$y), frequency = frequency(model$model$y), 
              names = colnames(model$model$y)))
  } 
