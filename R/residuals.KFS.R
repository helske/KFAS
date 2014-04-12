#' Extract Residuals of KFS output
#' 
#' @details For object of class KFS, several types of residuals can be computed:
#' 
#' \itemize{
#' \item 'recursive': One-step ahead prediction residuals
#' \deqn{v_{t,i}),}
#' with residuals being undefined in diffuse phase. Only supported for fully Gaussian models.
#' 
#' \item 'response': Data minus fitted values, \eqn{y-E(y)}{y-E(y)}. 
#' 
#' \item 'pearson':  \deqn{(y_{t,i}-\theta_{t,i})/\sqrt{V(\mu)_{t,i}}, \quad i=1,\ldots,p,t=1,\ldots,n,}{(y[t,i]-\theta[t,i])V(\mu)[t,i]^(-0.5), i=1,\ldots,p, t=1,\ldots,n,}
#'                   where \eqn{V(\mu_{t,i})}{V(\mu[t,i])} is the variance function of the model.
#'
#' \item 'state':  Residuals based on the smoothed disturbance terms \eqn{\eta} are defined as
#' \deqn{\hat \eta_t, \quad t=1,\ldots,n,}{L^{-1}[t] \eta[t], t=1,\ldots,n}.
#' 
#' \item 'deviance': Deviance residuals.
#' }
#' @S3method residuals KFS
#' @method residuals KFS
#' @param object KFS object
#' @param type Character string defining the type of residuals.
#' @param ... Ignored.

residuals.KFS <- function(object, type = c("recursive","deviance", "pearson", "response","state"), ...) {
  
  type <- match.arg(type)
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
  
  if ((type == "recursive" || type == "state") && any(object$model$distribution != "gaussian")) 
    stop("Recursive and state residuals are only supported for fully gaussian models.")
  
  series <- 
    switch(type,
           recursive = {
             if (is.null(object[["a",exact=TRUE]])) 
               stop("KFS object needs to contain filtered estimates of states. ")
             series <- object$v
             series[1:(object$d - 1), ] <- NA
             series[object$d, 1:object$j] <- NA
             series
           },
           response = {
             series <- object$model$y
             if(sum(bins<-object$model$distribution=="binomial")>0)
               series[, bins]<-
               series[, bins]/
               object$model$u[, bins]    
             series-fitted(object)
           },
           state = {
             if (is.null(object$etahat)) {
               stop("KFS object needs to contain smoothed estimates of state disturbances eta.")
             } else {
               object$etahat
             }           
           },
           pearson = {
             series <- object$model$y
             if(sum(bins<-object$model$distribution=="binomial")>0)
               series[, bins]<-
               series[, bins]/
               object$model$u[, bins]           
             (series - fitted(object))/sqrt(variance(object))
           },
           deviance = {  
             series <- object$model$y
             if(sum(bins<-object$model$distribution=="binomial")>0)
               series[, bins]<-
               series[, bins]/
               object$model$u[, bins]
             
             for (i in 1:attr(object$model, "p")) 
               series[, i] <- ifelse(series[, i]>object$muhat[, i],1,-1)*
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
                             -2*(log(ifelse(series[, i] == 0, 1, 
                                            series[,  i]/object$muhat[, i])) - 
                                   (series[, i] - object$muhat[, i])/object$muhat[, i]),
                           `negative binomial` = 
                             2*(series[, i] * log(pmax(1, series[, i])/object$muhat[, i]) - 
                                  (series[, i] + object$model$u[, i]) * 
                                  log((series[, i] + object$model$u[, i])/(object$muhat[, i] + 
                                                                                     object$model$u[, i])))))   
             
             series
           })
 ts(drop(series),start=start(object$model$y),frequency=frequency(object$model$y),names=colnames(object$model$y))
} 
