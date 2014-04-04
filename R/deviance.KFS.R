#' Deviance of a State Space Model
#' 
#' Returns the deviance of a object of class \code{KFS}.
#' @S3method deviance KFS
#' @method deviance KFS
#' @param object An object of class \code{KFS}.
#' @param \dots Ignored.
#' @return The value of the deviance extracted from object.
deviance.KFS<-function(object,...){
 sum(residuals(object,type="deviance")^2,na.rm=TRUE)
}