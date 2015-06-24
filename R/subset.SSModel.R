#' @method subset<- SSModel
#' @export
#' @rdname Extract.SSModel
`subset<-.SSModel` <- function(x, element, states, etas, series, times, 
  ..., value) {
  
  .Deprecated(new="[<-.SSModel")
  x[element, states, etas, series, times, ...] <- value
  x
}
#' @rdname Extract.SSModel
#' @export
`subset<-` <- function(x, ..., value) { 
  .Deprecated(new="[<-")
  UseMethod("subset<-")   
}

#' @method subset SSModel
#' @export 
#' @rdname Extract.SSModel
#' @param ... ignored.
subset.SSModel <- function(x, element, states, etas, series, times, ...) {  
  .Deprecated(new="[.SSModel")
  x[element, states, etas, series, times, ...]  
  
} 
