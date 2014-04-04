#' Test whether object is a valid \code{SSModel} object
#'
#' Function \code{is.SSModel} tests whether the object is a valid \code{SSModel} object.
#'
#' Note that the validity of the values in \code{y} and \code{Z} are not tested. 
#' These can contain NA values (but not infinite values),  with condition that when \code{Z[i,,t]} 
#' contains NA value, the corresponding \code{y[t,i]} must also have NA value. In this case 
#' \code{Z[i,,t]} is not referenced in filtering and smoothing, and algorithms works properly. 
#' Note also that this does result NA values in \code{thetahat}, so it could be beneficial to use 
#' for example zeroes in place of NA values in Z, making first sure that the above condition is met.
#'  
#' @export
#' @rdname checkModel
#' @aliases is.SSModel
#' @param object An object to be tested.
#' @param na.check Test the system matrices for NA and infinite values. Default is \code{FALSE}.
#' @param return.logical If \code{FALSE}, error is given if the the model is not a 
#' valid \code{SSModel} object. Otherwise logical value is returned. Defaults to \code{FALSE}.
#' @return Logical value or nothing, depending on the value of \code{return.logical}.
is.SSModel <- function(object, na.check = FALSE, return.logical = TRUE) {
    
    p <- attr(object, "p")
    m <- attr(object, "m")
    k <- attr(object, "k")
    n <- attr(object, "n")
    one <- as.integer(1)
    
    if (return.logical) {
        x <- inherits(object, "SSModel") && all(c("y", "Z", "H", "T", "R", "Q", "a1", "P1", "P1inf", "u", "distribution", "tol", 
            "call") %in% names(object)) && all(object$distribution %in% c("gaussian", "poisson", "binomial", "gamma", "negative binomial")) && 
            identical(dim(object$y), c(n, p)) && (identical(dim(object$Z), c(p, m, n)) || identical(dim(object$Z), c(p, m, one))) && 
            (object$H == "Omitted" || identical(dim(object$H), c(p, p, n)) || identical(dim(object$H), c(p, p, one))) && (identical(dim(object$T), 
            c(m, m, n)) || identical(dim(object$T), c(m, m, one))) && (identical(dim(object$R), c(m, k, n)) || identical(dim(object$R), 
            c(m, k, one))) && (identical(dim(object$Q), c(k, k, n)) || identical(dim(object$Q), c(k, k, one))) && identical(dim(object$a1), 
            c(m, one)) && identical(dim(object$P1), c(m, m)) && identical(dim(object$P1inf), c(m, m)) && (object$u == "Omitted" || 
            identical(dim(object$u), dim(object$y)))
        if (na.check) 
            x <- x && !any(sapply(c("H", "u", "T", "R", "Q", "a1", "P1", "P1inf"), function(x) any(is.na(object[[x]])) | any(is.infinite(object[[x]]))))
        x
    } else {
        if (!inherits(object, "SSModel")) 
            stop("Object is not of class 'SSModel'")
        
        if (!all(object$distribution %in% c("gaussian", "poisson", "binomial", "gamma", "negative binomial"))) 
            stop("The distributions of the observations are not valid. Possible choices are 'gaussian', 'poisson', 'binomial', 'gamma' and ,'negative binomial'.")
        
        if (na.check == TRUE && any(sapply(c("H", "u", "T", "R", "Q", "a1", "P1", "P1inf"), function(x) any(is.na(object[[x]])) | 
            any(is.infinite(object[[x]]))))) 
            stop("System matrices (excluding Z) contain NA or infinite values.")
        
        components <- c("y", "Z", "H", "T", "R", "Q", "a1", "P1", "P1inf", "u", "distribution", "tol", "call")
        if (!all(components %in% names(object))) 
            stop(paste("Model is not a proper object of class 'SSModel'. Following components are missing: ", paste(components[!(components %in% 
                names(object))], collapse = ", ")))
        if (!(identical(dim(object$y), c(n, p)) && (identical(dim(object$Z), c(p, m, n)) || identical(dim(object$Z), c(p, m, one))) && 
            (object$H == "Omitted" || identical(dim(object$H), c(p, p, n)) || identical(dim(object$H), c(p, p, one))) && (identical(dim(object$T), 
            c(m, m, n)) || identical(dim(object$T), c(m, m, one))) && (identical(dim(object$R), c(m, k, n)) || identical(dim(object$R), 
            c(m, k, one))) && (identical(dim(object$Q), c(k, k, n)) || identical(dim(object$Q), c(k, k, one))) && identical(dim(object$a1), 
            c(m, one)) && identical(dim(object$P1), c(m, m)) && identical(dim(object$P1inf), c(m, m)) && (object$u == "Omitted" || 
            identical(dim(object$u), dim(object$y))))) 
            stop("Model is not a proper object of class 'SSModel'. Check dimensions of system matrices.")
    }
    
} 
