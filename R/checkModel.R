#' Test whether object is a valid \code{SSModel} object
#'
#' Function \code{is.SSModel} tests whether the object is a valid \code{SSModel}
#' object.
#'
#' Note that the validity of the values in \code{y} and \code{Z} are not tested.
#' These can contain NA values (but not infinite values),  with condition that
#' when \code{Z[i,,t]} contains NA value, the corresponding \code{y[t,i]} must
#' also have NA value. In this case \code{Z[i,,t]} is not referenced in
#' filtering and smoothing, and algorithms works properly.
#'
#'
#' @export
#' @rdname checkModel
#' @aliases is.SSModel
#' @param object An object to be tested.
#' @param na.check Test the system matrices for \code{NA} and infinite values.
#' Also checks for large values (> 1e7) in covariance matrices \code{H}
#' and \code{Q} which could cause large rounding errors in filtering.
#' Positive semidefiniteness of these matrices is not checked. Default
#'   is \code{FALSE}.
#' @param return.logical If \code{FALSE} (default), an error is given if the the model is not
#'   a valid \code{SSModel} object. Otherwise logical value is returned.
#' @return Logical value or nothing, depending on the value of
#'   \code{return.logical}.
#' @examples
#' model <- SSModel(rnorm(10) ~ 1)
#' is.SSModel(model)
#' model['H'] <- 1
#' is.SSModel(model)
#' model$H[] <- 1
#' is.SSModel(model)
#' model$H[,,1] <- 1
#' is.SSModel(model)
#' model$H <- 1
#' is.SSModel(model)
is.SSModel <- function(object, na.check = FALSE, return.logical = TRUE) {

  tol <- 1e7
  p <- attr(object, "p")
  m <- attr(object, "m")
  k <- attr(object, "k")
  n <- attr(object, "n")
  tv <- unname(attr(object, "tv"))

  if (return.logical) {
    x <- inherits(object, "SSModel") &&
      all(c("y", "Z", "H", "T", "R", "Q", "a1","P1", "P1inf", "u",
            "distribution", "tol", "call") %in% names(object)) &&
      all(object$distribution %in%
            c("gaussian", "poisson", "binomial", "gamma",
              "negative binomial")) &&
      is.integer(c(p,m,k,n,tv)) &&
      identical(dim(object$y), c(n, p)) &&
      identical(dim(object$Z), c(p, m, (n - 1L) * tv[1] + 1L)) &&
      (identical(object$H, "Omitted") ||
         identical(dim(object$H), c(p, p, (n - 1L) * tv[2] + 1L))) &&
      identical(dim(object$T), c(m, m, (n - 1L) * tv[3] + 1L)) &&
      identical(dim(object$R), c(m, k, (n - 1L) * tv[4] + 1L)) &&
      identical(dim(object$Q), c(k, k, (n - 1L) * tv[5] + 1L)) &&
      identical(dim(object$a1), c(m, 1L)) &&
      identical(dim(object$P1), c(m, m)) &&
      identical(dim(object$P1inf), c(m, m)) &&
      (identical(object$u, "Omitted") ||
         identical(dim(object$u), dim(object$y))) &&
      all(diag(object$P1inf) %in% c(0,1)) &&
      all(object$P1inf[col(diag(m)) != row(diag(m))] == 0)
    if (na.check) {
      x <- x && !any(sapply(c("H", "u", "T", "R", "Q", "a1", "P1", "P1inf"),
                            function(x) any(is.na(object[[x]])) ||
                              any(is.infinite(object[[x]])))) &&
        max(object$Q) <= tol &&
        ifelse(identical(object$u, "Omitted"), max(object$H) <= tol, TRUE)
    }
    x
  } else {
    if (!inherits(object, "SSModel")) {
      stop("Object is not of class 'SSModel'")
    }

    if (!is.integer(c(p, m, k, n, tv))) {
      stop(paste0("Storage mode of some of the model attributes 'p', 'k', ",
                  "'m', 'n', 'tv' is not integer."))
    }

    components <- c("y", "Z", "H", "T", "R", "Q", "a1", "P1", "P1inf", "u",
                    "distribution", "tol", "call")
    if (!all(components %in% names(object))) {
      stop(paste("Model is not a proper object of class 'SSModel'.
                   Following componentsare missing: ",
                 paste(components[!(components %in% names(object))],
                       collapse = ", ")))
    }
    if (!all(object$distribution %in%
             c("gaussian", "poisson", "binomial", "gamma",
               "negative binomial"))) {
      stop(paste0("The distributions of the observations are not valid. ",
                  "Possible choices are 'gaussian', 'poisson', 'binomial', ",
                  "'gamma' and ,'negative binomial'."))
    }

    if (!(identical(dim(object$y), c(n, p)) &&
          identical(dim(object$Z), c(p, m, (n - 1L) * tv[1] + 1L)) &&
          (identical(object$H, "Omitted") ||
           identical(dim(object$H), c(p, p, (n - 1L) * tv[2] + 1L))) &&
          identical(dim(object$T), c(m, m, (n - 1L) * tv[3] + 1L)) &&
          identical(dim(object$R), c(m, k, (n - 1L) * tv[4] + 1L)) &&
          identical(dim(object$Q), c(k, k, (n - 1L) * tv[5] + 1L)) &&
          identical(dim(object$a1), c(m, 1L)) &&
          identical(dim(object$P1), c(m, m)) &&
          identical(dim(object$P1inf), c(m, m)) &&
          (identical(object$u, "Omitted") || identical(dim(object$u),
                                                       dim(object$y))))) {
      stop(paste0("Model is not a proper object of class 'SSModel'. ",
                  "Check dimensions of system matrices."))
    }
    if (na.check == TRUE &&
        (any(sapply(c("H", "u", "T", "R", "Q", "a1", "P1", "P1inf"),
                    function(x) any(is.na(object[[x]])) ||
                    any(is.infinite(object[[x]])))) ||
         max(object$Q) > tol || ifelse(identical(object$u, "Omitted"),
                                       max(object$H) > tol, FALSE))) {
      stop(paste0("System matrices (excluding Z) contain NA or infinite ",
                  "values, covariance matrices contain values larger ",
                  "than ",tol))
    }
    if (!all(diag(object$P1inf) %in% c(0,1)) ||
        !all(object$P1inf[col(diag(m)) != row(diag(m))] == 0)) {
      stop(paste0("Matrix P1inf is not a diagonal matrix with zeros and ",
                  "ones on diagonal."))
    }
  }
}
