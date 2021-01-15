#' LDL Decomposition of a Matrix
#'
#' Function \code{ldl} computes the LDL decomposition of a positive semidefinite matrix.
#'
#' @export
#' @param x Symmetrix matrix.
#' @param tol Tolerance parameter for LDL decomposition, determines which
#' diagonal values are counted as zero. Same value is used in isSymmetric function. 
#' Default is \code{max(100, max(abs(diag(as.matrix(x))))) * .Machine$double.eps}.
#' @return Transformed matrix with D in diagonal, L in strictly lower diagonal
#' and zeros on upper diagonal.
#' @examples
#' # Positive semidefinite matrix, example matrix taken from ?chol
#' x <- matrix(c(1:5, (1:5)^2), 5, 2)
#' x <- cbind(x, x[, 1] + 3*x[, 2])
#' m <- crossprod(x)
#' l <- ldl(m, tol = 1e-8) # arm64 Mac setup in CRAN fails with default tolerance
#' d <- diag(diag(l))
#' diag(l) <- 1
#' all.equal(l %*% d %*% t(l), m, tol = 1e-15)
ldl <- function(x, tol) {
  
  if(length(x) == 1) {
    if(x >= 0) return(matrix(x))
    if(x < 0) stop("Matrix x is not positive semidefinite.")
  }
  
  if (missing(tol)) {
    tol <- max(100, max(abs(diag(as.matrix(x))))) * .Machine$double.eps
  }
  if (!isSymmetric(x, tol = tol))
    stop("Matrix is not symmetric!")
  out <- .Fortran(fldl, x = x, as.integer(dim(x)[1]), tol = tol, info = integer(1))
  if (out$info != 0)
    stop("Matrix x is not positive semidefinite.")
  out$x
}
