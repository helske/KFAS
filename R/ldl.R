#' LDL Decomposition of a Matrix
#'
#' Function \code{ldl} computes the LDL decomposition of a positive semidefinite matrix.
#'
#' @export
#' @param x Symmetrix matrix.
#' @param tol Tolerance parameter for LDL decomposition, determines which
#' diagonal values are counted as zero. Same value is used in isSymmetric function.
#' @return Transformed matrix with D in diagonal, L in strictly lower diagonal 
#' and zeros on upper diagonal.
ldl <- 
  function(x, tol = max(abs(diag(x))) * .Machine$double.eps) {
    if (!isSymmetric(x, tol = tol)) 
        stop("Matrix is not symmetric!")
    out <- .Fortran(fldl, x = x, as.integer(dim(x)[1]), tol = tol, info = integer(1))
    out$x
} 
