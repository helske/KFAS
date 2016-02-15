#' @rdname SSModel
#' @export
SSMcustom <- function(Z, T, R, Q, a1, P1, P1inf, index, n = 1) {
  if (missing(index))
    index <- 1
  p <- length(index)
  if (length(Z) == 1 && p == 1) {
    dim(Z) <- c(1, 1, 1)
    m <- 1
  } else {
    if ((length(Z) == 1) || !(dim(Z)[1] == p) || !dim(Z)[3] %in% c(1, NA, n))
      stop("Misspecified Z, argument Z must be a (p x m) matrix, (p x m x 1) or (p x m x n) array, where p is the number of time series, m is the number of states.")
    m <- dim(Z)[2]
    dim(Z) <- c(p, m, (n - 1) * (max(dim(Z)[3], 0, na.rm = TRUE) > 1) + 1)
  }
  if (length(T) == 1 && m == 1) {
    dim(T) <- c(1, 1, 1)
  } else {
    if ((length(T) == 1) || any(dim(T)[1:2] != m) || !dim(T)[3] %in% c(1, NA, n))
      stop("Misspecified T, argument T must be a (m x m) matrix, (m x m x 1) or (m x m x n) array, where m is the number of states.")
    dim(T) <- c(m, m, (n - 1) * (max(dim(T)[3], 0, na.rm = TRUE) > 1) + 1)
  }
  if (length(Q) == 1) {
    dim(Q) <- c(1, 1, 1)
    k <- 1
  } else {
    if (!identical(dim(Q)[1], dim(Q)[2]) || dim(Q)[1] > m || !dim(Q)[3] %in%
        c(1, NA, n))
      stop("Misspecified Q, argument Q must be a (k x k) matrix, (k x k x 1) or (k x k x n) array, where k<=m is the number of disturbances eta, and m is the number of states.")
    k <- dim(Q)[1]
    dim(Q) <- c(k, k, (n - 1) * (max(dim(Q)[3], 0, na.rm = TRUE) > 1) + 1)
  }
  if (missing(R)) {
    R <- diag(m)[, 1:k, drop = FALSE]
    dim(R) <- c(m, k, 1)
  } else {
    if (all(c(length(R), k, m) == 1)) {
      dim(R) <- c(1, 1, 1)
    } else {
      if ((length(R) == 1) || !(dim(R)[1] == m) || dim(R)[2] != k || !dim(R)[3] %in% c(1, NA, n))
        stop("Misspecified R, argument R must be a (m x k) matrix, (m x k x 1) or (m x k x n) array, where k<=m is the number of disturbances eta, and m is the number of states.")
      dim(R) <- c(m, k, (n - 1) * (max(dim(R)[3], 0, na.rm = TRUE) > 1) + 1)
    }
  }
  if (missing(a1)) {
    a1 <- matrix(0, m, 1)
  } else {
    if (length(a1) <= m) {
      a1 <- matrix(a1, m, 1)
    } else stop("Misspecified a1, argument a1 must be a vector of length m, or (m x 1) matrix, where m is the number of state_names and 1<=t<=m.")
  }
  if (missing(P1)) {
    P1 <- matrix(0, m, m)
  } else {
    if (length(P1) == 1 && m == 1) {
      dim(P1) <- c(1, 1)
    } else {
      if (any(dim(P1)[1:2] != m))
        stop("Misspecified P1, argument P1 must be (m x m) matrix, where m is the number of states. ")
    }
  }
  if (missing(P1inf)) {
    P1inf <- matrix(0, m, m)
  } else {
    if (length(P1inf) == 1 && m == 1) {
      dim(P1inf) <- c(1, 1)
    } else {
      if (any(dim(P1inf)[1:2] != m))
        stop("Misspecified P1inf, argument P1inf must be a (m x m) matrix, where m is the number of states..")
    }
  }
  diag(P1inf)[diag(P1) > 0 || is.na(diag(P1))] <- 0
  state_names <- paste0("custom", 1:m)
  list(index = index, m = m, k = k, p = p, n = n, Z = Z, T = T, R = R, Q = Q, a1 = a1,
    P1 = P1, P1inf = P1inf, tvz = dim(Z)[3] > 1, tvt = dim(T)[3] > 1, tvr = dim(R)[3] >
      1, tvq = dim(Q)[3] > 1, state_names = state_names)
}
