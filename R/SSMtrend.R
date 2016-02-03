#' @rdname SSModel
#' @export
SSMtrend <- function(degree = 1, Q, type, index, a1, P1, P1inf, n = 1, ynames) {
  if (missing(index))
    index <- 1
  p <- length(index)
  if (!missing(ynames) && !is.null(ynames)) {
    ynames <- paste0(".", ynames)
  } else ynames <- ""
  if (missing(type)) {
    type <- 1L
  } else {
    type <- pmatch(x = type, table = c("distinct", "common"))
    if (is.na(type))
      stop("type must be 'distinct' or 'common'.")
  }
  if (!(length(degree) == 1 & degree > 0 & abs(degree - round(degree)) == 0))
    stop("Degree of the trend component must be positive integer. ")
  m <- ((p - 1) * (type == 1) + 1) * degree
  Z <- matrix(0, p, m)
  T <- matrix(0, m, m)
  if (type == 2) {
    Z[, 1] <- 1
    p <- 1
  } else {
    for (i in 1:p) Z[i, (i - 1) * degree + 1] <- 1
  }
  state_names <- switch(as.character(degree), 
    "1" = paste0("level", ynames), 
    "2" = paste0(c("level", "slope"), rep(ynames, each = degree)), 
    paste0("trend", rep(1:degree), rep(ynames, each = degree)))
  dxm <- 1 + 0:(m - 1) * (m + 1)
  T[dxm] <- 1
  if (degree > 1)
    T[dxm[-m] + m] <- rep(c(rep(1, degree - 1), 0), length = m - 1)
  if (missing(a1)) {
    a1 <- matrix(0, m, 1)
  } else {
    if (length(a1) != m || any(dim(a1) != c(m, 1)))
      stop("a1 must be a (m x 1) matrix where m is the number of states. ")
    a1 <- matrix(a1, m, 1)
  }
  if (missing(P1)) {
    P1 <- matrix(0, m, m)
  } else {
    if (length(P1) > 1 && any(dim(P1) != m))
      stop("P1 must be a (m x m) matrix where m is the number of states. ")
    P1 <- matrix(P1, m, m)
  }
  if (missing(P1inf)) {
    P1inf <- diag(m)
  } else {
    if (length(P1inf) > 1 && any(dim(P1inf) != m))
      stop("P1inf must be a (m x m) diagonal matrix where m is the number of states. ")
    P1inf <- matrix(P1inf, m, m)
  }
  diag(P1inf)[diag(P1) > 0 || is.na(diag(P1))] <- 0
  if (missing(Q)) {
    k <- 0
    Qm <- R <- NULL
    tvq <- 0
  } else {
    if (type == 1) {
      if (!is.list(Q)) {
        if (degree > 1) {
          stop("Q must be a list of length degree, which contains (p x p) matrices, (p x p x 1), or (p x p x n) arrays, where p is the number of series. ")
        } else Q <- list(Q)
      }
      tvq <- max(unlist(sapply(lapply(Q, dim), "[", 3)) > 1, 0, na.rm = TRUE)
      Qm <- array(0, c(m, m, tvq * (n - 1) + 1))
      if (!is.list(Q))
        stop("Q must be a list of length degree.")
      for (i in 1:degree) {
        if(length(Q[[i]]) != 1 && (is.null(dim(Q[[i]])) || any(dim(Q[[i]])[1:2]!=p) || !(max(dim(Q[[i]])[3], 1, na.rm = TRUE) %in% c(1, n))))
          stop("Q must be a list of length degree, which contains (p x p) matrices, (p x p x 1), or (p x p x n) arrays, where p is the number of series. ")
        Qm[seq(from = i, by = degree, length = p), seq(from = i, by = degree, length = p), ] <- Q[[i]]
      }
      k <- dim(Qm)[1]
      R <- diag(k)
    } else {
      if (is.list(Q) || (length(Q) != degree && is.null(dim(Q))) || (any(dim(Q)[1:2] !=
          degree) || !(max(1, dim(Q)[3], na.rm = TRUE) %in% c(1, n, NA))))
        stop("Misspecified Q, argument Q must be a vector of length d, (d x d) matrix, or (d x d x 1)/(d x d x n) array where d is the degree of the trend.")
      if (length(Q) == degree)
        Q <- diag(drop(Q), degree)
      tvq <- max(dim(Q)[3] == n, 0, na.rm = TRUE)
      Qm <- Q
      k <- dim(Qm)[1]
      R <- diag(k)
    }
  }
  list(index = index, m = m, k = k, Z = Z, T = T, R = R, Q = Qm, a1 = a1, P1 = P1,
    P1inf = P1inf, tvq = tvq, tvr = 0, tvz = 0, state_names = state_names)
}
