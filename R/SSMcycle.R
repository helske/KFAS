#' @rdname SSModel
#' @export
SSMcycle <- function(period, Q, type, index, a1, P1, P1inf, n = 1, ynames) {
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
  if (!(length(period) == 1 & period > 0))
    stop("Period of the cycle component must be larger than 0. ")
  lambda <- 2 * pi/period
  m <- 2 * ((p - 1) * (type == 1) + 1)
  Z <- matrix(0, p, m)
  T <- matrix(0, m, m)
  Z_univariate <- matrix(c(1, 0), 1, 2)
  T_univariate <- matrix(c(cos(lambda), -sin(lambda), sin(lambda), cos(lambda)),
    2, 2)
  if (type != 2) {
    for (i in 1:p) {
      Z[i, ((i - 1) * 2 + 1):(i * 2)] <- Z_univariate
      T[((i - 1) * 2 + 1):(i * 2), ((i - 1) * 2 + 1):(i * 2)] <- T_univariate
    }
  } else {
    Z <- matrix(Z_univariate, nrow = p, ncol = m, byrow = TRUE)
    T <- T_univariate
  }
  state_names <- paste0(c("cycle", "cycle*"), rep(ynames, each = 2))
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
      if (length(Q) != 1 && (is.null(dim(Q)) || any(dim(Q)[1:2] != p) || !(max(dim(Q)[3], 1, na.rm = TRUE) %in% c(1, n))))
        stop("Misspecified Q, argument Q must be (p x p) matrix, (p x p x 1), or (p x p x n) array where m is the number of time series.")
      tvq <- max(dim(Q)[3] == n, 0, na.rm = TRUE)
      Qm <- array(0, c(m, m, tvq * (n - 1) + 1))
      if (tvq) {
        for (i in 1:(tvq * (n - 1) + 1)) Qm[cbind(rep(1:(p * 2), p), rep(1:2,
          p^2) + rep(0:(p - 1) * 2, each = p * 2), i)] <- rep(Q[, , i], each = 2)
      } else Qm[cbind(rep(1:(p * 2), p), rep(1:2, p^2) + rep(0:(p - 1) * 2, each = p *
          2), 1)] <- rep(Q, each = 2)
    } else {
      if (length(Q) != 1 && (is.null(dim(Q)) || any(dim(Q)[1:2] != 1) || !(max(dim(Q)[3], 1, na.rm = TRUE) %in% c(1, n))))
        stop("Misspecified Q, argument Q must be a scalar, (1 x 1) matrix, or (1 x 1 x 1)/(1 x 1 x n) array.")
      tvq <- max(dim(Q)[3] == n, 0, na.rm = TRUE)
      Qm <- array(0, c(m, m, tvq * (n - 1) + 1))
      if (tvq) {
        for (i in 1:(tvq * (n - 1) + 1)) Qm[cbind(1:2, 1:2, i)] <- rep(Q[1,
          1, i], 2)
      } else Qm[cbind(1:2, 1:2, 1)] <- rep(Q, 2)
    }
    k <- dim(Qm)[1]
    R <- diag(k)
  }
  list(index = index, m = m, k = k, Z = Z, T = T, R = R, Q = Qm, a1 = a1, P1 = P1,
    P1inf = P1inf, tvq = tvq, tvr = 0, tvz = 0, state_names = state_names, period = period)
}
