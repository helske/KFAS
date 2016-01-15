#' @rdname SSModel
#' @seealso \code{artransform}
#' @export
SSMarima <- function(ar = NULL, ma = NULL, d = 0, Q, stationary = TRUE,
  index, n = 1, ynames) {

  # Q is either a p times p matrix or scalar (in univariate case)
  if (!is.null(ar) && stationary && !all(Mod(polyroot(c(1, -ar))) > 1))
    stop("ARIMA part is non-stationary.")
  if (missing(index))
    index <- 1
  p <- length(index)
  if (!missing(ynames) && !is.null(ynames)) {
    ynames <- paste0(".", ynames)
  } else ynames <- ""
  if (missing(Q)) {
    Q <- diag(p)
  } else {
    if (length(Q) == 1)
      Q <- matrix(Q)
    if (any(dim(Q)[1:2] != p) || length(dim(Q)) > 2)
      stop("Misspecified Q, argument Q must be (p x p) matrix where p is the number of series.")
  }
  ar_length <- length(ar)
  ma_length <- length(ma)
  d <- max(d, 0)
  m1 <- max(ar_length, ma_length + 1) + d
  k <- p
  Z_univariate <- matrix(0, 1, m1)
  T_univariate <- P1inf_univariate <- matrix(0, m1, m1)
  R_univariate <- matrix(0, m1, 1)
  Z_univariate[1, 1:(d + 1)] <- 1
  if (d > 0) {
    T_univariate[1:d, 1:d][upper.tri(T_univariate[1:d, 1:d], diag = TRUE)] <- 1
    T_univariate[1:d, (d + 1)] <- 1
    P1inf_univariate[1:d, 1:d] <- diag(1, d)
  }
  if (ar_length > 0)
    T_univariate[(d + 1):(d + ar_length), d + 1] <- ar
  if (m1 > (d + 1))
    T_univariate[(d + 1):(m1 - 1), (d + 2):m1] <- diag(1, max(ar_length, ma_length +
        1) - 1)
  R_univariate[d + 1, 1] <- 1
  if (ma_length > 0)
    R_univariate[(d + 2):(d + 1 + ma_length)] <- ma
  m <- p * m1
  Z <- matrix(0, p, m)
  T <- P1 <- P1inf <- matrix(0, m, m)
  R <- matrix(0, m, p)
  for (i in 1:p) {
    Z[i, ((i - 1) * m1 + 1):(i * m1)] <- Z_univariate
    T[((i - 1) * m1 + 1):(i * m1), ((i - 1) * m1 + 1):(i * m1)] <- T_univariate
    R[((i - 1) * m1 + 1):(i * m1), i] <- R_univariate
    P1inf[((i - 1) * m1 + 1):(i * m1), ((i - 1) * m1 + 1):(i * m1)] <- P1inf_univariate
  }
  if (stationary) {
    nd <- which(diag(P1inf) == 0)
    mnd <- length(nd)
    temp <- try(solve(a = diag(mnd^2) - matrix(kronecker(T[nd, nd], T[nd,
      nd]), mnd^2, mnd^2), b = c(R[nd, , drop = FALSE] %*%
          Q %*% t(R[nd, , drop = FALSE]))), TRUE)
    if (class(temp) == "try-error") {
      stop("ARIMA part is numerically too close to non-stationarity.")
    } else P1[nd, nd] <- temp
  } else diag(P1inf) <- 1
  state_names <- paste0(rep(paste0("arima", 1:m1), p), rep(ynames, each = m1))
  list(index = index, m = m, k = k, Z = Z, T = T, R = R, Q = Q, a1 = matrix(0,
    m, 1), P1 = P1, P1inf = P1inf, tvq = 0, tvr = 0, tvz = 0, state_names = state_names)
}
