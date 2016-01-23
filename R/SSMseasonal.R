#' @rdname SSModel
#' @export
SSMseasonal <- function(period, Q, sea.type = c("dummy", "trigonometric"),
  type, index, a1, P1, P1inf, n = 1, ynames) {
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
  sea.type <- match.arg(arg = sea.type)
  if (!(length(period) == 1 & period > 1))
    stop("Period of the seasonal component must be larger than 1.")
  period <- floor(period)
  m1 <- period - 1
  Z_univariate <- matrix(0, 1, m1)
  T_univariate <- matrix(0, m1, m1)
  if (sea.type == "dummy") {
    Z_univariate[1, 1] <- 1
    state_names <- paste0(rep(paste0("sea_dummy", 1:(period - 1)), each = 1),
      rep(ynames, each = period - 1))
    T_univariate[1, ] <- -1
    T_univariate[cbind(2:m1, 1:(m1 - 1))] <- 1
  } else {
    Z_univariate[1, ] <- rep(c(1, 0), length.out = period - 1)
    state_names <- paste0(rep(c("sea_trig", "sea_trig*"), each = 1, length.out = (period -
        1)), rep(1:floor(period/2), each = 2, length.out = (period - 1)), rep(ynames,
          each = period - 1))
    lambda <- 2 * pi * 1:floor((period - 1)/2)/period
    T_univariate[cbind(1:m1, 1:m1)] <- rep(c(cos(lambda), -1), each = 2, length = m1)
    T_univariate[which((col(T_univariate) - row(T_univariate)) == 1)[seq(from = 1,
      by = 2, length = length(lambda))]] <- sin(lambda)
    T_univariate[which((col(T_univariate) - row(T_univariate)) == -1)[seq(from = 1,
      by = 2, length = length(lambda))]] <- -sin(lambda)
  }
  m <- ((p - 1) * (type != 2) + 1) * (period - 1)
  T <- matrix(0, m, m)
  Z <- matrix(0, p, m)
  if (type != 2) {
    for (i in 1:p) {
      Z[i, ((i - 1) * m1 + 1):(i * m1)] <- Z_univariate
      T[((i - 1) * m1 + 1):(i * m1), ((i - 1) * m1 + 1):(i * m1)] <- T_univariate
    }
  } else {
    Z <- matrix(Z_univariate, nrow = p, ncol = m, byrow = TRUE)
    T <- T_univariate
  }
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
        stop("Misspecified Q, argument Q must be (p x p) matrix, (p x p x 1), or (p x p x n) array where p is the number of time series.")
    } else {
      if (length(Q) != 1 && (is.null(dim(Q)) || any(dim(Q)[1:2] != 1) || !(max(dim(Q)[3], 1, na.rm = TRUE) %in% c(1, n))))
        stop("Misspecified Q, argument Q must be a scalar, (1 x 1) matrix, or (1 x 1 x 1)/(1 x 1 x n) array.")
    }
    tvq <- max(dim(Q)[3] == n, 0, na.rm = TRUE)

    if (sea.type == "dummy") {
      if (type == 1) {
        Qm <- array(Q, c(p, p, tvq * (n - 1) + 1))
        k <- p
        R <- matrix(0, m, k)
        R[cbind(seq(1, m, period - 1), 1:k)] <- 1
      } else {
        Qm <- array(Q, c(1, 1, tvq * (n - 1) + 1))
        k <- 1
        R <- diag(m)[, 1, drop = FALSE]
      }
    } else {
      Qm <- array(0, c(m, m, tvq * (n - 1) + 1))
      if (type == 1) {
        if (tvq == 1) {
          for (i in 1:(tvq * (n - 1) + 1)){
            Qm[cbind(rep(1:(p * (period - 1)), p),
              rep(1:(period - 1), p^2) + rep(0:(p - 1) * (period - 1),
                each = p * (period - 1)), i)] <- rep(Q[, , i], each = (period - 1))
          }
        } else {
          Qm[cbind(rep(1:(p * (period - 1)), p),
            rep(1:(period - 1), p^2) + rep(0:(p - 1) * (period - 1),
              each = p * (period - 1)), 1)] <- rep(Q, each = (period - 1))
        }
      } else {
        if (tvq == 1) {
          for (i in 1:(tvq * (n - 1) + 1))
            Qm[cbind(1:(period - 1), 1:(period - 1), i)] <- Q[, , i]
        } else Qm[cbind(1:(period - 1), 1:(period - 1), 1)] <- Q
      }
      k <- m
      R <- diag(k)
    }
  }
  list(index = index, m = m, k = k, Z = Z, T = T, R = R, Q = Qm, a1 = a1, P1 = P1,
    P1inf = P1inf, tvq = tvq, tvr = 0, tvz = 0, state_names = state_names, period = period,
    sea.type = sea.type)
}
