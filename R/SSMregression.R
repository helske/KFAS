#' @rdname SSModel
#' @export
SSMregression <-  function(rformula, data, type, Q, index, R, a1, P1,
  P1inf, n = 1, ynames,remove.intercept=TRUE) {
  if (missing(index))
    index <- 1
  p <- length(index)
  if (!missing(ynames) && !is.null(ynames)) {
    ynames <- paste0(".", ynames)
  } else ynames <- ""
  if (missing(data)) {
    data <- environment(rformula)
  }
  if (missing(type)) {
    type <- 1L
  } else {
    type <- pmatch(x = type, table = c("distinct", "common"))
    if (is.na(type))
      stop("type must be 'distinct' or 'common'.")
  }
  old_option <- getOption("na.action")
  options(na.action = "na.pass")
  # case 1, input is formula
  if (inherits(rformula, "formula")) {
    if(remove.intercept){
      rformula <- update.formula(rformula, ~. + 1)
      X <- model.matrix(rformula, data = data)
      X <- X[, -(colnames(X) == "(Intercept)"), drop = FALSE]
    } else {
      if (length(attr(terms(rformula, data = data), "term.labels")) == 0 &&
          attr(terms(rformula,  data = data), "intercept") == 1) {
        X <- matrix(1, nrow = n, ncol = 1)
        colnames(X) <- "(Intercept)"
      } else  X <- model.matrix(rformula, data = data)

    }
    Xnames <- colnames(X)
    dims <- dim(X)
    if (missing(n)) {
      n <- dims[1]
    } else {
      if (n != dims[1])
        stop("Length of the series and covariates differ.")
    }
    if (any(is.na(X)))
      warning("Missing covariate values.")
    m <- dims[2] + dims[2] * ((p - 1) * (type != 2))
    Z <- array(0, c(p, m, n))
    if (type == 2) {
      for (i in 1:m) Z[, i, ] <- rep(X[, i], each = p)
    } else {
      for (i in 1:p) Z[i, ((i - 1) * dims[2] + 1):(i * dims[2]), ] <- t(X)
    }
    state_names <- paste0(rep(Xnames, times = (p - 1) * (type == 1) + 1), rep(ynames,
      each = dims[2]))
  } else {
    if (length(rformula) != p)
      stop("Length of the formula list is not equal to the number of series.")
    X <- vector("list", length = p)
    if (is.list(data) && !is.data.frame(data)) {
      if (length(data) != p)
        stop("Length of the data list is not equal to the number of series.")
      for (i in 1:p) {
        if(remove.intercept){
          rformula[[i]] <- update.formula(rformula[[i]], ~. + 1)
          X[[i]] <- model.matrix(rformula[[i]], data = data[[i]])
          X[[i]] <- X[[i]][, -(colnames(X[[i]]) == "(Intercept)"), drop = FALSE]
        }else {
          if (length(attr(terms(rformula[[i]], data = data[[i]]), "term.labels")) == 0 &&
              attr(terms(rformula[[i]],  data = data[[i]]), "intercept") == 1) {
            X[[i]] <- matrix(1, nrow = n, ncol = 1)
            colnames(X[[i]]) <- "(Intercept)"
          } else X[[i]] <- model.matrix(rformula[[i]], data = data[[i]])
        }
      }
    } else {
      for (i in 1:p) {
        if(remove.intercept){
          rformula[[i]] <- update.formula(rformula[[i]], ~. + 1)
          X[[i]] <- model.matrix(rformula[[i]], data = data)
          X[[i]] <- X[[i]][, -(colnames(X[[i]]) == "(Intercept)"), drop = FALSE]
        }else {
          if (length(attr(terms(rformula[[i]], data = data), "term.labels")) == 0 &&
              attr(terms(rformula[[i]],  data = data), "intercept") == 1) {
            X[[i]] <- matrix(1, nrow = n, ncol = 1)
            colnames(X[[i]]) <- "(Intercept)"
          } else X[[i]] <- model.matrix(rformula[[i]], data = data)
        }
      }
    }
    if (any(sapply(X, is.na)))
      warning("Missing values in X.")
    dims <- sapply(X, dim)
    if (missing(n)) {
      n <- dims[1, 1]
    } else {
      if (any(dims[1, ] != n))
        stop("Length of the series and covariates differ.")
    }
    if (type == 2 & length(unique(dims[2, ])) > 1)
      stop("Unequal number of covariates for different series.")
    if (type == 1) {
      m <- sum(dims[2, ])
      m_cumsum <- c(0, cumsum(dims[2, ]))
      Z <- array(0, dim = c(p, m, n))
      state_names <- NULL
      for (i in 1:p) {
        state_names <- c(state_names, paste0(colnames(X[[i]]), ynames[i]))
        Z[i, (m_cumsum[i] + 1):m_cumsum[i + 1], ] <- t(X[[i]])
      }
    } else {
      dims <- dims[, 1]
      m <- dims[2]
      Z <- array(0, c(p, m, n))
      for (j in 1:p) {
        for (i in 1:dims[2]) {
          Z[j, i, ] <- X[[j]][, i]
        }
      }
      state_names <- colnames(X[[1]])  #paste0(rep('beta',m),1:m)
    }
  }
  T <- diag(m)
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
    Q <- NULL
    tvq <- 0
  } else {
    if (length(Q) == 1)
      Q <- matrix(Q)
    if (!identical(dim(Q)[1], dim(Q)[2]) || isTRUE(dim(Q)[1] > m) || !(max(dim(Q)[3],
      1, na.rm = TRUE) %in% c(1, n)))
      stop(paste0("Misspecified Q, argument Q must be (k x k) matrix, (k x k x 1), or ",
        "(k x k x n) array where m is the number of disturbance terms."))
    k <- dim(Q)[1]
    tvq <- max(dim(Q)[3] == n, 0, na.rm = TRUE)
  }
  if (missing(R)) {
    tvr <- 0
    if (k > 0) {
      R <- diag(m)[, 1:k, drop = FALSE]
    } else R <- NULL
  } else {
    if (isTRUE(!(dim(R)[1] == m)) || isTRUE(dim(R)[2] != k) ||
        !(max(dim(R)[3], 1, na.rm = TRUE) %in% c(1, n)))
      stop(paste0("Misspecified R, argument R must be (m x k) matrix, (m x k x 1), or ",
        "(m x k x n) array where m is the number of states and k is the number of disturbance terms."))
    tvr <- max(dim(R)[3] == n, 0, na.rm = TRUE)
  }
  if (dim(unique(Z,MARGIN=3))[[3]] == 1) {
    Z <- Z[, , 1, drop = FALSE]
    tvz <- 0
  }
  else tvz <- 1
  options(na.action = old_option)
  list(index = index, m = m, k = k, Z = Z, T = T, R = R, Q = Q,
    a1 = a1, P1 = P1, P1inf = P1inf, tvq = tvq, tvr = tvr,
    tvz = tvz, state_names = state_names)
}
