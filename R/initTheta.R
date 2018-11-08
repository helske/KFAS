# Default initial values for linear predictor theta
#' @importFrom stats qlogis
initTheta <- function(y, u, distribution) {
  # time series division is much slower than matrix division
  y <- unclass(y)
  u <- unclass(u)
  if (any(ind <- distribution == "poisson")) {
    x <- y[, ind]/u[, ind]
    x[x < 0.1 | is.na(x)] <- 0.1
    y[, ind] <- log(x)
  }
  if (any(ind <- distribution == "binomial")) {
    y[, ind] <- qlogis((ifelse(is.na(y[, ind]), 0.5, y[, ind]) + 0.5)/(u[, ind] + 1))
  }
  if (any(ind <- distribution == "gamma")) {
    x <- y[, ind]
    x[is.na(x) | x < 1] <- 1
    y[, ind] <- log(x)
  }
  if (any(ind <- distribution == "negative binomial")) {
    x <- y[, ind]
    x[is.na(x) | x < 1/6] <- 1/6
    y[, ind] <- log(x)
  }
  y
}
