# Default initial values for linear predictor theta
initTheta <- function(y, u, distribution) {
  ind<-distribution == "poisson"
  if (any(ind)) {
    x <- y[, ind]/u[, ind]
    x[x < 0.1 | is.na(x)] <- 0.1
    y[, ind] <- log(x)
  }
  ind<-distribution == "binomial"
  if (any(ind)) {
    y[, ind] <- qlogis((ifelse(is.na(y[, ind]), 0.5, y[, ind]) + 0.5)/(u[, ind] + 1))
  }
  ind<-distribution == "gamma"
  if (any(ind)) {
    x <- y[, ind]
    x[is.na(x) | x < 1] <- 1
    y[, ind] <- log(x)
  }
  ind <- distribution == "negative binomial"
  if (any(ind)) {
    x <- y[, ind]
    x[is.na(x) | x < 1/6] <- 1/6
    y[, ind] <- log(x)
  }
  y
} 
