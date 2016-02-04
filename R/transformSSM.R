#' Transform Multivariate State Space Model for Sequential Processing
#' 
#' \code{transformSSM} transforms the general multivariate Gaussian state space model
#' to form suitable for sequential processing.
#' 
#' @details As all the functions in KFAS use univariate approach i.e. sequential processing, 
#'   the covariance matrix \eqn{H_t}{H[t]} of the observation equation needs to be 
#'   either diagonal or zero matrix. Function \code{transformSSM} performs either 
#'   the LDL decomposition of \eqn{H_t}{H[t]}, or augments the state vector with
#'   the disturbances of the observation equation.
#'   
#'   In case of a LDL decomposition, the new \eqn{H_t}{H[t]} contains the diagonal part of the 
#'   decomposition, whereas observations \eqn{y_t}{y[t]} and system matrices \eqn{Z_t}{Z[t]} are 
#'   multiplied with the inverse of \eqn{L_t}{L[t]}. Note that although the state estimates and 
#'   their error covariances obtained by Kalman filtering and smoothing are identical with those 
#'   obtained from ordinary multivariate filtering, the one-step-ahead errors 
#'   \eqn{v_t}{v[t]} and their variances \eqn{F_t}{F[t]} do differ. The typical 
#'   multivariate versions can be obtained from output of \code{\link{KFS}}
#'   using \code{\link{mvInnovations}} function.
#'   
#'   
#' @export
#' @param object State space model object from function \code{\link{SSModel}}.
#' @param type Option \code{"ldl"} performs LDL decomposition for covariance matrix \eqn{H_t}{H[t]},
#'   and multiplies the observation equation with the \eqn{L_t^{-1}}{L[t]^-1}, so \eqn{\epsilon_t^* 
#'   \sim N(0,D_t)}{\epsilon[t]* ~ N(0,D[t])}. Option \code{"augment"} adds 
#'   \eqn{\epsilon_t}{\epsilon[t]} to the state vector, so \eqn{Q_t}{Q[t]} becomes block diagonal 
#'   with blocks \eqn{Q_t}{Q[t]} and \eqn{H_t}{H[t]}.
#' @return \item{model}{Transformed model.}
transformSSM <- function(object, type = c("ldl", "augment")) {
  
  if (any(object$distribution != "gaussian")) 
    stop("Nothing to transform as matrix H is not defined for non-gaussian model.")
  is.SSModel(object, return.logical = FALSE)
  type <- match.arg(type, choices = c("ldl", "augment"))
  p <- attr(object, "p")
  n <- attr(object, "n")
  m <- attr(object, "m")
  r <- attr(object, "k")
  tv <- attr(object,"tv")
  tvh <- tv[2]    
  if (type == "ldl") {
    if (p > 1) { #do nothing for univariate series
      yt <- t(object$y)
      ymiss <- is.na(yt)
      tv[1] <- max(tv[1], tv[2])       
      if (sum(ymiss) > 0) {
        # find unique combinations of missing observations 
        # even though H (and Z) becomes time varying in case of partial missingness, 
        # no need to compute Cholesky for all t
        positions <- unique(ymiss, MARGIN = 2)
        nh <- dim(positions)[2]
        tv[1:2] <- 1
        Z <- array(object$Z, dim = c(p, m, n))
      } else {
        Z <- array(object$Z, dim = c(p, m, (n - 1) * tv[1] + 1))
        positions <- rep(FALSE, p)
        nh <- 1
      }
      positions <- as.matrix(positions)
      H <- array(object$H, c(p, p, n))
      if (tvh) { 
        # if H was already time varying, compute cholesky decompositions for all t
        nh <- n
        hchol <- 1:n
        uniqs <- 1:n
      } else {
        hchol <- rep(0, n)
        uniqs <- numeric(nh)
        for (i in 1:nh) {
          nhn <- which(colSums(ymiss == positions[, i]) == p)
          hchol[nhn] <- i
          uniqs[i] <- nhn[1]
        }
      }
      ichols <- H[, , uniqs, drop = FALSE]  
      ydims <- as.integer(colSums(!ymiss))
      yobs <- array(1:p, c(p, n))
      if (sum(ymiss) > 0) #positions of missing observations 
        for (i in 1:n) {
          if (ydims[i] != p && ydims[i] != 0) 
            yobs[1:ydims[i], i] <- yobs[!ymiss[, i], i]
          if (ydims[i] < p) 
            yobs[(ydims[i] + 1):p, i] <- NA
        }
      unidim <- ydims[uniqs]
      hobs <- yobs[, uniqs, drop = FALSE]
      storage.mode(yobs) <- storage.mode(hobs) <- storage.mode(hchol) <- "integer"
      # compute the transformations for y, Z and H
      out <- .Fortran(fldlssm, NAOK = TRUE, yt = yt, ydims = ydims, yobs = yobs, 
        tv = as.integer(tv), Zt = Z, p = p, m = m, 
        n = n, ichols = ichols, nh = as.integer(nh), hchol = hchol, 
        unidim = as.integer(unidim), info = as.integer(0), hobs = hobs, 
        tol = max(100, max(abs(apply(object$H, 3, diag)))) * .Machine$double.eps)
      if(out$info!=0){
        stop(switch(as.character(out$info),
          "1" = "LDL decomposition of H failed.",
          "2" = "Computing the inverse of L failed."                     
        ))  
      }        
      H <- array(0, c(p, p, ((n - 1) * tv[2] + 1)))
      for (i in 1:((n - 1) * tv[2] + 1)) #construct diagonal H
        diag(H[, , i]) <- diag(out$ichols[, , out$hchol[i]])
      attry <- attributes(object$y)
      object$y <- t(out$yt)
      attributes(object$y) <- attry
      object$Z <- out$Z
      object$H <- H
      attr(object, "tv") <- as.integer(tv)
    }
  } else { #augmentation, add new states for epsilon disturbances and set H=0
    T <- array(object$T, dim = c(m, m, (n - 1) * tv[3] + 1))
    R <- array(object$R, dim = c(m, r, (n - 1) * tv[4] + 1))
    Q <- array(object$Q, dim = c(r, r, (n - 1) * tv[5] + 1))
    H <- array(object$H, c(p, p, (n - 1) * tv[2] + 1))
    Z <- array(object$Z, dim = c(p, m, (n - 1) * tv[1] + 1))
    r2 <- r + p
    m2 <- m + p
    tv[5] <- max(tv[c(2, 5)])
    Qt2 <- array(0, c(r2, r2, 1 + (n - 1) * tv[5]))
    Qt2[1:r, 1:r, ] <- Q
    if (tv[2]) {
      Qt2[(r + 1):r2, (r + 1):r2, -n] <- H[, , -1]
    } else Qt2[(r + 1):r2, (r + 1):r2, ] <- H
    Zt2 <- array(0, c(p, m2, (n - 1) * tv[1] + 1))
    Zt2[1:p, 1:m, ] <- Z
    Zt2[1:p, (m + 1):m2, ] <- diag(p)
    Tt2 <- array(0, c(m2, m2, (n - 1) * tv[3] + 1))
    Tt2[1:m, 1:m, ] <- T
    Rt2 <- array(0, c(m2, r + p, (n - 1) * tv[4] + 1))
    Rt2[1:m, 1:r, ] <- R
    Rt2[(m + 1):m2, (r + 1):r2, ] <- diag(p)
    P12 <- P1inf2 <- matrix(0, m2, m2)
    P12[1:m, 1:m] <- object$P1
    P1inf2[1:m, 1:m] <- object$P1inf
    P12[(m + 1):m2, (m + 1):m2] <- H[, , 1]
    a12 <- matrix(0, m2, 1)
    a12[1:m, ] <- object$a1
    object$Z <- Zt2
    object$H <- array(0, c(p, p, 1))
    object$T <- Tt2
    object$R <- Rt2
    object$Q <- Qt2        
    if (p == 1) {
      rownames(a12) <- c(rownames(object$a1), "eps")
    } else {
      rownames(a12) <- c(rownames(object$a1),paste0(rep("eps.", p), 1:p))
    }
    object$a1 <- a12
    object$P1 <- P12
    object$P1inf <- P1inf2
    attr(object, "m") <- m2
    attr(object, "k") <- r2
    attr(object, "tv") <- as.integer(tv)
    
  }   
  invisible(object)
} 
