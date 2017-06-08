# functions for standardized residuals

varianceFilter <- function(object) {
  
  vars <- object$model$y
  for (i in 1:length(object$model$distribution)){
    vars[, i] <- switch(object$model$distribution[i], 
      gaussian = object$model$u[,i], 
      poisson = object$m[, i], 
      binomial = object$m[, i] * (1 - object$m[, i]) / object$model$u[, i], 
      gamma = object$m[, i]^2 / object$model$u[,i], 
      `negative binomial` = object$m[, i] + object$m[, i]^2 / object$model$u[, i])
  }
  vars
}
varianceSmoother <- function(object) {
  
  vars <- object$model$y
  for (i in 1:length(object$model$distribution)){
    vars[, i] <- switch(object$model$distribution[i], 
      gaussian = object$model$u[,i], 
      poisson = object$muhat[, i], 
      binomial = object$muhat[, i] * (1 - object$muhat[, i]) / object$model$u[, i], 
      gamma = object$muhat[, i]^2 / object$model$u[,i], 
      `negative binomial` = object$muhat[, i] + object$muhat[, i]^2 / object$model$u[, i])
  }
  vars
}

recursive_standardized <- function(object, stype) {
  
  if(any(object$model$distribution !=  "gaussian") && !("m" %in% names(object)))
    stop("KFS object does not contain filtered means. ")
  if (all(object$model$distribution ==  "gaussian") && !("v" %in% names(object)))
    stop("KFS object does not contain prediction errors v. ")
  
  p<-attr(object$model,"p")
  
  if(all(object$model$distribution ==  "gaussian")){
    if(stype == "cholesky" || p == 1){
      res <- object$v / sqrt(t(object$F))
      if (object$d > 0) {
        res[1:(object$d - 1), ] <- NA
        res[object$d, 1:object$j] <- NA
      } 
    } else {
      tmp <- mvInnovations(object)[c("v", "F")]
      res <- tmp$v / sqrt(t(apply(tmp$F, 3, diag)))
      res[1:object$d, ] <- NA      
    }
  } else {
    d<-KFS(approxSSM(object$model),filtering="state",smoothing="none")$d      
    vars<-varianceFilter(object)
    res<-object$model$y  
    if (sum(bins <- object$model$distribution == "binomial") > 0)
      res[, bins] <- res[, bins]/object$model$u[, bins]
    res <- res-object[["m", exact = TRUE]]   
    
    if(stype == "cholesky" & p > 1){
      for(i in (d + 1):attr(object$model, "n") ){
        yobs <- !is.na(res[i, ])
        if(sum(yobs) > 0){
          L <- ldl(diag(vars[i, yobs], p) + object$P_mu[yobs, yobs, i])
          D <- diag(L)
          diag(L) <- 1
          res[i, yobs] <- (1 / sqrt(D) * 
              backsolve(L, diag(nrow(L)), upper.tri = FALSE)) %*% res[i, yobs]
        }        
      }
    } else {        
      for(i in (d + 1):attr(object$model, "n") ){
        yobs <- !is.na(res[i, ])
        if(sum(yobs) > 0){
          res[i, ] <- res[i, ] / sqrt(vars[i, ] + 
              object$P_mu[, , i][1 + 0:(p - 1) * (p + 1)])
        }
      }      
    }
    if(d > 0){
      res[1:d, ] <- NA                                    
    }
    
  }
  res
}

pearson_standardized <- function(object, stype) {
  
  n <- attr(object$model, "n")
  
  if(!("muhat" %in% names(object)))
    stop("KFS object needs to contain smoothed means. ")
  p <- attr(object$model,"p")
  
  if (all(object$model$distribution == "gaussian")) {
    
    tv <- dim(object$model$H)[3] > 1
    res <- object$model$y - object$muhat 
    
    if(stype == "cholesky" & p > 1){
      for(t in 1:n){
        yobs <- !is.na(res[t,])
        if(sum(yobs) > 0){
          L <- ldl(object$model$H[yobs, yobs, (t - 1) * tv + 1] - object$V_mu[yobs, yobs, t])
          D <- diag(L)
          diag(L) <- 1
          pos <- D > sqrt(.Machine$double.eps) * max(D, 0)
          res[t, yobs][pos] <- (1 / sqrt(D[pos]) * 
              backsolve(L[pos, pos], diag(sum(pos)), upper.tri = FALSE)) %*% res[t, yobs][pos]
          res[t, yobs][!pos] <- NA
        }
      }
    } else {
      for(t in 1:n){
        yobs <- !is.na(res[t,])
        if(sum(yobs) > 0){
          D <- sqrt((object$model$H[,,(t-1)*tv+1] - object$V_mu[,,t])[1 + 0:(p - 1) * (p + 1)])
          pos <- D > sqrt(.Machine$double.eps) * max(D, 0)
          res[t, pos] <- res[t, pos] / D[pos]
          res[t, !pos] <- NA
        }
      }
    }
    
  } else {
    res<-object$model$y
    vars<-varianceSmoother(object)     
    if (sum(bins <- object$model$distribution == "binomial") > 0) 
      res[, bins] <- res[, bins]/object$model$u[, bins]
    res <- res - object$muhat  
    if(stype == "cholesky" & p > 1){
      for(t in 1:n){
        yobs <- !is.na(res[t, ])
        if(sum(yobs) > 0){
          L <- ldl(diag(vars[t,yobs], p) - object$V_mu[yobs,yobs,t])
          D <- diag(L)
          diag(L) <- 1
          pos <- D > sqrt(.Machine$double.eps) * max(D, 0)
          res[t, yobs][pos] <- (1 / sqrt(D[pos]) * 
              backsolve(L[pos, pos], diag(sum(pos)), upper.tri = FALSE)) %*% res[t, yobs][pos]
          res[t, yobs][!pos] <- NA
        }
      }       
    } else {
      for(t in 1:n){
        yobs <- !is.na(res[t, ])
        if(sum(yobs) > 0){
          D <- sqrt(vars[t,] - object$V_mu[,,t][1 + 0:(p - 1) * (p + 1)])
          pos <- D > sqrt(.Machine$double.eps) * max(D, 0)
          res[t, pos] <- res[t, pos] / D[pos]
          res[t, !pos] <- NA
        }
      }
    }
  }
  res
}

state_standardized <- function(object, stype) {
  
  if (!("etahat" %in% names(object)))
    stop("KFS object needs to contain smoothed estimates of state disturbances eta.")
  
  k <- attr(object$model, "k")
  n <- attr(object$model, "n")
  eta <- object$etahat
  tvq <- attr(object$model, "tv")[5]
  
  if(stype == "cholesky" & k > 1){
    for (i in 1:(n - 1)){
      L <- ldl(object$model$Q[, , i * tvq + 1] - object$V_eta[, , i])
      D <- diag(L)
      diag(L) <- 1
      pos <- D > sqrt(.Machine$double.eps) * max(D, 0)
      eta[i, pos] <- (1 / sqrt(D[pos]) * 
          backsolve(L[pos, pos], diag(sum(pos)), upper.tri = FALSE)) %*% eta[i, pos]
      eta[i, !pos] <- NA
    }
  } else {
    for (i in 1:(n - 1)){        
      D <- sqrt((object$model$Q[, , i * tvq + 1] - object$V_eta[, , i])[1 + 0:(k - 1) * (k + 1)])
      pos <- D > sqrt(.Machine$double.eps) * max(D, 0)
      eta[i, pos] <- eta[i, pos] / D[pos]
      eta[i, !pos] <- NA
    }
  }
  eta[n, ] <- 0
  eta
  
}