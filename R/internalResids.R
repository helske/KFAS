# functions for standardized residuals

varianceFilter <- function(object) {
  
  vars <- object$model$y
  for (i in 1:length(object$model$distribution)){
    vars[, i] <- switch(object$model$distribution[i], 
      gaussian = object$model$u[,i], 
      poisson = object$m[, i], 
      binomial = object$m[, i] * (1 - object$m[, i])/object$model$u[, i], 
      gamma = object$m[, i]^2/object$model$u[,i], 
      `negative binomial` = object$m[, i] + object$m[, i]^2/object$model$u[, i])
  }
  vars
}
varianceSmoother <- function(object) {
  
  vars <- object$model$y
  for (i in 1:length(object$model$distribution)){
    vars[, i] <- switch(object$model$distribution[i], 
      gaussian = object$model$u[,i], 
      poisson = object$muhat[, i], 
      binomial = object$muhat[, i] * (1 - object$muhat[, i])/object$model$u[, i], 
      gamma = object$muhat[, i]^2/object$model$u[,i], 
      `negative binomial` = object$muhat[, i] + object$muhat[, i]^2/object$model$u[, i])
  }
  vars
}

recursive_standardized <- function(object,stype) {
  
  if(any(object$model$distribution !=  "gaussian") && is.null(object[["m", exact = TRUE]]))
    stop("KFS object does not contain filtered means. ")
  if (all(object$model$distribution ==  "gaussian") && is.null(object[["v", exact = TRUE]])) 
    stop("KFS object does not contain prediction errors v. ")
  
  p<-attr(object$model,"p")
  
  if(all(object$model$distribution ==  "gaussian")){
    if(stype=="cholesky" || attr(object$model,"p")==1){
      res <- object$v/sqrt(t(object$F))
      if (object$d > 0) {
        res[1:(object$d - 1), ] <- NA
        res[object$d, 1:object$j] <- NA
      } 
    } else {
      tmp <- mvInnovations(object)[c("v","F")]
      res <- tmp$v/sqrt(t(apply(tmp$F,3,diag)))
      res[1:object$d, ] <- NA      
    }
  } else {
    d<-KFS(approxSSM(object$model),filtering="state",smoothing="none")$d      
    vars<-varianceFilter(object)
    res<-object$model$y  
    if (sum(bins <- object$model$distribution == "binomial") > 0)
      res[, bins] <- res[, bins]/object$model$u[, bins]
    res <- res-object[["m", exact = TRUE]]   
    
    if(stype=="cholesky" || p==1){
      for(i in (d+1):attr(object$model, "n") ){
        yobs<-which(!is.na(res[i,]))
        if(length(yobs)>0){
          tmp <- chol(diag(vars[i,yobs],p)+object$P_mu[yobs,yobs,i])
          res[i,yobs] <- res[i,yobs]%*%solve(tmp)
        }        
      }
    } else {        
      for(i in (d+1):attr(object$model, "n") ){
        yobs<-which(!is.na(res[i,]))
        if(length(yobs)>0){
          res[i,yobs] <- res[i,yobs]/sqrt(vars[i,yobs]+diag(object$P_mu[yobs,yobs,i]))
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
  
  if(is.null(object$muhat))
    stop("KFS object needs to contain smoothed means. ")
  p<-attr(object$model,"p")
  
  if (all(object$model$distribution == "gaussian")) {
    
    tv<-dim(object$model$H)[3] > 1
    res<-object$model$y-object$muhat 
    
    if(stype=="cholesky" || p==1){
      for(t in 1:n){
        yobs<-which(!is.na(res[t,]))
        if(length(yobs)>0){
          tmp <- chol(object$model$H[yobs,yobs,(t-1)*tv+1]-object$V_mu[yobs,yobs,t])
          res[t,yobs] <- res[t,yobs]%*%solve(tmp)
        }
      }
    } else {
      for(t in 1:n){
        yobs<-which(!is.na(res[t,]))
        if(length(yobs)>0){
          res[t,yobs] <- res[t,yobs]/sqrt(diag(object$model$H[yobs,yobs,(t-1)*tv+1]-object$V_mu[yobs,yobs,t]))
        }
      }
    }
    
  } else {
    res<-object$model$y
    vars<-varianceSmoother(object)     
    if (sum(bins <- object$model$distribution == "binomial") > 0) 
      res[, bins] <- res[, bins]/object$model$u[, bins]
    res<-res-object$muhat  
    if(stype == "cholesky" || p==1){
      for(t in 1:n){
        yobs<-which(!is.na(res[t,]))
        if(length(yobs)>0){
          tmp <- chol(diag(vars[t,yobs],p)-object$V_mu[yobs,yobs,t])
          res[t,yobs] <- res[t,yobs]%*%solve(tmp)
        }
      }       
    } else {
      for(t in 1:n){
        yobs<-which(!is.na(res[t,]))
        if(length(yobs)>0){
          tmp <- chol(diag(vars[t,yobs],p)-object$V_mu[yobs,yobs,t])
          res[t,yobs] <- res[t,yobs]/sqrt(vars[t,yobs]-diag(object$V_mu[yobs,yobs,t]))
        }
      }
      
    }
    
  }
  res
}

state_standardized <- function(object, stype) {
  
  if (is.null(object$etahat))
    stop("KFS object needs to contain smoothed estimates of state disturbances eta.")
  
  k <- attr(object$model, "k")
  n <- attr(object$model, "n")
  eta <- object$etahat
  if(stype=="cholesky" || k==1){
    if (dim(object$model$Q)[3] == 1) {
      z <- which(object$model$Q[, , 1][1 + 0:(k - 1) * (k + 1)] > 0)   
      if(length(z)>0)
        for (i in 1:(n - 1)) {
          eta[i, z] <- eta[i, z] %*% solve(chol(object$model$Q[z, z, 1] - object$V_eta[z,z, i]))          
        }
      
    } else {
      for (i in 1:(n - 1)){
        z <- which(object$model$Q[, , i][1 + 0:(k - 1) * (k + 1)] > 0) 
        if(length(z)>0)
          eta[i, z] <- eta[i, z] %*% solve(chol(object$model$Q[z, z, i] - object$V_eta[z,z, i]))        
        
      }
    }    
  } else {
    if (dim(object$model$Q)[3] == 1) {
      z <- which(object$model$Q[, , 1][1 + 0:(k - 1) * (k + 1)] > 0)   
      if(length(z)>0)
        for (i in 1:(n - 1)) {
          eta[i, z] <- eta[i, z]/sqrt(diag(object$model$Q[z, z, 1] - object$V_eta[z,z, i]))
        }
      
    } else {
      for (i in 1:(n - 1)){        
        z <- which(object$model$Q[, , i][1 + 0:(k - 1) * (k + 1)] > 0)  
        if(length(z)>0)
          eta[i, z] <- eta[i, z]/sqrt(diag(object$model$Q[z, z, i] - object$V_eta[z,z, i]))        
        
      }
    }
    
  }
  eta[n, ] <- 0
  eta
  
}


deviance_standardized <- function(object) {
  
  if (all(object$model$distribution == "gaussian")) {
    w <- matrix(apply(object$model$H, 3, diag), attr(object$model, "n"), 
      attr(object$model, "p"), byrow = TRUE)
  } else {
    w <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
    for (i in 1:attr(object$model, "p")) 
      w[, i] <- 
      switch(object$model$distribution[i], 
        gaussian = object$model$u[, i], 
        poisson = 1, 
        binomial = 1, 
        gamma = 1/object$model$u[, i], 
        `negative binomial` = 1)
  }
  series <- object$model$y
  if (sum(bins <- object$model$distribution == "binomial") > 0) 
    series[, bins] <- series[, bins]/object$model$u[, bins]
  for (i in 1:attr(object$model, "p")) 
    series[, i] <- ifelse(series[, i] >  object$muhat[, i], 1, -1) * 
    sqrt(switch(object$model$distribution[i],  
      gaussian = (series[, i] - object$muhat[, i])^2, 
      poisson = 2 * (series[, i] * log(ifelse(series[, i] == 0, 1, 
        series[, i]/object$muhat[, i])) -   
          series[, i] + object$muhat[, i]),
      binomial = 2 * object$model$u[, i] *
        (series[, i] * log(ifelse(series[, i] == 0, 1, 
          series[, i]/object$muhat[, i])) + 
            (1 - series[, i]) * log(ifelse(series[, i] == 1 | object$muhat[, i] == 1, 
              1, (1 - series[, i])/(1 - object$muhat[, i])))), 
      gamma = -2 * (log(ifelse(object$model$y[, i] == 0, 1, 
        object$model$y[, i]/object$muhat[, i])) - 
          (object$model$y[, i] - object$muhat[, i])/object$muhat[,  i]), 
      `negative binomial` = 2 * (object$model$y[, i] * log(pmax(1, object$model$y[, i])/object$muhat[, i]) 
        - (object$model$y[, i] + object$model$u[, i]) * 
          log((object$model$y[, i] + 
              object$model$u[, i])/(object$muhat[, i] +
                  object$model$u[, i])))))
  series/sqrt(w * (1 - hatvalues(object)))
}
