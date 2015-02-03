mvInnovations <- function(x){
  if(any(x$model$distribution!="gaussian"))
    stop("Function is only compatible with fully Gaussian models.")
  if(is.null(x[["a", exact = TRUE]]))
    stop("Function needs filtered estimates of states and their covariances.")
  out<-.Fortran(fmvfilter, NAOK = TRUE, attr(x$model, "tv")[1], 
                x$model$Z, attr(x$model, "p"), attr(x$model, "m"), attr(x$model, "n"),
                x$d, x$a[1:attr(x$model, "n"),], x$P[,,1:attr(x$model, "n")], 
                x$Pinf[,,1:x$d], v=x$model$y, 
                F=array(x$model$H,c(attr(x$model, "p"),attr(x$model, "p"),attr(x$model, "n"))),
                Finf=array(0,c(attr(x$model, "p"),attr(x$model, "p"),x$d)))
  out[c("v","F","Finf")]
}