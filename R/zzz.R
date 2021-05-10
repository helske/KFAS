.onAttach <- function(libname, pkgname) {
  note <- "Please cite KFAS in publications by using: \n
  Jouni Helske (2017). KFAS: Exponential Family State Space Models in R. Journal of Statistical Software, 78(10), 1-39. doi:10.18637/jss.v078.i10."
  packageStartupMessage(note)
}
