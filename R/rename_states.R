#' Rename the States of SSModel Object
#'
#' A simple function for renaming the states of \code{\link{SSModel}} object. 
#' Note that since KFAS version 1.2.3 the auxiliary functions such as 
#' \code{\link{SSMtrend}} have argument \code{state_names} which can be used to 
#' overwrite the default state names when building the model with \code{\link{SSModel}}.
#' 
#' @param model Object of class SSModel
#' @param state_names Character vector giving new names for the states.
#' @return Original model with dimnames corresponding to states renamed.
#' @export
#' @examples 
#' custom_model <- SSModel(1:10 ~ -1 + 
#' SSMcustom(Z = 1, T = 1, R = 1, Q = 1, P1inf = 1), H = 1) 
#' custom_model <- rename_states(custom_model, "level")
#' ll_model <- SSModel(1:10 ~ SSMtrend(1, Q = 1), H = 1)
#' test_these <- c("y", "Z", "H", "T", "R", "Q", "a1", "P1", "P1inf")
#' identical(custom_model[test_these], ll_model[test_these])
#' 
rename_states <- function(model, state_names) {
  rownames(model$a1) <- rownames(model$T) <- colnames(model$T) <- colnames(model$Z) <-
    rownames(model$R) <- rownames(model$P1) <- colnames(model$P1) <- rownames(model$P1inf) <-
    colnames(model$P1inf) <- state_names
  model
}