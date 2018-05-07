context("Comparison with lme4")

test_that("lmms works properly",{
  data("sleepstudy", package = "lme4")
  y_split <- split(sleepstudy["Reaction"], sleepstudy["Subject"])
  p <- length(y_split)
  y <- matrix(unlist(y_split), ncol = p,
    dimnames = list(NULL, paste("Subject", names(y_split))))
  
  P1 <- as.matrix(Matrix::.bdiag(replicate(p, matrix(NA, 2, 2), simplify = FALSE)))
  model_lmm <- SSModel(y ~ - 1 +
      SSMregression(~ Days, type = "common", remove.intercept = FALSE) +
      SSMregression(~ Days, remove.intercept = FALSE, P1 = P1),
    H = diag(NA, p), data = data.frame(Days = 0:9))
  
  update_lmm <- function(pars, model) {
    P1 <- diag(exp(pars[1:2]))
    P1[1, 2] <- pars[3]
    P1 <- crossprod(P1)
    model["P1", states = 3:38] <- 
      as.matrix(Matrix::.bdiag(replicate(p, P1, simplify = FALSE)))
    model["H"] <- diag(exp(pars[4]), p)
    model
  }
  
  fit_lmm <- fitSSM(model_lmm, c(1, 1, 1, 5), update_lmm, method = "BFGS")
  
  fm1 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
  
  expect_equal(lme4::fixef(fm1), KFS(fit_lmm$model)$alphahat[1,1:2])
  
  expect_equivalent(logLik(fit_lmm$model), unclass(logLik(fm1)))
  
  X <- model.matrix(fm1)
  expect_equal(logLik(fit_lmm$model, marginal = TRUE) - logLik(fit_lmm$model), 
    sum(log(diag(chol(t(X) %*% X)))))
})