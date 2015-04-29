context("Miscellancelous")

test_that("hatvalues works properly",{
model <- SSModel(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
out <- KFS(model, filtering = "state", smoothing = "none")
# estimate sigma2
model["H"] <- mean(c(out$v[1:out$d][out$Finf==0]^2/out$F[1:out$d][out$Finf==0],
                      out$v[-(1:out$d)]^2/out$F[-(1:out$d)]))
expect_equal(c(hatvalues(KFS(model))), 
             hatvalues(lm(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)),
             check.attributes = FALSE, tol=1e-6)
})

test_that("ldl works properly",{
  x <- matrix(c(1:5, (1:5)^2), 5, 2)
  x <- cbind(x, x[, 1] + 3*x[, 2])
  m <- crossprod(x)
  l <- ldl(m)
  d<-diag(diag(x))
  diag(l)<-1
  expect_equal(l%*%d%*%t(l), m)
})

