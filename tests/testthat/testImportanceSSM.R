test_that("importanceSSM works",{
  tol <- 1e-10
  data(sexratio)
  model <- SSModel(Male ~ SSMtrend(1, Q = 0.1), u = sexratio[,"Total"], data = sexratio,
                   distribution = "binomial")
  set.seed(1)
  expect_warning(imp <- importanceSSM(model, nsim = 25, antithetics = TRUE),NA)
  expect_equal(sum(imp$w), 100.1432, tol=1e-6)
  expect_equal(var(imp$s[1,1,]), 0.0002033084, tol=1e-6)
  expect_equal(var(imp$s[261,1,]), 6.560228e-05, tol=1e-6)
  set.seed(1)
  expect_warning(out <- KFS(model, nsim = 25, smoothing=c("mean","state","signal")),NA)

  muhat<-colSums(imp$w/sum(imp$w)*plogis(t(imp$s[,1,])))
  expect_equal(muhat,c(out$mu),tol=tol)
  expect_equal(colSums(imp$w/sum(imp$w)*plogis(t(imp$s[,1,]))^2) - muhat^2,c(out$V_mu),tol=tol)

  set.seed(1)
  expect_warning(imp <- importanceSSM(model, nsim = 25, filtered = TRUE, antithetics = TRUE),NA)
  set.seed(1)
  expect_warning(out <- KFS(model, nsim = 25, filtering=c("mean","state","signal")),NA)
  expect_equal(sum(imp$w), 26122.95, tol=1e-6)

  mu<-sum(imp$w[1,]/sum(imp$w[1,])*plogis(imp$s[1,1,]))
  expect_equal(mu,out$m[1],tol=tol)
  expect_equal(sum(imp$w[1,]/sum(imp$w[1,])*plogis(imp$s[1,1,])^2) - mu^2,out$P_mu[1],tol=tol)
  mu<-sum(imp$w[261,]/sum(imp$w[261,])*plogis(imp$s[261,1,]))
  expect_equal(mu,out$m[261],tol=tol)
  expect_equal(sum(imp$w[261,]/sum(imp$w[261,])*plogis(imp$s[261,1,])^2) - mu^2,out$P_mu[261],tol=tol)
})
