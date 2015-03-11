test_that("importanceSSM works",{
  
data(sexratio)
model <- SSModel(Male ~ SSMtrend(1, Q = 1.107652e-06), u = sexratio[,"Total"], data = sexratio,
                 distribution = "binomial")
set.seed(1)
expect_that(imp <- importanceSSM(model, nsim = 25, antithetics = TRUE),not(gives_warning()))
expect_equal(sum(imp$w), 100.0004, tol=1e-6)
expect_equal(var(imp$s[1,1,]), 8.566762e-06, tol=1e-6)
expect_equal(var(imp$s[261,1,]), 8.223318e-06, tol=1e-6)

expect_that(imp <- importanceSSM(model, nsim = 25, filtered = TRUE, antithetics = TRUE),not(gives_warning()))
expect_equal(sum(imp$w), 26100.06, tol=1e-6)
expect_identical(unique(imp$s[1,1,]), 0)
expect_equal(var(imp$s[261,1,]), 1.071424e-05, tol=1e-7)
})
