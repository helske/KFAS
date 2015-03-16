test_that("predict.SSModel works",{
  tol <- 1e-10
  
  expect_that(model1 <- SSModel(VanKilled ~ SSMtrend(1, 0.001) +
                                  SSMseasonal(period = 12, sea.type = "dummy", Q = 0),
                                data = Seatbelts, distribution = "poisson"),not(gives_warning()))
  expect_that(model2 <- SSModel(VanKilled ~ law + SSMtrend(1, 0.001) +
                                  SSMseasonal(period = 12, sea.type = "dummy", Q = 0),
                                data = Seatbelts, distribution = "poisson"),not(gives_warning()))
  set.seed(1)
  expect_that(predict(model1, n.ahead = 24, states = "trend", interval = "prediction", nsim=5), 
              not(gives_warning()))
  set.seed(1)
  expect_that(predict(model2, states = c("trend", "regression"),  interval = "prediction", nsim=5,
                      newdata = SSModel(ts(rep(NA,24),frequency=12) ~ law + SSMtrend(1, 0.001) +
                                          SSMseasonal(period = 12, sea.type = "dummy", Q = 0),
                                        data = data.frame(law=rep(1,24)), distribution = "poisson")),
              not(gives_warning()))
  
  set.seed(1)
  expect_that(predict(model2, states = c("trend", "regression"), nsim=5,
                      newdata = SSModel(ts(rep(NA,24),frequency=12) ~ law + SSMtrend(1, 0.001) +
                                          SSMseasonal(period = 12, sea.type = "dummy", Q = 0),
                                        data = data.frame(law=rep(1,24)), distribution = "poisson")),
              not(gives_warning()))
  expect_that(predict(model2, states = c("trend", "regression"),  interval = "confidence", nsim=5,
                      newdata = SSModel(ts(rep(NA,24),frequency=12) ~ law + SSMtrend(1, 0.001) +
                                          SSMseasonal(period = 12, sea.type = "dummy", Q = 0),
                                        data = data.frame(law=rep(1,24)), distribution = "poisson")),
              not(gives_warning()))
  expect_that(predict(model2, states = c("trend", "regression"),  interval = "confidence",
                      newdata = SSModel(ts(rep(NA,24),frequency=12) ~ law + SSMtrend(1, 0.001) +
                                          SSMseasonal(period = 12, sea.type = "dummy", Q = 0),
                                        data = data.frame(law=rep(1,24)), distribution = "poisson")),
              not(gives_warning()))
})
