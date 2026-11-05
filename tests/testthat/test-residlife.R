test_that("exponential distribution satisfies the memoryless property for MRL", {
  rate <- 0.2
  expected_mrl <- 1 / rate  # 5.0
  
  for (t0 in c(0, 5, 10, 25)) {
    res <- residlife(
      values = t0,
      distribution = "exponential",
      parameters = c(rate = rate),
      type = "mean"
    )
    expect_equal(res, expected_mrl, tolerance = 1e-4)
  }
})

test_that("exponential distribution median residual life is invariant to elapsed time", {
  rate <- 0.2
  expected_median <- log(2) / rate  # ~3.4657
  
  for (t0 in c(0, 5, 15)) {
    res <- residlife(
      values = t0,
      distribution = "exponential",
      parameters = c(rate = rate),
      p = 0.5,
      type = "percentile"
    )
    expect_equal(res, expected_median, tolerance = 1e-4)
  }
})

test_that("percentiles are monotonically ordered", {
  p25 <- residlife(values = 2, distribution = "exponential", parameters = c(rate = 0.5), p = 0.25, type = "percentile")
  p50 <- residlife(values = 2, distribution = "exponential", parameters = c(rate = 0.5), p = 0.50, type = "percentile")
  p75 <- residlife(values = 2, distribution = "exponential", parameters = c(rate = 0.5), p = 0.75, type = "percentile")
  
  expect_lt(p25, p50)
  expect_lt(p50, p75)
})

test_that("weibull reduces to exponential when shape = 1", {
  scale_val <- 5
  
  # For shape = 1, MRL must equal scale across all t0 (memoryless)
  for (t0 in c(0, 2, 10)) {
    res_mrl <- residlife(
      values = t0,
      distribution = "weibull",
      parameters = c(shape = 1, scale = scale_val),
      type = "mean"
    )
    expect_equal(res_mrl, scale_val, tolerance = 1e-4)
  }
})

test_that("weibull accepts flexible parameter ordering and unnamed vectors", {
  scale_val <- 4
  shape_val <- 2
  
  res_named_std <- residlife(values = 1, distribution = "weibull", parameters = c(shape = shape_val, scale = scale_val))
  res_named_rev <- residlife(values = 1, distribution = "weibull", parameters = c(scale = scale_val, shape = shape_val))
  res_unnamed   <- residlife(values = 1, distribution = "weibull", parameters = c(shape_val, scale_val))
  
  expect_equal(res_named_rev, res_named_std)
  expect_equal(res_unnamed, res_named_std)
})

test_that("weibull MRL at t0 = 0 matches theoretical expectation", {
  shape_val <- 1.5
  scale_val <- 3.0
  expected_mean <- scale_val * gamma(1 + 1 / shape_val)
  
  res <- residlife(
    values = 0,
    distribution = "weibull",
    parameters = c(shape = shape_val, scale = scale_val),
    type = "mean"
  )
  expect_equal(res, expected_mean, tolerance = 1e-4)
})
