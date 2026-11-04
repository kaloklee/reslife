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