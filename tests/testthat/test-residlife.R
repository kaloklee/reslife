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
test_that("gamma reduces to exponential when shape = 1", {
  rate_val <- 0.25
  scale_val <- 1 / rate_val
  
  # For shape = 1, MRL must equal scale (1/rate) across all t0 (memoryless)
  for (t0 in c(0, 2, 8)) {
    res_mrl_scale <- residlife(
      values = t0,
      distribution = "gamma",
      parameters = c(shape = 1, scale = scale_val),
      type = "mean"
    )
    res_mrl_rate <- residlife(
      values = t0,
      distribution = "gamma",
      parameters = c(shape = 1, rate = rate_val),
      type = "mean"
    )
    expect_equal(res_mrl_scale, scale_val, tolerance = 1e-4)
    expect_equal(res_mrl_rate, scale_val, tolerance = 1e-4)
  }
})

test_that("gamma accepts flexible parameter ordering and unnamed vectors", {
  shape_val <- 2.0
  scale_val <- 3.0
  rate_val  <- 1 / scale_val
  
  res_named_scale <- residlife(values = 2, distribution = "gamma", parameters = c(shape = shape_val, scale = scale_val))
  res_named_rate  <- residlife(values = 2, distribution = "gamma", parameters = c(rate = rate_val, shape = shape_val))
  res_unnamed     <- residlife(values = 2, distribution = "gamma", parameters = c(shape_val, scale_val))
  
  expect_equal(res_named_rate, res_named_scale, tolerance = 1e-4)
  expect_equal(res_unnamed, res_named_scale, tolerance = 1e-4)
})

test_that("gamma MRL at t0 = 0 matches unconditional expectation a * lambda", {
  # Case 1: shape > 1
  shape_gt1 <- 2.5
  scale_val <- 4.0
  res_gt1 <- residlife(values = 0, distribution = "gamma", parameters = c(shape = shape_gt1, scale = scale_val), type = "mean")
  expect_equal(res_gt1, shape_gt1 * scale_val, tolerance = 1e-4)
  
  # Case 2: shape < 1
  shape_lt1 <- 0.5
  res_lt1 <- residlife(values = 0, distribution = "gamma", parameters = c(shape = shape_lt1, scale = scale_val), type = "mean")
  expect_equal(res_lt1, shape_lt1 * scale_val, tolerance = 1e-4)
})

test_that("gompertz accepts flexible parameter ordering and scale/rate conversion", {
  shape_val <- 0.05
  rate_val  <- 0.01
  scale_val <- 1 / rate_val
  
  res_rate  <- residlife(values = c(10, 20), distribution = "gompertz", parameters = c(shape = shape_val, rate = rate_val))
  res_scale <- residlife(values = c(10, 20), distribution = "gompertz", parameters = c(scale = scale_val, shape = shape_val))
  res_unnamed <- residlife(values = c(10, 20), distribution = "gompertz", parameters = c(shape_val, rate_val))
  
  expect_equal(res_scale, res_rate, tolerance = 1e-4)
  expect_equal(res_unnamed, res_rate, tolerance = 1e-4)
})

test_that("gompertz MRL matches numerical integration of survival function", {
  shape_val <- 0.04
  rate_val  <- 0.02
  t0 <- 15
  
  # Numerical integration of S(u)/S(t0) from t0 to Inf
  surv_gompertz <- function(u) {
    exp((rate_val / shape_val) * (1 - exp(shape_val * u)))
  }
  st0 <- surv_gompertz(t0)
  expected_mrl <- integrate(function(u) surv_gompertz(u) / st0, lower = t0, upper = Inf)$value
  
  res <- residlife(
    values = t0,
    distribution = "gompertz",
    parameters = c(shape = shape_val, rate = rate_val),
    type = "mean"
  )
  expect_equal(res, expected_mrl, tolerance = 1e-4)
})

test_that("gengamma.orig reduces to exponential when b = 1 and k = 1", {
  scale_val <- 4.0
  
  # When shape (b) = 1 and k = 1, Stacy's gengamma reduces to Exponential(rate = 1/scale)
  # MRL must equal scale across all evaluation points
  for (t0 in c(0, 1.5, 5)) {
    res_mrl <- residlife(
      values = t0,
      distribution = "gengamma.orig",
      parameters = c(shape = 1, scale = scale_val, k = 1),
      type = "mean"
    )
    expect_equal(res_mrl, scale_val, tolerance = 1e-4)
  }
})

test_that("gengamma.orig reduces to weibull when k = 1", {
  shape_val <- 1.8
  scale_val <- 3.5
  t_vals <- c(0, 1, 3)
  
  res_gengamma <- residlife(
    values = t_vals,
    distribution = "gengamma.orig",
    parameters = c(shape = shape_val, scale = scale_val, k = 1),
    type = "all"
  )
  
  res_weibull <- residlife(
    values = t_vals,
    distribution = "weibull",
    parameters = c(shape = shape_val, scale = scale_val),
    type = "all"
  )
  
  expect_equal(res_gengamma, res_weibull, tolerance = 1e-4)
})

test_that("gengamma.orig accepts flexible parameter ordering and unnamed vectors", {
  b_val <- 1.5
  a_val <- 2.5
  k_val <- 0.8
  
  res_std <- residlife(
    values = 2,
    distribution = "gengamma.orig",
    parameters = c(shape = b_val, scale = a_val, k = k_val)
  )
  res_reordered <- residlife(
    values = 2,
    distribution = "gengamma.orig",
    parameters = c(k = k_val, scale = a_val, shape = b_val)
  )
  res_unnamed <- residlife(
    values = 2,
    distribution = "gengamma.orig",
    parameters = c(b_val, a_val, k_val)
  )
  
  expect_equal(res_reordered, res_std, tolerance = 1e-4)
  expect_equal(res_unnamed, res_std, tolerance = 1e-4)
})

test_that("gengamma matches gengamma.orig via parameter transformation", {
  # Parameter equivalences:
  # k = 1 / Q^2  => Q = 1 / sqrt(k)
  # b = abs(Q) / sigma => sigma = Q / b
  # a = exp(mu - (log(Q^(-2)) * sigma) / abs(Q)) => mu = log(a) + (log(Q^(-2)) * sigma) / abs(Q)
  b_val <- 1.5
  a_val <- 3.0
  k_val <- 0.8
  
  Q_val <- 1 / sqrt(k_val)
  sigma_val <- Q_val / b_val
  mu_val <- log(a_val) + (log(Q_val^(-2)) * sigma_val) / abs(Q_val)
  
  t_vals <- c(0, 1.2, 4.0)
  
  res_orig <- residlife(
    values = t_vals,
    distribution = "gengamma.orig",
    parameters = c(shape = b_val, scale = a_val, k = k_val),
    type = "all"
  )
  
  res_prentice <- residlife(
    values = t_vals,
    distribution = "gengamma",
    parameters = c(mu = mu_val, sigma = sigma_val, Q = Q_val),
    type = "all"
  )
  
  expect_equal(res_prentice, res_orig, tolerance = 1e-4)
})

test_that("gengamma accepts flexible parameter ordering and unnamed vectors", {
  mu_val    <- 1.2
  sigma_val <- 0.6
  Q_val     <- 0.9
  
  res_std <- residlife(
    values = 1.5,
    distribution = "gengamma",
    parameters = c(mu = mu_val, sigma = sigma_val, Q = Q_val)
  )
  res_reordered <- residlife(
    values = 1.5,
    distribution = "gengamma",
    parameters = c(Q = Q_val, mu = mu_val, sigma = sigma_val)
  )
  res_unnamed <- residlife(
    values = 1.5,
    distribution = "gengamma",
    parameters = c(mu_val, sigma_val, Q_val)
  )
  
  expect_equal(res_reordered, res_std, tolerance = 1e-4)
  expect_equal(res_unnamed, res_std, tolerance = 1e-4)
})

test_that("gengamma reduces to weibull when Q = 1", {
  # In Prentice form, Q = 1 corresponds to Weibull with:
  # shape = 1 / sigma
  # scale = exp(mu)
  mu_val <- 1.5
  sigma_val <- 0.8
  
  shape_weibull <- 1 / sigma_val
  scale_weibull <- exp(mu_val)
  t_vals <- c(0.5, 2.0)
  
  res_gengamma <- residlife(
    values = t_vals,
    distribution = "gengamma",
    parameters = c(mu = mu_val, sigma = sigma_val, Q = 1),
    type = "all"
  )
  
  res_weibull <- residlife(
    values = t_vals,
    distribution = "weibull",
    parameters = c(shape = shape_weibull, scale = scale_weibull),
    type = "all"
  )
  
  expect_equal(res_gengamma, res_weibull, tolerance = 1e-4)
})