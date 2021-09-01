data("pollution_data")

test_that("ev_trawl_fit_errors", {
  n <- 10
  depth <- 5
  data <- seq_len(n)

  testthat::expect_error(ev_trawl_fit(data, depth, "qwe"))
  testthat::expect_error(ev_trawl_fit(data, depth, "GMM", mode = "tmp"))
  testthat::expect_error(ev_trawl_fit(data, depth, "GMM", bounds = "tmp"))
})

test_that("ev_trawl_fit__does_not_smoke", {
  n <- 500
  depth <- 3
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]

  fit_params <- ev_trawl_fit(data, depth, "GMM")
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- ev_trawl_fit(data, depth, "PL")
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- ev_trawl_fit(data, depth, "GMM", mode = "full")
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- ev_trawl_fit(data, depth, "PL", mode = "full")
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))
})

test_that("sub_sample_fit__does_not_smoke", {
  n <- 1000
  sub_length <- 500
  depth <- 3
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]

  fit_params <- sub_sample_fit(
    data,
    sample_length = sub_length, depth = depth, method = "PL", trials = 1,
    mode = "full", type = "exp", bounds = "multiplier", parallel = F
  )
  fit_params <- fit_params$estimators[1, ]
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- sub_sample_fit(
    data,
    sample_length = sub_length, depth = depth, method = "GMM", trials = 1,
    mode = "full", type = "exp", bounds = "multiplier", parallel = F
  )
  fit_params <- fit_params$estimators[1, ]
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- sub_sample_fit(
    data,
    sample_length = sub_length, depth = depth, method = "PL", trials = 1,
    mode = "two-stage", type = "exp", bounds = "multiplier", parallel = F
  )
  fit_params <- fit_params$estimators[1, ]
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- sub_sample_fit(
    data,
    sample_length = sub_length, depth = depth, method = "GMM", trials = 1,
    mode = "two-stage", type = "exp", bounds = "multiplier", parallel = F
  )
  fit_params <- fit_params$estimators[1, ]
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))
})

test_that("sub_sample_fit__parallel__does_not_smoke", {
  testthat::skip_on_os("windows")
  n <- 1000
  sub_length <- 500
  depth <- 3
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]

  fit_params <- sub_sample_fit(
    data,
    sample_length = sub_length, depth = depth, method = "PL", trials = 1,
    mode = "full", type = "exp", bounds = "multiplier", parallel = T
  )
  fit_params <- fit_params$estimators[1, ]
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- sub_sample_fit(
    data,
    sample_length = sub_length, depth = depth, method = "GMM", trials = 1,
    mode = "full", type = "exp", bounds = "multiplier", parallel = T
  )
  fit_params <- fit_params$estimators[1, ]
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- sub_sample_fit(
    data,
    sample_length = sub_length, depth = depth, method = "PL", trials = 1,
    mode = "two-stage", type = "exp", bounds = "multiplier", parallel = T
  )
  fit_params <- fit_params$estimators[1, ]
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- sub_sample_fit(
    data,
    sample_length = sub_length, depth = depth, method = "GMM", trials = 1,
    mode = "two-stage", type = "exp", bounds = "multiplier", parallel = T
  )
  fit_params <- fit_params$estimators[1, ]
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))
})


test_that("fit__does_not_smoke", {
  n <- 1000
  depth <- 3
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]

  fit_params <- fit(
    data,
    depth = depth, method = "PL", mode = "full",
    type = "exp", bounds = "multiplier", parallel = F
  )
  fit_params <- fit_params$estimators
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- fit(
    data,
    depth = depth, method = "PL", mode = "full",
    type = "exp", bounds = "config", parallel = F
  )
  fit_params <- fit_params$estimators
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))
})


test_that("fit__parallel__does_not_smoke", {
  testthat::skip_on_os("windows")
  n <- 1000
  depth <- 3
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]

  fit_params <- fit(
    data,
    depth = depth, method = "PL", mode = "full",
    type = "exp", bounds = "multiplier", parallel = T
  )
  fit_params <- fit_params$estimators
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))

  fit_params <- fit(
    data,
    depth = depth, method = "PL", mode = "full",
    type = "exp", bounds = "config", parallel = T
  )
  fit_params <- fit_params$estimators
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:4] > 0))
})
