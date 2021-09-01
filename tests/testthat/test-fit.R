data("pollution_data")

test_that("ev_trawl_fit_errors", {
  n <- 2000
  depth <- 5
  data <- seq_len(n)

  testthat::expect_error(ev_trawl_fit(data, depth, "qwe"))
  testthat::expect_error(ev_trawl_fit(data, depth, "GMM", mode = "tmp"))
  testthat::expect_error(ev_trawl_fit(data, depth, "GMM", bounds = "tmp"))
})

test_that("ev_trawl_fit_does_not_smoke", {
  n <- 2000
  depth <- 5
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]

  fit_params <- ev_trawl_fit(data, depth, "GMM")
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:3] > 0))

  fit_params <- ev_trawl_fit(data, depth, "PL")
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:3] > 0))

  fit_params <- ev_trawl_fit(data, depth, "GMM", mode = "full")
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:3] > 0))

  fit_params <- ev_trawl_fit(data, depth, "PL", mode = "full")
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:3] > 0))
})


test_that("ev_trawl_fit_does_not_smoke", {
  n <- 2000
  depth <- 5
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]

  fit_params <- fit(data, depth, "PL", parallel = T)
  testthat::expect_equal(length(fit_params), 4)
  testthat::expect_true(all(fit_params[2:3] > 0))
})
