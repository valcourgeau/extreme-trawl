test_that("trawl_gamma_error", {
  params <- c(0.8, 3)
  h_max <- 50
  h_collection <- seq_len(h_max)
  one_vals <- vapply(h_collection, function(h) {
    gamma_trawl$trawl_b_one(params, h)
  }, 1)
  testthat::expect_error(gamma_trawl$trawl_b_one(-params, h))
  testthat::expect_error(gamma_trawl$trawl_b_one(-params, h))
  testthat::expect_error(gamma_trawl$trawl_b_two(-params, h))
  testthat::expect_error(gamma_trawl$trawl_b_three(-params, h))
})

test_that("trawl_gamma_one", {
  params <- c(0.8, 3)
  h_max <- 50
  h_collection <- seq_len(h_max)
  one_vals <- vapply(h_collection, function(h) {
    gamma_trawl$trawl_b_one(params, h)
  }, 1)
  testthat::expect_equal(length(one_vals), h_max)
  testthat::expect_true(all(one_vals > 0))
})

test_that("trawl_gamma_two", {
  params <- c(0.8, 3)
  h_max <- 50
  h_collection <- seq_len(h_max)
  two_vals <- vapply(h_collection, function(h) {
    gamma_trawl$trawl_b_two(params, h)
  }, 1)
  testthat::expect_equal(length(two_vals), h_max)
  testthat::expect_true(all(two_vals > 0))
})

test_that("trawl_gamma_three", {
  params <- c(0.8, 3)
  h_max <- 50
  h_collection <- seq_len(h_max)
  three_vals <- vapply(h_collection, function(h) {
    gamma_trawl$trawl_b_three(params, h)
  }, 1)
  one_vals <- vapply(h_collection, function(h) {
    gamma_trawl$trawl_b_one(params, h)
  }, 1)
  testthat::expect_equal(length(three_vals), h_max)
  testthat::expect_true(all(three_vals > 0))
  testthat::expect_equal(one_vals, three_vals, tol = 1e-4)
})

test_that("trawl_gamma_completement", {
  params <- c(0.8, 3)
  h_max <- 50
  h_collection <- seq_len(h_max)
  should_equal_to_one <- vapply(h_collection, function(h) {
    gamma_trawl$trawl_b_one(params, h) +
      gamma_trawl$trawl_b_two(params, h)
  }, 1)
  testthat::expect_equal(length(should_equal_to_one), h_max)
  testthat::expect_equal(
    should_equal_to_one, rep(params[1] / (params[2] - 1), h_max),
    tolerance = 1e-3
  )
})

test_that("trawl_gamma_cfg", {
  cfg <- gamma_trawl$config()
  testthat::expect_equal(cfg$n_params, 2)
  testthat::expect_equal(length(cfg$lower), cfg$n_params)
  testthat::expect_equal(length(cfg$upper), cfg$n_params)
  testthat::expect_true(all(cfg$lower > 0))
  testthat::expect_true(all(cfg$upper < Inf))
})
