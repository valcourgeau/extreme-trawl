
test_that("trawl_exponential_error", {
  rho <- 0.2
  h_max <- 50
  h_collection <- seq_len(h_max)
  one_vals <- vapply(h_collection, function(h) {
    exponential_trawl$trawl_b_one(rho, h)
  }, 1)
  testthat::expect_error(exponential_trawl$trawl_b_one(c(rho, rho), h))
  testthat::expect_error(exponential_trawl$trawl_b_two(c(rho, rho), h))
  testthat::expect_error(exponential_trawl$trawl_b_three(c(rho, rho), h))
})

test_that("trawl_exponential_one", {
  rho <- 0.2
  h_max <- 50
  h_collection <- seq_len(h_max)
  one_vals <- vapply(h_collection, function(h) {
    exponential_trawl$trawl_b_one(rho, h)
  }, 1)
  testthat::expect_equal(length(one_vals), h_max)
  testthat::expect_true(all(one_vals > 0))
})

test_that("trawl_exponential_two", {
  rho <- 0.2
  h_max <- 50
  h_collection <- seq_len(h_max)
  two_vals <- vapply(h_collection, function(h) {
    exponential_trawl$trawl_b_two(rho, h)
  }, 1)
  testthat::expect_equal(length(two_vals), h_max)
  testthat::expect_true(all(two_vals > 0))
})

test_that("trawl_exponential_three", {
  rho <- 0.2
  h_max <- 50
  h_collection <- seq_len(h_max)
  three_vals <- vapply(h_collection, function(h) {
    exponential_trawl$trawl_b_three(rho, h)
  }, 1)
  one_vals <- vapply(h_collection, function(h) {
    exponential_trawl$trawl_b_one(rho, h)
  }, 1)
  testthat::expect_equal(length(three_vals), h_max)
  testthat::expect_true(all(three_vals > 0))
  testthat::expect_equal(one_vals, three_vals, tol = 1e-4)
})


test_that("trawl_exponential_completement", {
  rho <- 0.2
  h_max <- 50
  h_collection <- seq_len(h_max)
  should_equal_to_one <- vapply(h_collection, function(h) {
    exponential_trawl$trawl_b_one(rho, h) +
      exponential_trawl$trawl_b_two(rho, h)
  }, 1)
  testthat::expect_equal(length(should_equal_to_one), h_max)
  testthat::expect_equal(
    should_equal_to_one, rep(1 / rho, h_max),
    tolerance = 1e-3
  )
})

test_that("trawl_exponential_cfg", {
  cfg <- exponential_trawl$config()
  testthat::expect_equal(cfg$n_params, 1)
  testthat::expect_equal(length(cfg$lower), cfg$n_params)
  testthat::expect_equal(length(cfg$upper), cfg$n_params)
  testthat::expect_true(all(cfg$lower > 0))
  testthat::expect_true(all(cfg$upper < Inf))
})
