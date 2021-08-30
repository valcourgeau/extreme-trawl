test_that("autocovariance_matrix__values", {
  d <- 5
  n <- 1000
  vals <- matrix(rnorm(n * d), ncol = d)
  k <- 10
  autocov_vals <- autocovariance_matrix(vals, k)
  testthat::expect_equal(dim(autocov_vals), c(k + 1, d, d))
})

test_that("autocovariance_matrix__errors", {
  testthat::expect_error(autocovariance_matrix(vals, 0))
  testthat::expect_error(autocovariance_matrix(vals, -1))
})

test_that("make_hac__matrix", {
  set.seed(43)
  d <- 5
  n <- 10000
  vals <- matrix(rnorm(n * d), ncol = d)
  k <- 10
  autocov_vals <- autocovariance_matrix(vals, k)
  m_hac <- make_hac(autocov_vals)
  testthat::expect_equal(dim(m_hac), c(d, d))
  testthat::expect_equal(m_hac, diag(d), tolerance = 0.05)

  m_hac <- make_hac(autocov_vals, near_pd = T)
  testthat::expect_equal(dim(m_hac), c(d, d))
  testthat::expect_equal(as.vector(m_hac), as.vector(diag(d)), tolerance = 0.05)
})
