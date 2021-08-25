test_that("Custom likelihood - positive xi, only exceedance", {
  xi <- 1.
  kappa <- 9.
  sigma <- 1.
  n <- 30000

  set.seed(42)
  p_zero <- 1 - 1 / (1 + kappa)
  zeroes <- runif(n = n, min = 0, max = 1)
  test_samples <- evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma)
  test_samples[which(zeroes < p_zero)] <- 0.0

  cm_mle <- custom_likelihood(data = test_samples)
  pars <- c(xi, sigma, kappa)
  expect_true(abs(cm_mle(pars)) < 1e10)
})

test_that("Custom likelihood - positive xi", {
  xi <- .1
  kappa <- 9.
  sigma <- 1.
  n <- 30000

  set.seed(42)
  test_samples <- evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma)
  cm_mle <- custom_likelihood(data = test_samples)
  pars <- c(xi, sigma, kappa)
  testthat::expect_true(abs(cm_mle(pars)) < 1e10)
  testthat::expect_false(any(is.na(cm_mle(pars))))
})

test_that("Custom likelihood - matrix error", {
  xi <- .1
  kappa <- 9.
  sigma <- 1.
  n <- 30000

  set.seed(42)
  test_samples <- evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma)
  matrix_test_samples <- cbind(test_samples, test_samples)
  testthat::expect_error(custom_likelihood(data = matrix_test_samples))
})


test_that("Custom MLE - positive xi", {
  xi <- 1.
  kappa <- 9.
  sigma <- 1.
  n <- 30000
  p_zero <- 1 - 1. / (1. + kappa)

  set.seed(42)
  zeroes <- runif(n = n, min = 0, max = 1)
  test_samples <- evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma)
  test_samples[which(zeroes < p_zero)] <- 0.0
  cm_mle <- custom_marginal_mle(data = test_samples)
  testthat::expect_equal(cm_mle, c(xi, sigma, kappa), tolerance = 5e-2)
})


test_that("Composite MLE - score ACF", {
  xi <- 1.
  kappa <- 9.
  sigma <- 1.
  n <- 30000
  p_zero <- 1 - 1. / (1. + kappa)

  set.seed(42)
  zeroes <- runif(n = n, min = 0, max = 1)
  test_samples <- evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma)
  test_samples[which(zeroes < p_zero)] <- 0.0

  k <- 3
  cm_hac <- composite_marginal_hac(
    data = test_samples, params = c(xi, sigma, kappa),
    k = k, 1000
  )
  testthat::expect_equal(dim(cm_hac), c(k, k))
  testthat::expect_true(all(diag(cm_hac) > 0))
})



test_that("Composite MLE - positive/negative xi", {
  kappa <- 4.
  p_zero <- 1 - 1. / (1. + kappa)
  n <- 50000
  xi_seq <- seq(from = -.3, to = .3, length.out = 4)
  sigma_seq <- seq(from = .1, to = 10, length.out = 4)
  for (xi in xi_seq) {
    for (sigma in sigma_seq) {
      set.seed(42)
      zeroes <- runif(n = n, min = 0, max = 1)
      test_samples <- evir::rgpd(n = n, mu = 0, xi = xi, beta = sigma)
      test_samples[which(zeroes < p_zero)] <- 0.0
      cm_mle <- composite_marginal_mle(data = test_samples)
      testthat::expect_equal(
        cm_mle / c(xi, sigma, kappa), rep(1, 3),
        tolerance = .50
      )
    }
  }
})
