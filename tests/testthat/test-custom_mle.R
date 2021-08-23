test_that("Custom likelihood - positive xi, only exceedance", {
  xi <- 1.
  kappa <- 9.
  sigma <- 1.
  n <- 100000

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
  n <- 100000

  set.seed(42)
  test_samples <- evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma)
  cm_mle <- custom_likelihood(data = test_samples)
  pars <- c(xi, sigma, kappa)
  expect_true(abs(cm_mle(pars)) < 1e10)
})


test_that("Custom MLE - positive xi", {
  xi <- 1.
  kappa <- 9.
  sigma <- 1.
  n <- 100000
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
  n <- 100000
  p_zero <- 1 - 1. / (1. + kappa)

  set.seed(42)
  zeroes <- runif(n = n, min = 0, max = 1)
  test_samples <- evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma)
  test_samples[which(zeroes < p_zero)] <- 0.0

  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 20000

  print(composite_marginal_hac(
    data = pollution_data[, test_column], params = c(xi, sigma, kappa),
    k = 3, 1000
  ))
})



test_that("Composite MLE - positive/negative xi", {
  kappa <- 9.
  p_zero <- 1 - 1. / (1. + kappa)
  n <- 100000
  xi_seq <- seq(from = -.3, to = .3, length.out = 4)
  sigma_seq <- seq(from = .1, to = 10, length.out = 4)
  for (xi in xi_seq) {
    for (sigma in sigma_seq) {
      set.seed(42)
      zeroes <- runif(n = n, min = 0, max = 1)
      test_samples <- evir::rgpd(n = n, mu = 0, xi = xi, beta = sigma)
      test_samples[which(zeroes < p_zero)] <- 0.0
      cm_mle <- composite_marginal_mle(data = test_samples)
      print(cm_mle)
      testthat::expect_equal(
        cm_mle / c(xi, sigma, kappa),
        rep(1, 3),
        tolerance = .30
      )
      # 30% tolerance
    }
  }
})
