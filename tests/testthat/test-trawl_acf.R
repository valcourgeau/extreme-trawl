test_that("trawl acf", {
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  a <- trawl_autocorrelation$acf_trawl_single(
    h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
  )
  b <- trawl_autocorrelation$cpp_acf_trawl_single(
    h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
  )

  testthat::expect_equal(a, b, tolerance = 1e-3)
})

test_that("acf trawl acf", {
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  a <- trawl_autocorrelation$acf_trawl(
    h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
  )
  b <- trawl_autocorrelation$cpp_acf_trawl(
    h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
  )

  testthat::expect_equal(a, b, tolerance = 1e-3)
})

test_that("acf_trawl_single__time_trial", {
  time_divisor <- 1000000
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  trials <- 50

  time_old <- microbenchmark::microbenchmark(
    trawl_autocorrelation$acf_trawl_single(
      h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
    ),
    times = trials
  )$time / time_divisor
  time_new <- microbenchmark::microbenchmark(
    trawl_autocorrelation$cpp_acf_trawl_single(
      h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
    ),
    times = trials
  )$time / time_divisor
  # 10 times as fast
  testthat::expect_equal(
    mean(time_new) / mean(time_old), .10,
    tolerance = .05
  )
})



test_that("trawl_autocorrelation$acf_trawl__time_trial", {
  time_divisor <- 1000000
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  h_collection <- 1:20

  time_old <- microbenchmark::microbenchmark(
    trawl_autocorrelation$acf_trawl(
      h = h_collection, alpha = alpha, beta = beta, kappa = kappa, rho = rho
    ),
    times = 5
  )$time / time_divisor
  time_new <- microbenchmark::microbenchmark(
    trawl_autocorrelation$cpp_acf_trawl(
      h = h_collection, alpha = alpha, beta = beta, kappa = kappa, rho = rho
    ),
    times = 5
  )$time / time_divisor

  testthat::expect_equal(
    time_new / time_old, rep(0., length(time_old)),
    tolerance = .15
  )
})


test_that("trawl_autocorrelation$acf_trawl__vals", {
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  h_collection <- 1:15

  r_only_acf <- trawl_autocorrelation$acf_trawl(
    h = h_collection, alpha = alpha, beta = beta,
    kappa = kappa, rho = rho, cov = T
  )

  testthat::expect_false(any(is.na(r_only_acf > 0)))
  testthat::expect_true(all(r_only_acf > 0))
  testthat::expect_true(all(diff(r_only_acf) < 0))
})

test_that("trawl_autocorrelation$acf_trawl__vals", {
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  h_collection <- 1:15

  cpp_acf <- trawl_autocorrelation$acf_trawl(
    h = h_collection, alpha = alpha, beta = beta,
    kappa = kappa, rho = rho, cov = T
  )

  testthat::expect_false(any(is.na(cpp_acf > 0)))
  testthat::expect_true(all(cpp_acf > 0))
  testthat::expect_true(all(diff(cpp_acf) < 0))
})
