test_that("trawl acf", {
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  a <- acf_trawl_single(
    h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
  )
  b <- cpp_acf_trawl_single(
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

  a <- acf_trawl(
    h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
  )
  b <- cpp_acf_trawl(
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
    acf_trawl_single(
      h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
    ),
    times = trials
  )$time / time_divisor
  time_new <- microbenchmark::microbenchmark(
    cpp_acf_trawl_single(
      h = h, alpha = alpha, beta = beta, kappa = kappa, rho = rho
    ),
    times = trials
  )$time / time_divisor
  # 10 times as fast
  testthat::expect_lte(mean(time_new) / mean(time_old), 1)
})


test_that("acf_trawl__time_trial", {
  time_divisor <- 1000000
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  h_collection <- 1:20

  time_old <- microbenchmark::microbenchmark(
    acf_trawl(
      h = h_collection, alpha = alpha, beta = beta, kappa = kappa, rho = rho
    ),
    times = 5
  )$time / time_divisor
  time_new <- microbenchmark::microbenchmark(
    cpp_acf_trawl(
      h = h_collection, alpha = alpha, beta = beta, kappa = kappa, rho = rho
    ),
    times = 5
  )$time / time_divisor

  testthat::expect_equal(
    time_new / time_old, rep(0., length(time_old)),
    tolerance = .15
  )
})

test_that("acf_trawl__vals", {
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  h_collection <- 1:15

  r_only_acf <- acf_trawl(
    h = h_collection, alpha = alpha, beta = beta,
    kappa = kappa, rho = rho, cov = T
  )

  testthat::expect_false(any(is.na(r_only_acf > 0)))
  testthat::expect_true(all(r_only_acf > 0))
  testthat::expect_true(all(diff(r_only_acf) < 0))
})


test_that("acf_trawl_revisited_num_approx__vals", {
  h_max <- 5
  alpha <- 1
  beta <- 1
  kappa <- 19.
  rho <- .2

  h_collection <- seq_len(h_max)

  cpp_acf <- cpp_acf_trawl(
    h = h_collection, alpha = alpha, beta = beta,
    kappa = kappa, rho = rho, cov = T
  )
  acf_trawl_vals <- acf_trawl_revisited_num_approx(
    h_collection,
    alpha = alpha, beta = beta, kappa = kappa, rho = rho
  )
  testthat::expect_equal(unlist(acf_trawl_vals), unlist(cpp_acf), tol = 1e-3)
})

test_that("cross_moment__vals", {
  h <- 1
  alpha <- 1
  beta <- 1
  kappa <- 19.
  rho <- .2
  type <- "exp"

  h_collection <- 1:15

  cm_trawls <- cross_moment_trawls(
    h = h, alpha = alpha, beta = beta, rho = rho, kappa = kappa, type = type
  )
  end_seq <- 50
  delta <- 0.1
  trawl_fct <- get_trawl_functions(type = type)
  b_1_func <- trawl_fct[[1]]
  b_2_func <- trawl_fct[[2]]
  b_3_func <- trawl_fct[[3]]
  a_total <- b_1_func(param = rho, h = h) + b_2_func(param = rho, h = h)
  b_0_minus_h <- -alpha * b_3_func(param = rho, h = h) / a_total
  b_0_h <- -alpha * b_2_func(param = rho, h = h) / a_total

  seq_kappa <- seq(kappa, kappa + end_seq, by = delta)
  cpp_cm_trawls <- cross_moment(
    xs = seq_kappa, beta = beta, delta = delta,
    b_oh = b_0_h, b_o_exc_h = b_0_minus_h
  )

  testthat::expect_equal(cpp_cm_trawls, cm_trawls, tol = 1e-3)
})
