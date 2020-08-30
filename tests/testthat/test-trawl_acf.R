test_that("trawl acf", {
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2
  # print(acf_trawl_inv(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho))
  a <- acf_trawl(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho)
  b <- AcfTrawl(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho)

  testthat::expect_equal(a, b, tolerance=1e-3)
})

test_that("trawl acf time", {
  TIME_DIVISOR <- 1000000
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  time_old <- microbenchmark::microbenchmark(acf_trawl(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho), times=5)$time / TIME_DIVISOR
  time_new <- microbenchmark::microbenchmark(AcfTrawl(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho), times=5)$time / TIME_DIVISOR

  testthat::expect_equal(time_new/time_old, rep(0., length(time_old)), tolerance=.15)
})

test_that("collection trawl acf", {
  TIME_DIVISOR <- 1000000
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2

  h_collection <- 1:20

  time_old <- microbenchmark::microbenchmark(
    acf_trawl_num_approx(h=h_collection, alpha=alpha, beta=beta, kappa=kappa, rho=rho),
    times=5)$time / TIME_DIVISOR
  time_new <- microbenchmark::microbenchmark(
    AcfTrawlCollection(h=h_collection, alpha=alpha, beta=beta, kappa=kappa, rho=rho),
    times=5)$time / TIME_DIVISOR

  testthat::expect_equal(time_new/time_old, rep(0., length(time_old)), tolerance=.15)
})



