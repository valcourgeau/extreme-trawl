test_that("trawl acf", {
  h <- 1
  alpha <- 2.
  beta <- 10.
  kappa <- 19.
  rho <- .2
  # print(acf_trawl_inv(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho))
  a <- acf_trawl(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho)
  b <- acf_trawl_revised(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho)
  # print(microbenchmark::microbenchmark(acf_trawl(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho), times=5))
  # print(microbenchmark::microbenchmark(acf_trawl_revised(h=h, alpha=alpha, beta=beta, kappa=kappa, rho=rho), times=5))

  testthat::expect_equal(a, b, tolerance=1e-3)
  # print(FirstMoment(c(1,2), .1, 10., -.3, -1.4))
  # print(CrossMoment(c(1,2), .1, 10., -.3, -1.4))
})
