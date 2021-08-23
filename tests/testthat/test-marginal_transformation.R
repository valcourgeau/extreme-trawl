test_that("MapInverse into Map is identity", {
  set.seed(42)
  gpd_sample <- evir::rgpd(n = 10000, mu = 0.0, beta = 5, xi = -0.2)
  tmp <- transformation_map(
    transformation_map_inverse(gpd_sample, params = c(-0.2, 5, 3)),
    params = c(-0.2, 5, 3)
  )
  expect_true(all(abs(sort(tmp) - sort(gpd_sample)) < 1e-10))
})


test_that("Map into MapInverse is identity", {
  set.seed(42)
  gpd_sample <- evir::rgpd(n = 10000, mu = 0.0, beta = 5, xi = -0.2)
  tmp <- transformation_map_inverse(
    transformation_map(gpd_sample, params = c(-0.2, 5, 3)),
    params = c(-0.2, 5, 3)
  )
  expect_true(all(abs(sort(tmp) - sort(gpd_sample)) < 1e-10))
})


test_that("MapInverse trf", {
  set.seed(42)
  xi <- -.2
  sigma <- 10.0
  kappa <- 10
  n <- 100000
  gpd_sample <- evir::rgpd(n = n, mu = 0.0, beta = sigma, xi = xi)
  tmp <- transformation_map_inverse(
    gpd_sample,
    params = c(xi, sigma, kappa)
  )

  pars <- as.numeric(evir::gpd(tmp, threshold = 0, method = "ml")$par.ests)
  testthat::expect_equal(pars, c(1, 1 + kappa), tolerance = 5e-2)
})

test_that("Map trf", {
  set.seed(42)
  xi <- -.2
  sigma <- 10.0
  kappa <- 10
  n <- 100000
  gpd_sample <- evir::rgpd(n = n, mu = 0.0, beta = 1 + kappa, xi = 1.)
  tmp <- transformation_map(
    gpd_sample,
    params = c(xi, sigma, kappa)
  )
  pars <- as.numeric(evir::gpd(tmp, threshold = 0, method = "ml")$par.ests)
  testthat::expect_equal(pars, c(xi, sigma), tolerance = 5e-2)
})
