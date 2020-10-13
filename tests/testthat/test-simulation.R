test_that("ExceedancesSimulation - simple", {
  set.seed(42)
  params <- c(.1, 1., 19, .20)
  n <- 10000
  vd <- 50
  type <- 'exp'

  exc <- ExceedancesSimulation(params=params, n=n, vanishing_depth=vd, type=type)
  gpd_fit_sim <- evir::gpd(exc$exceedances, threshold = 0.0, method = 'ml') # returns (exceedances, latent)

  testthat::expect_equal(unname(gpd_fit_sim$par.ests)[1], params[1], tolerance=1.96*unname(gpd_fit_sim$par.ses)[1])
  testthat::expect_equal(unname(gpd_fit_sim$par.ests)[2], params[2], tolerance=1.96*unname(gpd_fit_sim$par.ses)[2])

  testthat::expect_equal(mean(exc$exceedances > 0), .05, tolerance=8e-3)
  print(as.vector(acf(exc$exceedances, lag.max = 10, plot = F)$acf))
  print(TrawlAutocorrelation$AcfTrawlCollection(
    h=c(.01, 1:10), alpha=1, beta=1, kappa=params[3], rho=params[4],
    cov = F, type=type, delta = .1, end_seq=100))
})
