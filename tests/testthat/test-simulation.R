test_that("ExceedancesSimulation - simple", {
  set.seed(42)
  params <- c(.1, 1., 19, .2)
  n <- 20000
  vd <- 50
  type <- 'exp'

  exc <- ExceedancesSimulation(params=params, n=n, vanishing_depth=vd, type=type)
  print(evir::gpd(exc, threshold = 0.0, method = 'ml'))
  testthat::expect_equal(mean(exc > 0), .05, tolerance=8e-3)
  print(as.vector(acf(exc, lag.max = 40, plot = F)$acf))
  print(TrawlAutocorrelation$AcfTrawlCollection(h=c(.01, 1:10), alpha=1, beta=1, kappa=params[3], rho=params[4], cov = F, type=type))

})
