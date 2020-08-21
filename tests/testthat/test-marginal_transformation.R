test_that("MapInverse into Map is identity", {
  set.seed(42)
  gpd_sample <- evir::rgpd(n = 10000, mu = 0.0, beta = 5, xi = -0.2)
  tmp <- TransformationMap(
    TransformationMapInverse(gpd_sample, params = c(-0.2, 5, 3)),
    params = c(-0.2, 5, 3))
  expect_true(all(abs(sort(tmp)-sort(gpd_sample)) < 1e-10))
})

test_that("Map into MapInverse is identity", {
  set.seed(42)
  gpd_sample <- evir::rgpd(n = 10000, mu = 0.0, beta = 5, xi = -0.2)
  tmp <- TransformationMapInverse(
    TransformationMap(gpd_sample, params = c(-0.2, 5, 3)),
    params = c(-0.2, 5, 3))
  expect_true(all(abs(sort(tmp)-sort(gpd_sample)) < 1e-10))
})
